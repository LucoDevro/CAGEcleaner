#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.local_run import LocalRun
from cagecleaner.region_run import RegionRun
from cagecleaner.file_utils import _extractOneRegion

import logging
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.contrib.concurrent import thread_map

LOG = logging.getLogger()


class LocalRegionRun(LocalRun, RegionRun):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        return None
    
    def extractRegions(self):
        """
        This method extracts the genomic regions surrounding each cluster hit using the specified sequence margin.
        Regions at contig edges are discarded when applying the strict region flag.
        """
        
        # Extract regions in parallel
        regions = [r.to_dict() for i,r in self.binary_df.iterrows()]
        contig_ends = thread_map(lambda x: _extractOneRegion(x, self.margin, self.TEMP_GENOME_DIR, 
                                                             self.DEREP_IN_DIR, self.strict_regions), 
                                 regions, 
                                 max_workers = self.cores,
                                 leave = False,
                                 disable = self.no_progress)
        contig_end = sum(contig_ends)

        LOG.info(f'{contig_end} regions were at a contig end.')
        if self.strict_regions:
            LOG.info('These regions have been discarded from the analysis.')
            
        return None
    
    def mapDereplicationToBinary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).  
        
        Mutates:
            self.binary_df: pd.DataFrame: The binary table derived from a cblaster Session object.
        """
        
        # Parse the MMseqs clustering table
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / "derep_cluster.tsv"
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 names = ['representative', 'Scaffold'],
                                 header = None, index_col = 'Scaffold'
                                 )
        # Determine dereplication status
        derep_df['dereplication_status'] = derep_df.index == derep_df['representative']
        derep_df['dereplication_status'] = np.where(derep_df['dereplication_status'], 'dereplication_representative', 'redundant')
        # Add dereplication status columns to binary table by joining with the clustering table
        LOG.debug("Joining MMseqs2 clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.merge(derep_df, left_on = "Scaffold", right_index = True)
        
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])           
        LOG.info("Mapping done!")
        
        return None
    
    def run(self):
        
        LOG.info("--- STEP 1: Staging genomes for dereplication. ---")
        self.prepareGenomes()
        
        LOG.info("--- STEP 2: Extracting genomic regions. ---")
        self.extractRegions()
        
        LOG.info("--- STEP 3: Dereplicating. ---")
        self.dereplicateRegions()
        
        LOG.info("--- STEP 4: Mapping dereplication output to binary table. ---")
        self.mapDereplicationToBinary()
        
        LOG.info("--- STEP 5: Recovering hit diversity. ---")
        self.recoverHits()
        
        LOG.info("--- STEP 6: Filtering session file. ---")
        self.filterSession()
        
        LOG.info("--- STEP 7: Generating output files")
        self.generateOutput()
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
    
        return None
    
    