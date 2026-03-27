#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.remote_run import RemoteRun
from cagecleaner.region_run import RegionRun
from cagecleaner import util

import logging
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.contrib.concurrent import thread_map

LOG = logging.getLogger()


class RemoteRegionRun(RemoteRun, RegionRun):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        return None
    
    def downloadRegions(self) -> None:
        accessions = self.binary_df['Scaffold'].to_list()
        ranges = (self.binary_df['Start'].map(str) + ':' + self.binary_df['End'].map(str)).to_list()
        regions = list(zip(accessions, ranges))
        
        if self.download_workers > 2:
            max_workers = 2
            LOG.warning("You are requesting too many download workers by NCBI policy. Limiting the number to 2 for compliance.")
        else:
            max_workers = self.download_workers
        thread_map(lambda x: util._downloadOneRegion(x, self.DEREP_IN_DIR), regions, 
                   max_workers = max_workers,
                   leave = False,
                   disable = self.no_progress)
        
        return None
    
    def mapDereplicationToBinary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
        """
        # Read the MMseqs2 clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / "derep_cluster.tsv"
        # Convert to dataframe
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 names = ['representative', 'Scaffold'],
                                 header = None, index_col = 'Scaffold'
                                 )
        # Add dereplication status column based on whether the scaffold's ID is the same as the representative's
        # If there is a match between binary_df['assembly_file'] and derep_dfs['assembly'] (its index column), the representative and status is added.
        derep_df['dereplication_status'] = derep_df.index == derep_df['representative']
        derep_df['dereplication_status'] = np.where(derep_df['dereplication_status'], 'dereplication_representative', 'redundant')
        # Discard region coordinates from labels for proper merging later on
        derep_df.index = derep_df.index.str.split(':').str[0]
        derep_df['representative'] = derep_df['representative'].str.split(':').str[0]
        
        LOG.debug("Joining MMseqs2 clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.merge(derep_df, left_on = "Scaffold", right_index = True)
        
        # Sort by representative ID and then by dereplication status
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])
        LOG.info("Mapping done!")
            
        return None
    
    def run(self):
        
        LOG.info('--- STEP 1: Downloading genomic regions. ---')
        self.downloadRegions() # Download genomic regions using ncbi_acc_download
        
        LOG.info("--- STEP 2: Dereplicating ---")
        self.dereplicateRegions()  # Dereplicate.
        
        LOG.info("--- STEP 3: Mapping dereplication clustering to binary table ---")
        self.mapDereplicationToBinary()  # Map each row in the binary table with its representative
        
        LOG.info("--- STEP 4: Recovering hit diversity ---")
        self.recoverHits()  # Recover hits by content and score depending on user input
        
        LOG.info("--- STEP 5: Filtering original session file. ---")
        self.filterSession()  # Filter the original session file, retaining only the dereplicated hits.
        
        LOG.info("--- STEP 6: Generating output files ---")
        self.generateOutput()  # Generate all output files.
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
