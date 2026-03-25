#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.local_run import LocalRun
from cagecleaner.genome_run import GenomeRun
from cagecleaner import util

import logging
import os
import pandas as pd
from pathlib import Path


LOG = logging.getLogger()


class LocalGenomeRun(LocalRun, GenomeRun):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        return None
    
    def mapDereplicationToBinary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).  
        
        Mutates:
            self.binary_df: pd.DataFrame: The binary table derived from a cblaster Session object.
        """
        # Full genome dereplication using skDER
        def extractAssembly(file_path: str) -> str:
            return util.removeSuffixes(os.path.basename(file_path))
        
        def renameLabel(label: str) -> str:
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant',
                       'outside_cutoffs_requested': 'redundant'} # edge case that clusters by skani dist, but fails the clustering cutoffs
            return mapping[label]
        
        LOG.debug("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 converters = {'assembly': extractAssembly,
                                               'representative': extractAssembly,
                                               'dereplication_status': renameLabel},
                                 names = ['assembly', 'representative', 'dereplication_status'],
                                 usecols = [0,1,4], header = 0, index_col = 'assembly'
                                 )
        # Join with binary df on Organism column. 
        # Every Organism row is retained (left join).
        # If there is a match between binary_df['Organism'] and derep_df['assembly'] (index), the representative and status is added.
        LOG.debug("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(derep_df, on='Organism')
        
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])           
        LOG.info("Mapping done!")
        
        return None
    
    def run(self):
        
        LOG.info("--- STEP 1: Staging genomes for dereplication. ---")
        self.prepareGenomes()
        
        LOG.info("--- STEP 2: Dereplicating. ---")
        self.dereplicateGenomes()
        
        LOG.info("--- STEP 3: Mapping dereplication output to binary table. ---")
        self.mapDereplicationToBinary()
        
        LOG.info("--- STEP 4: Recovering hit diversity. ---")
        self.recoverHits()
        
        LOG.info("--- STEP 5: Filtering session file. ---")
        self.filterSession()
        
        LOG.info("--- STEP 6: Generating output files")
        self.generateOutput()
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
    