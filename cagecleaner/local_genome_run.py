#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.local_run import LocalRun
from cagecleaner.genome_run import GenomeRun
from cagecleaner.file_utils import remove_suffixes

import logging
import pandas as pd
from pathlib import Path


LOG = logging.getLogger()


class LocalGenomeRun(LocalRun, GenomeRun):
    """
    Subclass orchestrating the workflow for dereplication by whole-genome similarity using genomes from local sources.
    
    This class combines local genome file handling with whole-genome dereplication workflows.
    It stages the local genome assemblies for dereplication (converting genbanks to fastas), performs
    whole-genome dereplication using skDER, and integrates the dereplication results back into the binary table.
    Unmapped scaffolds (scaffolds for which no assembly file was found) are tracked and reported separately.
    
    Inherits from:
        LocalRun: Intermediary class providing local file handling utilities.
        GenomeRun: Intermediary class providing genome dereplication utilities.
    """
    
    def __init__(self, args):
        """
        Initialise a LocalGenomeRun instance.
        
        Calls the parent class constructor to set up common configuration, temporary directories,
        and logging infrastructure.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
        
        Returns:
            None
        """
        
        super().__init__(args)
        
        return None
    
    def join_dereplication_with_binary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).  
        
        Mutates:
            self.binary_df: pd.DataFrame: The binary table derived from a cblaster Session object.
            
        Join dereplication clustering results with the binary table.
        
        Reads the skDER clustering output file, converts it to a DataFrame, and joins it with
        the binary table based on assembly ID. This associates each genome in the binary table
        with its dereplication status (representative or redundant) and representative assembly.
        The resulting table is sorted by representative and dereplication status for clarity.
        
        Mutates:
            self.binary_df (pd.DataFrame): Updated in-place with additional columns for
                'representative' and 'dereplication_status', and sorted by these columns.
        
        Returns:
            None
        
        Raises:
            FileNotFoundError: If the skDER clustering file cannot be found at the expected path.
            RuntimeError: If the dereplication table is empty.
            RuntimeError: If the binary table is empty after joining with the dereplication table.
            
        Notes:
            This is the workflow-specific implementation of the abstract method inherited from its grandparent class Run.
        """
        # Full genome dereplication using skDER
        def extract_assembly_id(file_path: str) -> str: return remove_suffixes(Path(file_path).name)
        
        def rename_label(label: str) -> str:
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant',
                       'outside_cutoffs_requested': 'redundant'} # edge case that clusters by skani dist, but fails the clustering cutoffs
            return mapping[label]
        
        LOG.debug("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        try:
            derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                     converters = {'assembly': extract_assembly_id,
                                                   'representative': extract_assembly_id,
                                                   'dereplication_status': rename_label},
                                     names = ['assembly', 'representative', 'dereplication_status'],
                                     usecols = [0,1,4], header = 0, index_col = 'assembly'
                                     )
        except FileNotFoundError as err:
            LOG.critical(f'{err}')
            raise err
            
        if derep_df.empty:
            msg = "Dereplication table is empty."
            LOG.error(msg)
            raise RuntimeError(msg)
            
        # Join with binary df on Organism column. 
        # Every Organism row is retained (left join).
        # If there is a match between binary_df['Organism'] and derep_df['assembly'] (index), the representative and status is added.
        LOG.debug("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(derep_df, on='Organism')
        
        # Verify binary table is not empty
        if self.binary_df.empty:
            msg = "Binary table is empty after joining with the dereplication table!"
            LOG.error(msg)
            raise RuntimeError(msg)
        
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])           
        LOG.info("Mapping done!")
        
        return None
    
    def run(self):
        """
        Execute the complete local genome dereplication pipeline.
        
        Orchestrates all processing steps in sequence: stages local genome assemblies for dereplication,
        runs skDER dereplication on the staged genomes, joins dereplication clustering results with the binary table,
        recovers hit diversity information from the dereplication output, filters the original session
        based on dereplication results, and generates final output files with dereplication metadata.
        Cleans up temporary working directories upon completion.
        
        Returns:
            None
        """
        
        LOG.info("--- STEP 1: Staging genomes for dereplication. ---")
        self.prepare_genomes()
        
        LOG.info("--- STEP 2: Dereplicating. ---")
        self.dereplicate_genomes()
        
        LOG.info("--- STEP 3: Mapping dereplication output to binary table. ---")
        self.join_dereplication_with_binary()
        
        LOG.info("--- STEP 4: Recovering hit diversity. ---")
        self.recover_hits()
        
        LOG.info("--- STEP 5: Filtering session file. ---")
        self.filter_session()
        
        LOG.info("--- STEP 6: Generating output files")
        self.generate_output()
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
    