#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.local_run import LocalRun
from cagecleaner.region_run import RegionRun
from cagecleaner.file_utils import _extract_one_region

import logging
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.contrib.concurrent import thread_map

LOG = logging.getLogger()


class LocalRegionRun(LocalRun, RegionRun):
    """
    Subclass orchestrating the workflow for dereplication by region similarity using regions
    extracted from local sources.
    
    This class combines local file handling with region-based dereplication workflows.
    It extracts genomic regions of interest (with optional sequence margins) from local
    assembly files, performs MMseqs2-based sequence dereplication on the extracted regions,
    and integrates dereplication results back into the binary table. Handles contig edge
    cases where regions with sequence margins extend beyond scaffold boundaries according
    to user-specified behavior (keep but clip them (permissive), or discard them (strict)).
    
    Inherits from:
        LocalRun: Intermediary class providing local file handling utilities.
        RegionRun: Intermediary class providing region dereplication utilities.
    """
    def __init__(self, args):
        """
        Initialise a LocalRegionRun instance.
        
        Calls the parent class constructor to set up common configuration, temporary directories,
        and logging infrastructure.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
        
        Returns:
            None
        """
        
        super().__init__(args)
        
        return None
    
    
    def extract_regions(self):
        """
        Extract genomic regions surrounding cluster hits using sequence margins.
        
        Processes each cluster hit in the binary table to extract the genomic region
        with sequence margins from the local assembly files. Extraction is performed
        in parallel using multiple worker threads. Regions that extend beyond contig
        boundaries are treated as specified by the user (strict_regions flag).
        
        When strict_regions is enabled, regions at contig edges are excluded from downstream
        dereplication analysis. When disabled (permissive mode), such regions are retained
        but clipped to the contig boundaries.
        
        The extracted regions are written to DEREP_IN_DIR for use in the dereplication step.
        
        Mutates:
            Writes extracted region sequences to temporary files in DEREP_IN_DIR.
        
        Returns:
            None
        """
        
        regions = [r.to_dict() for i,r in self.binary_df.iterrows()]
        contig_ends = thread_map(lambda x: _extract_one_region(x, self.margin, self.TEMP_GENOME_DIR, 
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
    
    
    def join_dereplication_with_binary(self) -> None:
        """
        Join dereplication clustering results with the binary table.
        
        Reads the MMseqs2 clustering output file and joins it with the binary table based on
        scaffold ID and region coordinates (Start, End). This associates each extracted region
        in the binary table with its dereplication status (representative or redundant) and
        representative region identifier. The resulting table is sorted by representative and
        dereplication status for clarity.
        
        The clustering table is parsed to extract scaffold and coordinate information from
        compound region identifiers. Dereplication status is determined by comparing each
        region's identifier with its assigned representative.
        
        Mutates:
            self.binary_df (pd.DataFrame): Updated in-place with additional columns for
                'representative' and 'dereplication_status', and sorted by these columns.
                The temporary 'Region' column is removed after processing.
        
        Returns:
            None
            
        Notes:
            The Region temporary column in self.binary_df was added by joining with the MMseqs2 dereplication table.
            This is the workflow-specific implementation of the abstract method inherited from its grandparent class Run.
        """
        
        # Parse the MMseqs clustering table
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / "derep_cluster.tsv"
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 names = ['representative', 'Region'],
                                 header = None
                                 )

        # Determine dereplication status
        derep_df[['Scaffold', 'Start', 'End']] = derep_df['Region'].str.split(pat = "|", expand = True)
        derep_df[['Start', 'End']] = derep_df[['Start', 'End']].astype(int)
        derep_df['dereplication_status'] = derep_df['Region'] == derep_df['representative']
        derep_df['dereplication_status'] = np.where(derep_df['dereplication_status'], 'dereplication_representative', 'redundant')
        
        # Add dereplication status columns to binary table by joining with the clustering table
        LOG.debug("Joining MMseqs2 clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.merge(derep_df, on = ["Scaffold", 'Start', 'End'])
        self.binary_df.drop(columns = ['Region'], inplace = True)
        
        # Sort by representative ID and then by dereplication status
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])           
        LOG.info("Mapping done!")
        
        return None
    
    
    def run(self):
        """
        Execute the complete local region-based dereplication pipeline.
        
        Orchestrates all processing steps in sequence: stages local genome assemblies for region extraction,
        extracts genomic regions of interest from the staged genomes, runs MMseqs2 dereplication on the
        extracted regions, joins the dereplication clustering results with the binary table, recovers hit 
        diversity information from the dereplication output, filters the original session based on dereplication
        results, and generates final output files with dereplication metadata.
        
        Cleans up temporary working directories upon completion.
        
        Returns:
            None
        """
        
        
        LOG.info("--- STEP 1: Staging genomes for dereplication. ---")
        self.prepare_genomes()
        
        LOG.info("--- STEP 2: Extracting genomic regions. ---")
        self.extract_regions()
        
        LOG.info("--- STEP 3: Dereplicating. ---")
        self.dereplicate_regions()
        
        LOG.info("--- STEP 4: Mapping dereplication output to binary table. ---")
        self.join_dereplication_with_binary()
        
        LOG.info("--- STEP 5: Recovering hit diversity. ---")
        self.recover_hits()
        
        LOG.info("--- STEP 6: Filtering session file. ---")
        self.filter_session()
        
        LOG.info("--- STEP 7: Generating output files")
        self.generate_output()
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
    
        return None
    
    