#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.remote_run import RemoteRun
from cagecleaner.region_run import RegionRun
from cagecleaner.communication import fetch_contig_lengths, download_regions

import logging
import warnings
import pandas as pd
import numpy as np
from pathlib import Path

warnings.filterwarnings(action = "ignore", module = "Bio")
LOG = logging.getLogger()


class RemoteRegionRun(RemoteRun, RegionRun):
    """
    Subclass orchestrating the workflow for dereplication by region sequence similary using regions downloaded from NCBI.
    
    This class combines remote region downloading functionality with local region dereplication.
    It downloads genomic regions (with optional sequence margins) from NCBI based on scaffold IDs and
    coordinates in a binary table, performs MMseqs2-based sequence dereplication on the downloaded regions,
    and integrates dereplication results back into the binary table. Handles contig edge cases
    where regions with sequence margins extend beyond the scaffold boundaries according to user-specified behavior
    (keep but clip them (permissive), or discard them (strict)).
    
    It links NCBI assembly accession IDs to the scaffold IDs in the binary table, downloads the
    corresponding genome assemblies, maps the scaffold IDs back to the downloaded assemblies, performs
    whole-genome dereplication using skDER, and integrates dereplication results back into
    the binary table. Unmapped scaffolds (scaffolds for which no assemby file was found) are tracked
    and reported separately.
    
    Inherits from:
        RemoteRun: Intermediary class providing common downloading configurations.
        RegionRun: Intermediary class providing region dereplication utilities
    """
    
    def __init__(self, args):
        """
        Initialise a RemoteRegionRun instance.
        
        Calls the parent class constructor to set up common configuration, temporary directories,
        and logging infrastructure.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
        
        Returns:
            None
        """
        
        super().__init__(args)
        
        return None
    
    
    def fetch_regions(self) -> None:
        """
        Fetch genomic regions from NCBI based on scaffold coordinates with optional sequence margins.
        
        Extracts region coordinates from the binary table and applies user-specified margins
        (upstream and downstream sequence extensions). Handles regions that extend beyond contig boundaries
        by either discarding them (strict mode) or clipping them (permissive mode). Fetches contig
        lengths from NCBI using Entrez utilities. Logs statistics about discarded or clipped regions.
        Downloads the regions.
        
        Mutates:
            self.binary_df (pd.DataFrame): Adds 'Contig_length' column from NCBI data.
            self.DEREP_IN_DIR: Populates folder with downloaded gzipped genomic region FASTA files.
        
        Returns:
            None
        """
        ## Add sequence margins
        regions = self.binary_df[['Scaffold', 'Start', 'End']].copy()
        regions['Start'] -= self.margin
        regions['End'] += self.margin
        
        ## First fetch the length of each contig
        contig_ids = regions['Scaffold'].to_list()
        lengths_df = fetch_contig_lengths(contig_ids)
        regions = regions.merge(lengths_df, on = "Scaffold", how = 'inner') # Discard the ones without length
        
        # Report the ones without length
        regions_without = regions.merge(lengths_df, on = "Scaffold", how = 'left', indicator = True)
        regions_without = regions_without[regions_without['_merge'] == 'left_only'].drop(columns = ['_merge'])
        if not regions_without.empty:
            LOG.warning(f'Contigs for which no length could be retrieved: {" ".join(regions_without["Scaffold"].to_list())}')
            LOG.warning('These contigs will not be downloaded anyway.')
        
        # Add the contig length column to the extended binary for easy joining later 
        # when merging the dereplication status with the extended binary table
        self.binary_df = self.binary_df.merge(lengths_df, on = 'Scaffold', how = 'left')
        
        ## Treat contig edges upstream (put at zero, or remove)
        mask_up = regions['Start'] <= 0
        if self.strict_regions:
            regions = regions[~mask_up]
        else:
            regions['Start'] = regions['Start'].clip(lower = 1)
            
        ## Treat contig edges downstream
        mask_down = regions['End'] > regions['Contig_length']
        if self.strict_regions:
            regions = regions[~mask_down]
        else:
            regions['End'] = regions['End'].clip(upper = regions['Contig_length'])
            
        ## Some stats
        diff_upstream = mask_up.sum()
        diff_downstream = mask_down.sum()
        diff = diff_upstream + diff_downstream
        LOG.info(f'{diff} regions were at a contig end.')
        if self.strict_regions:
            LOG.info('These regions will not be downloaded in strict mode.')
        
        ## Download the regions
        download_regions(regions, 
                         directory = self.DEREP_IN_DIR, 
                         download_workers = self.download_workers,
                         progress_bar = self.no_progress)
        
        return None
    
    
    def join_dereplication_with_binary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
            
        Join MMseqs2 dereplication clustering results with the binary table.
        
        Reads the dereplication clustering output file and parses region coordinates. Assigns
        dereplication status ('dereplication_representative' or 'redundant') based on whether
        each region's accession matches its representative. Handles four possible cases of region boundary
        mismatches due to region boundary adjustments during preprocessing:
        
        1. No contig edges: Exact coordinate match after removing margins.
        2. Upstream edge: Match on scaffold and end coordinate (contig was clipped at the upstream edge).
        3. Downstream edge: Match on scaffold and start coordinate (contig was clipped at the downstream edge).
        4. Both edges: Interval containment where the original cluster region is within the dereplicated region
           (contig was clipped at both edges).
        
        Concatenates results from all cases and removes temporary helper columns. Sorts the final
        table by representative ID and dereplication status.
        
        Mutates:
            self.binary_df (pd.DataFrame): Joined with dereplication results; adds 'representative'
                                           and 'dereplication_status' columns; removes 'Contig_length'
                                           and 'Region' columns.
        
        Returns:
            None
            
        Notes:
            The Region temporary column in self.binary_df was added by joining with the MMseqs2 dereplication table.
            The Contig_length temporary column in self.binary_df was added by fetch_regions when the sequence margins
            were added on-the-fly when fetching the regions.
            This is the workflow-specific implementation of the abstract method inherited from its grandparent class Run.
        """
        # Read the MMseqs2 clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / "derep_cluster.tsv"
        # Convert to dataframe
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 names = ['representative', 'Region'],
                                 header = None)
        
        # Parse the dereplicated region coordinates (which includes the margin)
        derep_df[['Scaffold', 'Range']] = derep_df['Region'].str.split(pat = ":", expand = True)
        derep_df[['Start', 'End']] = derep_df['Range'].str.split(pat = "-", expand = True)
        derep_df[['Start', 'End']] = derep_df[['Start', 'End']].astype(int)
        derep_df.drop(columns = ['Range'], inplace = True)
        
        # Add dereplication status column based on whether the scaffold's ID is the same as the representative's
        # If there is a match between binary_df['assembly_file'] and derep_dfs['assembly'] (its index column), the representative and status is added.
        derep_df['dereplication_status'] = derep_df['Region'] == derep_df['representative']
        derep_df['dereplication_status'] = np.where(derep_df['dereplication_status'], 'dereplication_representative', 'redundant')
        
        # Subtract the region margin again to merge properly with the extended binary table
        binary_merged_two_edges = binary_merged_edge_up = binary_merged_edge_down = binary_merged_no_edges = pd.DataFrame()

        # Case 1: no contig edges
        derep_no_edges = derep_df.copy()
        key = ['Scaffold', 'Start', 'End']
        derep_no_edges['Start'] += self.margin
        derep_no_edges['End'] -= self.margin
        binary_merged_no_edges = self.binary_df.merge(derep_no_edges, on = key, how = 'inner')
        remaining_binary = self.binary_df.merge(derep_no_edges[key], on = key, how = 'left', indicator = True)
        remaining_binary = remaining_binary[remaining_binary['_merge'] == 'left_only'].drop(columns = ['_merge'])
        
        # Case 2: contig edge upstream
        if not remaining_binary.empty:
            derep_edge_up = derep_df.copy()
            key = ['Scaffold', 'End']
            derep_edge_up['End'] -= self.margin
            binary_merged_edge_up = remaining_binary.merge(derep_edge_up, on = key, how = 'inner')
            remaining_binary = remaining_binary.merge(binary_merged_edge_up[key], on = key, how = 'left', indicator = True)
            remaining_binary = remaining_binary[remaining_binary['_merge'] == 'left_only'].drop(columns = ['_merge'])
        
        # Case 3: contig edge downstream
        if not remaining_binary.empty:
            derep_edge_down = derep_df.copy()
            key = ['Scaffold', 'Start']
            derep_edge_down['Start'] += self.margin
            binary_merged_edge_down = remaining_binary.merge(derep_edge_down, on = key, how = 'inner')
            remaining_binary = remaining_binary.merge(binary_merged_edge_down[key], on = key, how = 'left', indicator = True)
            remaining_binary = remaining_binary[remaining_binary['_merge'] == 'left_only'].drop(columns = ['_merge'])
        
        # Case 4: two contig edges
        if not remaining_binary.empty:
            derep_two_edges = derep_df.copy()
            remaining_binary['Interval'] = [pd.Interval(left = crds.Start, right = crds.End, closed = 'both')
                                            for crds in remaining_binary[['Start', 'End']].itertuples()]
            derep_two_edges['Interval'] = [pd.Interval(left = drp.Start, right = drp.End, closed = 'both')
                                           for drp in derep_two_edges[['Start', 'End']].itertuples()]
            binary_merged_two_edges = remaining_binary.merge(derep_two_edges, on = 'Scaffold', how = 'inner',
                                                             suffixes = ['_bin', '_drep'])
            binary_merged_two_edges = binary_merged_two_edges[binary_merged_two_edges.apply(lambda row:
                                                                                            row['Interval_bin'] 
                                                                                            in row['Interval_drep'],
                                                                                            axis = 1)]

        # Concat the hits from all cases
        binary_merged = pd.concat([binary_merged_no_edges, binary_merged_edge_up, 
                                   binary_merged_edge_down, binary_merged_two_edges])
        self.binary_df = binary_merged
        
        # Sort by representative ID and then by dereplication status
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])
        self.binary_df.drop(columns = ['Contig_length', 'Region'], inplace = True)
        
        LOG.info("Mapping done!")
            
        return None
    
    
    def run(self):
        """
        Execute the complete remote region dereplication pipeline.
        
        Orchestrates all processing steps in sequence: fetches genomic regions from NCBI with
        optional sequence margins and contig boundary handling, performs MMseqs2-based dereplication on
        the downloaded regions, joins the dereplication results with the binary table, recovers hit
        diversity according to user parameters, filters the original session file, and generates
        output files. Cleans up the temporary directory upon completion.
        
        Returns:
            None
        """
        
        LOG.info('--- STEP 1: Downloading genomic regions. ---')
        self.fetch_regions() # Download genomic regions using ncbi_acc_download
        
        LOG.info("--- STEP 2: Dereplicating ---")
        self.dereplicate_regions()  # Dereplicate.
        
        LOG.info("--- STEP 3: Mapping dereplication clustering to binary table ---")
        self.join_dereplication_with_binary()  # Map each row in the binary table with its representative
        
        LOG.info("--- STEP 4: Recovering hit diversity ---")
        self.recover_hits()  # Recover hits by content and score depending on user input
        
        LOG.info("--- STEP 5: Filtering original session file. ---")
        self.filter_session()  # Filter the original session file, retaining only the dereplicated hits.
        
        LOG.info("--- STEP 6: Generating output files ---")
        self.generate_output()  # Generate all output files.
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
