#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.remote_run import RemoteRun
from cagecleaner.region_run import RegionRun
from cagecleaner.communication import _downloadOneRegion

import logging
import warnings
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from Bio import Entrez

warnings.filterwarnings(action = "ignore", module = "Bio")
LOG = logging.getLogger()


class RemoteRegionRun(RemoteRun, RegionRun):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        return None
    
    def downloadRegions(self) -> None:
        max_attempts = 3
        
        regions = self.binary_df[['Scaffold', 'Start', 'End']].copy()
        regions['Start'] -= self.margin
        regions['End'] += self.margin
        
        ## Treat contig edges upstream (put at zero, or remove)
        mask_up = regions['Start'] <= 0
        diff_upstream = mask_up.sum()
        if self.strict_regions:
            regions = regions[~mask_up]
        else:
            regions['Start'] = regions['Start'].clip(lower = 1)
            
        ## Treat contig edges downstream
        # First get the length of each contig
        LOG.debug('Retrieving contig lengths via Entrez')
        
        contig_ids = regions['Scaffold'].to_list()
        for attempt in range(max_attempts):
            try:
                with Entrez.efetch(db = 'nucleotide', id = contig_ids, rettype = 'docsum') as handle:
                    records = Entrez.read(handle)
            except:
                if attempt+1 < max_attempts:
                    LOG.warning(f'Error getting contig lengths in attempt {attempt+1}. Max. attempts: {max_attempts}')
                else:
                    LOG.error(f'Error getting contig lengths after {max_attempts} attempts. Exiting...')
                    sys.exit()
                    
        # Parse the response
        lengths = [int(s['Length']) for s in records]
        accessions_with_lengths = [str(s['AccessionVersion']) for s in records]
        lengths_df = pd.DataFrame({'Scaffold': accessions_with_lengths, 'Contig_length': lengths})
        lengths_df.drop_duplicates(inplace = True, ignore_index = True)
        regions = regions.merge(lengths_df, on = "Scaffold", how = 'inner') # Discard the ones without length
        
        # Report the ones without length
        regions_without = regions.merge(lengths_df, on = "Scaffold", how = 'left', indicator = True)
        regions_without = regions_without[regions_without['_merge'] == 'left_only'].drop(columns = ['_merge'])
        if not regions_without.empty:
            LOG.warning(f'Contigs for which no length could be retrieved: {" ".join(regions_without["Scaffold"].to_list())}')
            LOG.warning('These contigs have been discarded anyway.')
        
        # Add the contig length column to the extended binary for easy joining later 
        # when merging the dereplication status with the extended binary table
        self.binary_df = self.binary_df.merge(lengths_df, on = 'Scaffold', how = 'left')
        
        # Then treat the edges
        mask_down = regions['End'] > regions['Contig_length']
        diff_downstream = mask_down.sum()
        if self.strict_regions:
            regions = regions[~mask_down]
        else:
            regions['End'] = regions['End'].clip(upper = regions['Contig_length'])
            
        ## Some stats
        diff = diff_upstream + diff_downstream
        LOG.info(f'{diff} regions were at a contig end.')
        if self.strict_regions:
            LOG.info('These regions will not be downloaded.')
        
        ## Download the regions
        accessions = regions['Scaffold'].to_list()
        ranges = (regions['Start'].astype(str) + ':' + regions['End'].astype(str)).to_list()
        regions = list(zip(accessions, ranges))
        
        if self.download_workers > 2:
            max_workers = 2
            LOG.warning("You are requesting too many download workers by NCBI policy. Limiting the number to 2 for compliance.")
        else:
            max_workers = self.download_workers
        
        thread_map(lambda x: _downloadOneRegion(x, self.DEREP_IN_DIR), regions, 
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
        
        # Subtract the region margin again to merge properly with the ext. binary table
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
    
