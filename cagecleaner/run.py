#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Internal imports:
from cagecleaner.file_utils import removeSuffixes
from cagecleaner.utils import correctLayouts, generate_cblaster_session

# External libraries:
import pandas as pd
import os
import shutil
import logging
import sys
from tempfile import TemporaryDirectory
from abc import ABC, abstractmethod
from pathlib import Path
from copy import deepcopy
from scipy.stats import zscore
from random import choice
from cblaster.classes import Session
from cblaster.extract_clusters import get_sorted_cluster_hierarchies


LOG = logging.getLogger()


class Run(ABC):
    
    def __init__(self, args):
        
        super().__init__()
        
        ## Some defensive checks: ##
        assert args is not None, "No arguments were given. ArgParse object is None."
        assert args.session.exists(), "Provided input does not exist."
        assert args.genome_dir.exists() and args.genome_dir.is_dir(), "Provided genome directory does not exist or is not a directory."
        assert args.coverage >= 0 and args.coverage <= 100, "Coverage threshold should be a number between 0 and 100."
        assert args.zscore_outlier_threshold > 0, "Z-score threshold for recovery should be greater than zero."
        assert args.minimal_score_difference >= 0, "Minimal score difference for recovery cannot be smaller than zero."
        assert args.cores > 0, "Amount of CPU cores to use should be greater than zero."

        ## User-defined variables: ##
        # Parse the general arguments:
        self.cores: int = args.cores
        self.verbosity: int = args.verbosity
        self.no_progress: bool = args.no_progress
        
        # Parse IO arguments:
        self.keep_dereplication: bool = args.keep_dereplication
        self.keep_downloads: bool = args.keep_downloads
        self.keep_intermediate: bool = args.keep_intermediate
        
        # Define paths and check output folders
        self.OUT_DIR: Path = args.output_dir.resolve()  # Output directory
        Path(args.temp_dir).resolve().mkdir(parents = True, exist_ok = True)
        self.TEMP_DIR_CONTEXT: TemporaryDirectory = TemporaryDirectory(dir = args.temp_dir, delete = False)  # Temporary directory (within the given temporary directory)
        self.TEMP_DIR = Path(self.TEMP_DIR_CONTEXT.name)
        self.USER_GENOME_DIR: Path = args.genome_dir.resolve()  # Genome directory provided by the user
        self.DEREP_OUT_DIR: Path = self.TEMP_DIR / 'derep_out'  # Folder where skder will place its output. The directory will be made by the dereplication script when it is called.
        # Output directory should not exist yet.
        try:
            self.OUT_DIR.mkdir(parents = True)
        except FileExistsError:
            if args.force:
                LOG.warning('Output folder already exists, but it will be overwritten.')
            else:
                LOG.error('Output folder already exists! Rerun with -f to overwrite it.')
                sys.exit()
        
        # Parse bypass and excluded scaffolds/assemblies:
        self.bypass_organisms: set = {i.strip() for i in args.bypass_organisms.split(',')}  # Converts comma-separated sequence to a set.
        self.bypass_scaffolds: set = {i.strip() for i in args.bypass_scaffolds.split(',')}
        self.excluded_organisms: set = {i.strip() for i in args.excluded_organisms.split(',')}
        self.excluded_scaffolds: set = {i.strip() for i in args.excluded_scaffolds.split(',')}

        # Download arguments:
        self.download_workers: int = args.download_workers
        self.download_batch: int = args.download_batch
        
        # Dereplication arguments:
        self.strict_regions: bool = args.strict_regions
        self.identity: float = args.identity
        self.coverage: float = args.coverage
        self.low_mem: bool = args.low_mem
        self.margin: int = args.margin
        
        # Hit recovery arguments:
        self.no_recovery_by_content: bool = args.no_recovery_by_content
        self.no_recovery_by_score: bool = args.no_recovery_by_score
        self.zscore_outlier_threshold: float = args.zscore_outlier_threshold
        self.minimal_score_difference: float = args.minimal_score_difference
        
        # This variable will store the filtered session file, the end result.
        self.filtered_session: Session | None = None
        
        ## Define the cblaster session
        match args.session.suffix:
            # If supplied, just parse it
            case '.json':
                LOG.info('Found cblaster session. Parsing it.')
                self.session: Session = Session.from_file(args.session.resolve())
            # If not supplied, generate one from the alternative hit and cluster tables
            case _:
                if args.session.is_dir() and '.tsv' in [f.suffix for f in args.session.iterdir()]:
                    LOG.info('Found hit and cluster tables. Generating cblaster session from these.')
                    self.session: Session = generate_cblaster_session(args.session, args.mode)
                else:
                    LOG.critical('Invalid input mode! Exiting...')
                    sys.exit()
        
        ## Now get the information we want from the session file
        # First parse the binary table
        with (self.TEMP_DIR / 'binary.txt').open('w') as handle:
           self.session.format("binary", delimiter = "\t", fp = handle)
        self.binary_df = pd.read_table(self.TEMP_DIR / 'binary.txt', 
                                       sep = "\t", 
                                       converters= {'Organism': removeSuffixes})  # removeSuffixes only relevant in local mode. 
        os.remove(self.TEMP_DIR / 'binary.txt')
        
        # Then join the binary table with cluster numbers, layouts and strand positions
        all_scaffolds = [sc for _,sc,_ in get_sorted_cluster_hierarchies(self.session, max_clusters=None)]
        scaffold_number_map = pd.DataFrame([{'Scaffold': sc.accession,
                                             'Number': cl.number,
                                             'Start': cl.start,
                                             'End': cl.end,
                                             'Strand': tuple([sbj.strand for sbj in cl.subjects]),
                                             'Layout_group': tuple(cl.indices)}
                                            for sc in all_scaffolds for cl in sc.clusters])
        self.binary_df = self.binary_df.merge(scaffold_number_map, on = ['Scaffold', 'Start', 'End'])
        
        # Check whether layouts are truly different when seen from the complementary strand
        # and correct if necessary
        self.binary_df = correctLayouts(self.binary_df)
        
        return None
    
    @abstractmethod
    def mapDereplicationToBinary(self):
        pass
    
    def recoverHits(self) -> None:
        """
        This function groups the binary_df by representatives.
        This grouping is then further subdivided into groups based on the gene absence/presence (gene cluster content)
        Within each structural grouping, outliers are recoverd based on z-scores and a representative is picked (unless the dereplication representative is part of that structural subgroup).
        
        Mutates:
            self.binary_df: pd.DataFrame: Binary table derived from the cblaster Session object.
        """
        
        def recoverHitsByScore(df: pd.DataFrame) -> pd.DataFrame:
            """
            Auxiliary function that takes a dataframe (structural subgroup in this context) and calculates z_scores and minimal difference. 
            It then alters the OG binary table at the correct index and returns the updated grouping to recoverByContent()
            
            Input:
                df: pd.DataFrame: A pandas dataframe. In this case, a two-layer grouping from the binary table (first grouped on representative, and then on gene cluster content)
            
            Mutates:
                self.binary_df: pd.DataFrame: Binary table derived from the cblaster Session object.
                
            Returns:
                pd.DataFrame: An updated dataframe containing rows 'readded_by_score'
            """
            # Add a column with the z-scores based on the Score in each row:
            df['z_score'] = zscore(df['Score'])
            # Get the mean modal score. Mean because there migth be multiple modal values
            modal_score = df['Score'].mode().mean()
            # Alter the OG binary df at the indices where the score difference and the z_score pass the thresholds
            for index, row in df.iterrows():
                if (abs(row['Score'] - modal_score) >= self.minimal_score_difference) and (abs(row['z_score']) >= self.zscore_outlier_threshold):
                    # Alter the df that was passed as argument such that recoverByContent has this information.
                    # Only do this if that row is not already a dereplication representative, 
                    # in which case the hit is retained anyway.
                    df.at[index, 'dereplication_status'] = 'readded_by_score' if row['dereplication_status'] != 'dereplication_representative' else row['dereplication_status']
                    
                    # Alter the OG binary df:
                    self.binary_df.at[index, 'dereplication_status'] = 'readded_by_score' if row['dereplication_status'] != 'dereplication_representative' else row['dereplication_status']
        
            return df
        
        
        # If the user is not interested in recovering by content, skip this workflow.
        if self.no_recovery_by_content:
            LOG.info("Skipping hit recovery.")
            return None
        
        if self.no_recovery_by_score:
            LOG.info("Skipping hit recovery by score.")
        
        # Group by representative and layout group
        grouped_by_rep_layout = self.binary_df.groupby(['representative', 'Layout_group'])
        for rep_layout, group in grouped_by_rep_layout:
            LOG.debug(f"Recovering hits in the group of {rep_layout[0]} with layout {rep_layout[1]}.")
            LOG.debug(f"-> {len(group)} hits in this group")
            
            # Now we want to recover hits with outlier scores within this group of hits with same representative and cluster layout
            if not(self.no_recovery_by_score):
                group = recoverHitsByScore(group)  # This group now contains rows that are 'readded_by_score'
            
            # If there is a dereplication presentative in this group, we can skip it since it's already represented
            if 'dereplication_representative' in group['dereplication_status'].to_list():
                continue
            else:
                # If there is no dereplication representative, recover a random one (not already recovered by score), 
                # and change its status in the OG binary dataframe
                
                # First we exclude the ones that were already readded by score:
                group_without_score_recovered = group[group['dereplication_status'] != 'readded_by_score']
                # Update the OG binary table at a random index (choice()) from this group
                self.binary_df.at[choice(group_without_score_recovered.index), 'dereplication_status'] = 'readded_by_content'
    
        # Log some counts
        recovered_by_content = sum(self.binary_df['dereplication_status'] == 'readded_by_content')
        LOG.info(f"Total hits recovered by alternative gene cluster composition: {recovered_by_content}")
        if not(self.no_recovery_by_score):
            recovered_by_score = sum(self.binary_df['dereplication_status'] == 'readded_by_score')
            LOG.info(f"Total hits recovered by outlier hit score: {recovered_by_score}")
                
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status', 'Layout_group'])

        return None
    
    def filterSession(self) -> None:
        """
        Filter the original session file, keeping only those entries that are not marked as 'redundant' after dereplication.
        
        Mutates:
            self.filtered_session: Session: The filtered Session object.
        """

        dereplicated_scaffolds = [scaffold.strip() for scaffold in self.binary_df[self.binary_df['dereplication_status'] != 'redundant']['Scaffold'].to_list()]
        
        # Get a dictionary export of the session object
        session_dict = self.session.to_dict()
        
        # Make a deep copy to start carving out the new session
        filtered_session_dict = deepcopy(session_dict)
        
        ## Filtering the json session file
        # Remove all cluster hits that do not link with a dereplicated scaffold
        # Remove strains that have no hits left
        
        scaffolds_removed_total = 0  # Counter to keep track of deleted hits
        # Loop over the list of organisms. This list contains dictionaries
        # Going in reverse so the index doesn't change by popping items
        for org_idx, org in reversed(list(enumerate(session_dict['organisms']))):
            scaffolds_removed = 0  # Counter
            
            # Get the full name of the organism, with strain (if it is not empty):
            org_full_name = f"{org['name']} {org['strain']}".strip()
            LOG.debug(f"Carving out organism {org_full_name}")
            
            # First we check whether we need to bypass this assembly. If yes, don't pop this one and proceed with the next one
            if self.bypass_organisms != {''} and org_full_name in self.bypass_organisms:
                LOG.debug("-> Bypassing organism {org_full_name}")
                continue
            
            # Now we go to the scaffold level and loop over each scaffold associated with this organism
            for hit_idx, hit in reversed(list(enumerate(org['scaffolds']))):
                # Obtain the prefixed hit for proper exclusion (in local mode, the user provides scaffolds to exclude in the form of <assembly:scaffold> to prevent issues with duplicate scaffold IDs)
                prefixed_scaffold = org_full_name + ':' + hit['accession']
                
                # If we have to bypass it, jump to the next one in line
                if self.bypass_scaffolds != {''} and prefixed_scaffold.endswith(tuple(self.bypass_scaffolds)):
                    LOG.debug(f"-> Bypassing scaffold {hit['accession']}")
                    continue
                
                # Remove it when it's not in the list of dereplicated scaffolds
                elif hit['accession'] not in dereplicated_scaffolds:
                    filtered_session_dict['organisms'][org_idx]['scaffolds'].pop(hit_idx)
                    scaffolds_removed += 1  # Update counter
                    
            # Tell the user how many scaffolds were removed in this organism:
            LOG.debug(f"-> Removed {scaffolds_removed} scaffolds for {org_full_name}")
            
            # If the scaffold list for this organism is now empty, remove the entire organism:
            if not filtered_session_dict['organisms'][org_idx]['scaffolds']:
                filtered_session_dict['organisms'].pop(org_idx)
            
            # Update the total counter
            scaffolds_removed_total += scaffolds_removed  
        
        # Store the filtered session internally:
        self.filtered_session = Session.from_dict(filtered_session_dict)
        
        LOG.info("Filtering done.")
        LOG.info(f"{scaffolds_removed_total} redundant scaffolds have been removed.")
        LOG.info(f"{len(dereplicated_scaffolds)} scaffolds have been retained.")

        return None
    
    def generateOutput(self):
        """
        Generate all output files:
            - Filtered_session
            - Filtered_binary
            - Filtered_summary
            - List of retained cluster numbers
            - Cluster sizes for each representative genome
            - Extended binary table
        
        Optional:
            - skDER output
            - Downloaded genomes
        """
        # Generate the outputs
        # Session file
        LOG.debug("Writing filtered session file.")
        with open(self.OUT_DIR / "filtered_session.json", "w") as filtered_session_handle:
            self.filtered_session.to_json(fp = filtered_session_handle)
            
        # Binary table
        LOG.debug("Writing filtered binary table.")
        with open(self.OUT_DIR / "filtered_binary.txt", 'w') as filtered_binary_handle:
            self.filtered_session.format(form = "binary", delimiter = "\t", fp = filtered_binary_handle)
            
        # Summary file
        LOG.debug("Writing filtered summary file.")
        with open(self.OUT_DIR / "filtered_summary.txt", "w") as filtered_summary_handle:
            self.filtered_session.format(form = "summary", fp = filtered_summary_handle)    
            
        # List of cluster numbers
        LOG.debug("Writing list of retained cluster numbers.")
        filtered_cluster_numbers = [cluster['number'] 
                                    for organism in Session.to_dict(self.filtered_session)['organisms'] 
                                    for scaffold in organism['scaffolds'] 
                                    for cluster in scaffold['clusters']]
        with open(self.OUT_DIR / "retained_cluster_numbers.txt", "w") as numbers_handle:
            numbers_handle.write(','.join([str(nb) for nb in filtered_cluster_numbers]))
                
        # Cluster sizes:
        LOG.debug("Writing cluster sizes.")
        cluster_sizes = self.binary_df.groupby('representative').size().to_frame(name='cluster_size')
        cluster_sizes.to_csv(self.OUT_DIR / 'cluster_sizes.txt', sep='\t')
        
        # Keep temp output
        if self.keep_intermediate or self.keep_downloads:
            LOG.info('Copying downloaded genomes to output folder.')
            shutil.copytree(self.DEREP_IN_DIR, self.OUT_DIR / 'downloads', 
                            dirs_exist_ok = True, ignore_dangling_symlinks = True)
        if self.keep_intermediate or self.keep_dereplication:
            LOG.info("Copying dereplication results to output folder.")
            shutil.copytree(self.DEREP_OUT_DIR, self.OUT_DIR / 'dereplication', 
                            dirs_exist_ok = True, ignore_dangling_symlinks = True)
            
        # Extended binary:
        LOG.debug("Writing extended binary.")
        self.binary_df.to_csv(self.OUT_DIR / 'extended_binary.txt', sep='\t', index = False)
        
        LOG.info(f"Finished! Output files can be found in {self.OUT_DIR}.")
                
        return None
    
    
    @abstractmethod
    def run():
        pass
    
    