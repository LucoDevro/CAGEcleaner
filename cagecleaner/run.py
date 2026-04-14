#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.file_utils import remove_suffixes
from cagecleaner.utils import correct_layouts, generate_cblaster_session

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
    """
    Abstract base class orchestrating the complete CAGEcleaner dereplication workflow.
    
    Handles which workflows to call for dereplicating genome mining hits across all
    search modes (local/remote sources, genome/region-based dereplication). In all
    workflows, all hits, their metadata and their dereplication status are recorded
    in a cblaster-style binary table extended with additional columns.
    
    The typical high-level workflow is:
    1. Parse input session and create binary table
    2. Dereplicate sequences (implemented by subclasses)
    3. Map dereplication results to binary table
    4. Recover hits by score or content diversity
    5. Filter session and generate outputs
    
    Note:
        This is the abstract grandparent class with globally shared methods. Check out
        parent classes (GenomeRun, RegionRun, RemoteRun, LocalRun) for partially shared methods.
        Subclasses (LocalGenomeRun, RemoteGenomeRun, LocalRegionRun, RemoteRegionRun) provide concrete
        implementations and inherit from these parent classes through multiple inheritance.
        
    See Also:
        LocalGenomeRun: Full-genome dereplication for hits in local sequences.
        RemoteGenomeRun: Full-genome dereplication for hits in remote sequences.
        LocalRegionRun: Region-based dereplication for hits in local sequences.
        RemoteRegionRun: Region-based dereplication for hits in remote sequences.
    """
    
    def __init__(self, args):
        """
        Initialise a Run instance.
        
        Links the parsed arguments to the correct attributes, sets paths and controls output directories.
        Parses the initial binary table, either from a cblaster session, or from TSV tables.
        Generates layout groups for all hits and corrects them for strand location if necessary.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
                
        Raises:
            AssertionError: If input validation fails (missing files, invalid ranges, etc).
            FileExistsError: If output directory already exists (unless -f flag used).
            
        Returns:
            None
        """
        
        super().__init__()
        
        ## Some defensive checks: ##
        assert args is not None, "No arguments were given. ArgParse object is None."
        assert args.session.exists(), "Provided input does not exist."
        assert args.genome_dir.exists() and args.genome_dir.is_dir(), "Provided genome directory does not exist or is not a directory."
        assert args.coverage >= 0 and args.coverage <= 100, "Coverage threshold should be a number between 0 and 100."
        assert args.zscore_outlier_threshold > 0, "Z-score threshold for recovery should be greater than zero."
        assert args.minimal_score_difference >= 0, "Minimal score difference for recovery cannot be smaller than zero."
        assert args.cores > 0, "Amount of CPU cores to use should be greater than zero."
        
        ## Passing on the CLI arguments
        # Parse the general arguments:
        self.cores: int = args.cores
        self.verbosity: int = args.verbosity
        self.no_progress: bool = args.no_progress
        
        # Parse IO arguments:
        self.keep_dereplication: bool = args.keep_dereplication
        self.keep_downloads: bool = args.keep_downloads
        self.keep_intermediate: bool = args.keep_intermediate
        
        # Define paths and check output folder
        self.OUT_DIR: Path = args.output_dir.resolve()  # Output directory
        Path(args.temp_dir).resolve().mkdir(parents = True, exist_ok = True)
        self.TEMP_DIR_CONTEXT: TemporaryDirectory = TemporaryDirectory(dir = args.temp_dir, delete = False)  # Temporary directory (within the given temporary directory)
        self.TEMP_DIR = Path(self.TEMP_DIR_CONTEXT.name)
        self.USER_GENOME_DIR: Path = args.genome_dir.resolve()  # Genome directory provided by the user
        self.DEREP_OUT_DIR: Path = self.TEMP_DIR / 'derep_out'
        self.DEREP_IN_DIR: Path | None = None # Input directory for the dereplication.
        
        # Output directory should not exist yet, unless flagged.
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
        
        ## Initialise the binary table
        self._initialise_binary(args)
        
        return None
    
    
    def _initialise_binary(self, args):
        """
        Initialise the binary table from the provided input files (cblaster session json, or TSV tables).
        
        Recognises which session parsing mode to use from file extensions and presences. Always generates a 
        cblaster Session instance, which stored in the self.session attribute.
        Then parses the binary table generated from this Session and joins it with more cluster information
        retrieved using the cluster extraction function get_sorted_cluster_hierarchies.
        Generates cluster layout groups from this and corrects them for strand location.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
            
        Note:
            Exits when no valid input files have been found.
        
        Returns:
            None
        """
        ## Get a Session object from the inputs
        match args.session.suffix:
            # If session supplied, just parse it
            case '.json':
                LOG.info('Found cblaster session. Parsing it.')
                self.session: Session = Session.from_file(args.session.resolve())
            case _:
                # If not, generate one from the TSV tables
                if args.session.is_dir() and {'clusters.tsv', 'hits.tsv', 'queries.tsv'} <= {p.name for p in args.session.iterdir()}:
                    LOG.info('Found tsv tables. Generating cblaster session from these.')
                    self.session: Session = generate_cblaster_session(args.session, args.source)
                # Else, we have a problem...
                else:
                    LOG.critical('Invalid input mode! Exiting...')
                    sys.exit()
        
        ## Now get the information we want from the session
        # First generate and parse the binary table
        with (self.TEMP_DIR / 'binary.txt').open('w') as handle:
           self.session.format("binary", delimiter = "\t", fp = handle)
        self.binary_df = pd.read_table(self.TEMP_DIR / 'binary.txt', 
                                       sep = "\t", 
                                       converters= {'Organism': remove_suffixes})  # removeSuffixes only relevant in local mode. 
        os.remove(self.TEMP_DIR / 'binary.txt')
        
        # Then add columns with cluster numbers, layouts and strand positions to the binary table
        all_scaffolds = {sc for _,sc,_ in get_sorted_cluster_hierarchies(self.session, max_clusters=None)}
        scaffold_number_map = pd.DataFrame([{'Scaffold': sc.accession,
                                             'Number': cl.number,
                                             'Start': cl.start,
                                             'End': cl.end,
                                             'Strand': tuple([sbj.strand for sbj in cl.subjects]),
                                             'Layout_group': tuple(cl.indices)}
                                            for sc in all_scaffolds for cl in sc.clusters])
        self.binary_df = self.binary_df.merge(scaffold_number_map, on = ['Scaffold', 'Start', 'End'])
        self.binary_df = correct_layouts(self.binary_df) # Generate layout tuples
        
        return None
    
    
    @abstractmethod
    def join_dereplication_with_binary(self) -> None:
        """
        Joins dereplication results with the binary table.
        
        Parses the dereplication clustering table generated by skDER or MMseqs2.
        Each row gets tagged with:
        - representative: Which genome is the dereplication representative of this assembly?
        - dereplication_status: 'dereplication_representative' or 'redundant'
        
        Subclass implementations vary by mode:
        - Local: Join on genome file names (after suffix removal)
        - Remote: Join on assembly accession numbers
        
        This method must:
        1. Read dereplication output from self.DEREP_OUT_DIR
        2. Parse clustering results (which genomes are reps vs redundant)
        3. Join with self.binary_df on appropriate column
        4. Handle unmatched entries gracefully (log warnings, drop if needed)
        5. Update self.binary_df with new columns
        
        Mutates:
            self.binary_df: Adds 'representative' and 'dereplication_status' columns.
            
        Raises:
            FileNotFoundError: If dereplication output file not found.
            ValueError: If output format is unexpected.
            
        Expected Result:
            self.binary_df should now have columns:
            - representative: Genome ID of the dereplication representative
            - dereplication_status: 'dereplication_representative' | 'redundant'
            
        Returns:
            None
            
        Notes:
            This is an abstract method that is passed on as an abstract class by the parent classes.
            Only the specific classes define the method for their use case.
        """
        pass
    
    
    def recover_hits(self) -> None:
        """
        Recovers redundant hits based on cluster layout abd/or outlier homology scores.
        
        This method attempts to restore hits that were marked as redundant during
        dereplication if they represent genuinely distinct gene cluster compositions
        or if they have outlier homology scores. Recovery occurs
        at two levels:
        
        1. **Score-based recovery**: Within each representative genome's clusters,
           identify hits with outlier cblaster scores (Z-score based). These hits
           may represent functionally important evolutionary variants.
        
        2. **Content-based recovery**: For each cluster layout (set of homologs with a particular
           homolog order), select one hit at random to recover to represent that layout, unless the
           dereplication representative is a member of that set.
        
        Recovery is applied hierarchically: hits are grouped by representative
        and cluster layout, then outliers are identified within these groups.
        
        Args:
            None (uses class attributes)
        
        Mutates:
            self.binary_df: Updates the 'dereplication_status' column for recovered hits:
            - 'readded_by_score': Hit with outlier score in its group
            - 'readded_by_content': Hit representing unique gene cluster layout
            
        Returns:
            None
        """
        
        def recover_hits_by_score(df: pd.DataFrame) -> pd.DataFrame:
            """
            Recover hits with outlier cblaster scores within a structural subgroup.
            
            Auxiliary function that processes a single group of hits (defined by
            representative genome and gene cluster layout). Calculates z-scores
            and identifies outlier hits that differ significantly from the average
            of the modal scores in the group.
            
            Args:
                df (pd.DataFrame): Subset of binary_df for a group hits with common
                representative and cluster layout
            
            Mutates:
                self.binary_df: Updates dereplication_status for outlier hits.
                df: Adds z_score column; updates dereplication_status.
                
            Returns:
                df (pd.DataFrame): An updated dataframe containing rows 'readded_by_score'.
            
            Note:
                This function returns the updated subset to the parent function to
                facilitate it identifying hits to recovery by cluster layout
            """
            # Add a column with the z-scores based on the Score in each row:
            df['z_score'] = zscore(df['Score'])
            # Get the mean modal score. Mean because there migth be multiple modal values
            modal_score = df['Score'].mode().mean()
            # Alter the run's binary df at the indices where the score difference and the z_score pass the thresholds
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
                group = recover_hits_by_score(group)  # This group now contains rows that are 'readded_by_score'
            
            # If there is a dereplication presentative in this group, we can skip picking a layout representative
            # since it's already represented by the dereplication representative
            if 'dereplication_representative' in group['dereplication_status'].to_list():
                continue
            else:
                # If there is no dereplication representative, recover a random one (not already recovered by score), 
                # and change its status in the run's binary dataframe
                
                # First we exclude the ones that were already readded by score (as identified by recover_hits_by_score()):
                group_without_score_recovered = group[group['dereplication_status'] != 'readded_by_score']
                # Update the binary table at a random index (choice()) from this group
                self.binary_df.at[choice(group_without_score_recovered.index), 'dereplication_status'] = 'readded_by_content'
    
        # Log some counts
        recovered_by_content = sum(self.binary_df['dereplication_status'] == 'readded_by_content')
        LOG.info(f"Total hits recovered by alternative gene cluster composition: {recovered_by_content}")
        if not(self.no_recovery_by_score):
            recovered_by_score = sum(self.binary_df['dereplication_status'] == 'readded_by_score')
            LOG.info(f"Total hits recovered by outlier hit score: {recovered_by_score}")
                
        # Sort for a nice grouping in the extended binary table
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status', 'Layout_group'])

        return None
    
    
    def filter_session(self) -> None:
        """
        Filter the original cblaster session to retain only dereplicated hits.
        
        Removes all scaffolds with all hits marked as 'redundant' during dereplication 
        from the original cblaster session object. This produces a cleaned session containing
        only representative hits.
        
        The filtering operates at two levels:
        1. **Scaffold level**: Remove scaffolds without hits in the dereplicated set
        2. **Organism level**: Remove organisms with no remaining scaffolds
        
        This method respects bypass flags set by the user, allowing certain organisms
        or scaffolds to bypass filtering (i.e., always retained).
        
        Args:
            None (uses class attributes)
        
        Mutates:
            self.filtered_session (Session): Newly created filtered session object.
        
        Returns:
            None
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
                LOG.debug(f"-> Bypassing organism {org_full_name}")
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
    
    
    def generate_output(self):
        """
        Generate all output files from the filtered session and results.
        
        Creates the complete set of output files in the specified output directory,
        including the filtered cblaster session, binary table, and summary file.
        
        
        Output files generated (always):
        - filtered_session.json: Filtered cblaster session object (JSON format)
        - filtered_binary.txt: Binary presence/absence table (tab-separated)
        - filtered_summary.txt: Summary for each hit per organism (text file)
        - retained_cluster_numbers.txt: Comma-separated list of retained cluster IDs from the cblaster session
        - cluster_sizes.txt: Number of hits per representative (tab-separated)
        - extended_binary.txt: Binary table with dereplication status annotations
        
        Optional intermediate files (if keep flags enabled):
        - downloads/: Copy of downloaded sequences (remote mode only)
        - dereplication/: skDER or MMseqs2 dereplication output files
        
        Args:
            None (uses class attributes)
        
        Returns:
            None
        """
        ## Generate the outputs
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
            
        # Extended binary:
        LOG.debug("Writing extended binary.")
        self.binary_df.to_csv(self.OUT_DIR / 'extended_binary.txt', sep='\t', index = False)
        
        # Cluster sizes:
        LOG.debug("Writing cluster sizes.")
        cluster_sizes = self.binary_df.groupby('representative').size().to_frame(name='cluster_size')
        cluster_sizes.to_csv(self.OUT_DIR / 'cluster_sizes.txt', sep='\t')
        
        ## Copy over requested temporary files
        if self.keep_intermediate or self.keep_downloads:
            LOG.info('Copying downloaded genomes to output folder.')
            shutil.copytree(self.DEREP_IN_DIR, self.OUT_DIR / 'downloads', 
                            dirs_exist_ok = True, ignore_dangling_symlinks = True)
        if self.keep_intermediate or self.keep_dereplication:
            LOG.info("Copying dereplication results to output folder.")
            shutil.copytree(self.DEREP_OUT_DIR, self.OUT_DIR / 'dereplication', 
                            dirs_exist_ok = True, ignore_dangling_symlinks = True)
            
        LOG.info(f"Finished! Output files can be found in {self.OUT_DIR}.")
            
        return None
    
    
    @abstractmethod
    def run():
        """
        Execute the complete CAGEcleaner dereplication and recovery workflow.
        
        This is the main orchestration method that coordinates all workflow steps
        in the correct sequence. Subclasses must implement this to define their
        specific mode (local/remote, genome/region-based).
        
        The typical workflow executed by subclasses is:
        1. Prepare inputs: Validate and stage genomes/regions for dereplication
        2. Dereplication: Run skDER (genome) or MMseqs2 (region) clustering
        3. Integrate dereplication: Join dereplication output with binary table
        4. Hit recovery: Restore hits with outlier scores or unique gene content
        5. Session filtering: Remove redundant scaffolds from session
        6. Output: Generate all result files
        7. Cleanup: Remove temporary files
        
        Implementation details differ by mode:
        - Local/Genome: Use local sequences and skDER for full-genome ANI clustering
        - Local/Region: Use local sequences and MMseqs2 for region sequence clustering
        - Remote/Genome: Download genomes from NCBI Assembly and use skDER ANI clustering
        - Remote/Region: Download regions from NCBI Nucleotide and use MMseqs2 sequence clustering
        
        Args:
            None (uses class attributes)
            
        Mutates:
            self.binary_df: Updated through dereplication and recovery steps.
            self.filtered_session: Populated by filter_session().
            
        Returns:
            None
            
        Expected Results After Execution:
            - self.binary_df: Contains dereplication_status column with values:
                * 'dereplication_representative': Genome/region kept as representative
                * 'redundant': Marked for removal
                * 'readded_by_score': Recovered due to outlier score
                * 'readded_by_content': Recovered due to unique gene content
            
            - self.filtered_session: Filtered Session object ready for output
                (may be None if no hits/clusters identified)
            
            - OUT_DIR: Contains all output files (see generate_output for details)
            
            - TEMP_DIR: Contains intermediate files
        """
        pass
    
    