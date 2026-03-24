#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Internal imports:
from cagecleaner import util

# External libraries:
import pandas as pd
import os
import subprocess
import gzip
import shutil
import logging
from tempfile import TemporaryDirectory
from abc import ABC, abstractmethod
from pathlib import Path
from copy import deepcopy
from scipy.stats import zscore
from random import choice
from Bio import SeqIO
from cblaster.classes import Session
from importlib import resources


LOG = logging.getLogger()


class Run(ABC):
    """
    This is an abstract class representing a typical CAGEcleaner Run.
    A Run is either a LocalRun or a RemoteRun, depending on the mode in which cblaster was executed.
    It contains an initializer, shared by both LocalRun and RemoteRun, that internally parses and stores the arguments given through the command line.
    The static method fromArgs() is used to initialize a LocalRun or RemoteRun based on the mode of the provided session file.
        
    One can visualize the workflow as follows:
        
       -- Remote -> fetchAssemblyIDs() -> downloadGenomes() -> mapAssembliesToBinary() -> EXTENDED_BINARY -> |
      |                                                                                                      |
Run --|                                                                                                       -> dereplicate() -> mapDereplicationOutToBinary() -> recoverHits() -> filterSession() -> generateOutput()
      |                                                                                                      |
       -- Local/HMM ------------------------------------> prepareGenomes() --------------------------------> |
        
        With EXTENDED_BINARY containing links to the corresponding genome file for each row.
        (In local mode this link essentially already exists through the 'Organism' column)
        
    """                                                         
    def __init__(self, args):
        
        # Parse dereplication method
        self.regions: bool = (args.method == "regions")
        
        ## Some defensive checks: ##
        assert args is not None, "No arguments were given. ArgParse object is None."
        assert args.session_file.exists() and args.session_file.is_file(), "Provided session file does not exist or is not a file."
        assert args.genome_dir.exists() and args.genome_dir.is_dir(), "Provided genome directory does not exist or is not a directory."
        assert args.identity <= 100 and ((args.identity >= 82 and not(self.regions)) or (args.identity >= 0 and self.regions)), "Identity threshold should be between 82 and 100 in case of full-genome-based dereplication (see skani documentation) or between 0 and 100 in case of region-based dereplication."
        assert args.coverage >= 0 and args.coverage <= 100, "Coverage threshold should be a number between 0 and 100."
        assert args.zscore_outlier_threshold > 0, "Z-score threshold for recovery should be greater than zero."
        assert args.minimal_score_difference >= 0, "Minimal score difference for recovery cannot be smaller than zero."
        assert args.cores > 0, "Amount of CPU cores to use should be greater than zero."
        assert not(self.regions) or args.margin >= 0, "Region margin cannot be negative when dereplicating regions."

        ## User-defined variables: ##
        # Parse the general arguments:
        self.cores: int = args.cores
        self.verbosity: int = args.verbosity
        
        # Parse IO arguments:
        self.session: Session = Session.from_file(args.session_file.resolve())  # Stores the session file as a Session object
        self.OUT_DIR: Path = args.output_dir.resolve()  # Output directory
        Path(args.temp_dir).resolve().mkdir(parents = True, exist_ok = True)
        self.TEMP_DIR_CONTEXT: TemporaryDirectory = TemporaryDirectory(dir = args.temp_dir, delete = False)  # Temporary directory (within the given temporary directory)
        self.TEMP_DIR = Path(self.TEMP_DIR_CONTEXT.name)
        self.USER_GENOME_DIR: Path = args.genome_dir.resolve()  # Genome directory provided by the user
        self.keep_dereplication: bool = args.keep_dereplication
        self.keep_downloads: bool = args.keep_downloads
        self.keep_intermediate: bool = args.keep_intermediate
        
        # Parse bypass and excluded scaffolds/assemblies:
        self.bypass_organisms: set = {i.strip() for i in args.bypass_organisms.split(',')}  # Converts comma-separated sequence to a set.
        self.bypass_scaffolds: set = {i.strip() for i in args.bypass_scaffolds.split(',')}
        self.excluded_organisms: set = {i.strip() for i in args.excluded_organisms.split(',')}
        self.excluded_scaffolds: set = {i.strip() for i in args.excluded_scaffolds.split(',')}

        # Download arguments:
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
            
        ## Now some non-user defined variables follow: ##
        # Make a binary table file from the session object:
        with (self.TEMP_DIR / 'binary.txt').open('w') as handle:
           self.session.format("binary", delimiter = "\t", fp = handle)

        # Store it internally as a dataframe:     
        self.binary_df = pd.read_table(self.TEMP_DIR / 'binary.txt', 
                                       sep = "\t", 
                                       converters= {'Organism': util.removeSuffixes})  # removeSuffixes only relevant in local mode. 
                  
        # This variable will store the filtered session file, the end result.
        self.filtered_session: Session = None
        
        # Set directory for skDER output and Path to dereplication script
        self.DEREP_OUT_DIR: Path = self.TEMP_DIR / 'derep_out'  # Folder where skder will place its output. The directory will be made by the dereplication script when it is called.
        self.GENOME_DIR: Path = self.TEMP_DIR / 'genomes'  # Default path where genomes will be stored. In local mode this can change to USER_GENOME_DIR.
        self.REGION_DIR: Path = self.TEMP_DIR / 'regions' # Path where the genomic regions will be saved temporarily for region-based dereplication
        
        self.DEREPLICATE_SCRIPT: Path = Path(resources.files(__name__)) / 'dereplicate_assemblies.sh'  # Path to the dereplication script

        # Make working subdirectories. Check for presence and force if ok.
        # Genome directory can already exist
        self.GENOME_DIR.mkdir(parents = True, exist_ok = True)
        
        # Region directory should not exist yet.
        try:
            self.REGION_DIR.mkdir(parents = True)
        except FileExistsError:
            if args.force:
                LOG.warning('Region folder already exists, but it will be overwritten.')
            else:
                LOG.error('Region folder already exists! Rerun with -f to overwrite it.')
                sys.exit()
        
        # Output directory should not exist yet.
        try:
            self.OUT_DIR.mkdir(parents = True)
        except FileExistsError:
            if args.force:
                LOG.warning('Output folder already exists, but it will be overwritten.')
            else:
                LOG.error('Output folder already exists! Rerun with -f to overwrite it.')
                sys.exit()
        

    def dereplicate(self):
        """
        This method is the entry point for the dereplication.
        It calls the appropriate dereplication method for the sequences being dereplicated (full genomes vs. regions).
        """
        if self.regions:
            LOG.info("Extracting genomic regions.")
            self.extractRegions()
            LOG.info("Starting region-based dereplication.")
            self.dereplicateRegions()
        else:
            LOG.info("Starting full genome dereplication.")
            self.dereplicateGenomes()
            
        return None
        
    def dereplicateGenomes(self):
        """
        This method takes the path to a genome folder and dereplicates the genomes using skDER.
        skDER output is stored in TEMP_DIR/derep_out.
        """    
        # Define current working directory:
        # home = os.getcwd()
        # Navigate to the temp directory:
        # os.chdir(self.TEMP_DIR)
        # Initiate the dereplication script:
        LOG.info(f'Dereplicating genomes in {str(self.GENOME_DIR)} with identity cutoff of {str(self.identity)} % and coverage cutoff of {str(self.coverage)} %')
        
        LOG.info("Starting skDER")
        subprocess.run(['skder',
                        '-g', self.GENOME_DIR,
                        '-o', str(self.DEREP_OUT_DIR),
                        '-i', str(self.identity),
                        '-f', str(self.coverage),
                        '-c', str(self.cores),
                        '-d', "low_mem_"*self.low_mem + 'greedy',
                        '-n'
                        ],
                       check = True)
        
        LOG.info("Dereplication done!")
        
        extensions = {'.fna','.fa','.fasta','.fna.gz','.fa.gz','.fasta.gz'}
        paths = [str(p) for p in self.GENOME_DIR.iterdir() if extensions & set(p.suffixes)]
        before = len(paths)
        after = len(list((self.DEREP_OUT_DIR / 'Dereplicated_Representative_Genomes').iterdir()))
        LOG.info(f'{before} genomes were reduced to {after} genomes.')
        # subprocess.run(['bash', str(self.DEREPLICATE_SCRIPT),
        #                 str(self.identity), str(self.coverage), 
        #                 str(self.cores), str(self.GENOME_DIR), 'low_'*self.low_mem + 'mem'], 
        #                check = True)
        # Go back home
        # os.chdir(home)
        
        return None
    
    def extractRegions(self):
        """
        This method extracts the genomic regions surrounding each cluster hit using the specified sequence margin.
        Regions at contig edges are discarded when applying the strict region flag.
        """
        # Make temporary regions directory
        self.REGION_DIR.mkdir(parents = True, exist_ok = True)
        # Loop over all hits in the binary table
        contig_end = 0
        for _, row in self.binary_df.iterrows():
            # Get the ID and coordinates of the genomic region
            assembly = row['assembly_file']
            scaffold = row['Scaffold']
            begin = row['Start'] - self.margin
            end = row['End'] + self.margin
            with gzip.open(self.GENOME_DIR / assembly, "rt") as handle:
                seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
                # If strict, skip regions that are at a contig edge
                # Handle the annoying habit of any2fasta in the local model to omit version digits in scaffold IDs,
                # which results in scaffold KeyErrors.
                try:
                    scaffold_to_extract_from = seqs[scaffold]
                except KeyError:
                    scaffold_to_extract_from = seqs[scaffold.split('.')[0]]
                length = len(scaffold_to_extract_from)
                if end >= length or begin < 0:
                    contig_end += 1
                    if self.strict_regions:
                        continue
                    end = min(end, length)
                    begin = max(0, begin)
                # Extract genomic region from assembly
                region = scaffold_to_extract_from[begin:end]
                # Write in a new compressed fasta file
                with gzip.open(self.REGION_DIR / assembly, "wt") as out_handle:
                    SeqIO.write(region, out_handle, "fasta")

        LOG.info(f'{contig_end} regions were at a contig end.')
        if self.strict_regions:
            LOG.info('These regions have been discarded from the analysis.')
            
        return None
    
    def dereplicateRegions(self):
        """
        This method takes the path to a genomic regions folder and dereplicates them using MMseqs2.
        MMseqs2 output is stored in TEMP_DIR/derep_out.
        """
        subprocess.run(['mmseqs', 'easy-cluster',
                        *[str(p) for p in self.REGION_DIR.glob('*')],
                        str(self.DEREP_OUT_DIR),
                        str(self.DEREP_OUT_DIR / 'tmp'),
                        '--min-seq-id', str(self.identity/100),
                        '-c', str(self.coverage/100),
                        '--threads', str(self.cores),
                        '-v', str(0)],
                       check = True)
        
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
        
        def recoverHitsByScore(df: pd.DataFrame) -> pd.DataFrame():
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
                    # Alter the df that was passed as argument such that recoverByContent has this information. Only do this if that row is not already a dereplication representative, in which case the hit is retained anyway.
                    df.at[index, 'dereplication_status'] = 'readded_by_score' if row['dereplication_status'] != 'dereplication_representative' else row['dereplication_status']
                    # Alter the OG binary df:
                    self.binary_df.at[index, 'dereplication_status'] = 'readded_by_score' if row['dereplication_status'] != 'dereplication_representative' else row['dereplication_status']
        
            return df
        
        
        # If the user is not interested in recovering by content, skip this workflow.
        if self.no_recovery_by_content == True:
            LOG.info("Skipping hit recovery.")
        
        else:   
            if self.no_recovery_by_score == True:
                LOG.info("Skipping hit rcovery by score.")
            # Loop over each representative cluster:
            grouped_by_rep = self.binary_df.groupby('representative')  # Define the grouping
            for representative, cluster in grouped_by_rep:
                LOG.debug(f"Recovering hits in the group of {representative}.")
                LOG.debug(f"-> {len(cluster)} hits in this group")
                # Now we have to create subgroups within each group based on the amount of genes in each cluster:        
                # Loop over the grouping:
                grouped_by_content = cluster.groupby(self.session.queries)
                LOG.debug(f"-> {len(grouped_by_content)} subgroups based on gene cluster composition.")
                for _, group in grouped_by_content:
                    # Now we want to recover by cblaster score:
                    if self.no_recovery_by_score == False:
                        group = recoverHitsByScore(group)  # This group now contains rows that are 'readded_by_score'
                    # Check if a representative is in here. If yes, continue:
                    if 'dereplication_representative' in group['dereplication_status'].to_list():
                        continue
                    else:
                        # If there is no representative, pick a random one (that is not already recovered by score), 
                        # get its index, and change the status in the OG binary dataframe:
                        # First we exclude the ones that were already readded by score:
                        group = group[group['dereplication_status'] != 'readded_by_score']
                        # Update the OG binary table at a random index (choice()) from this group
                        self.binary_df.at[choice(group.index), 'dereplication_status'] = 'readded_by_content'
        
            if self.no_recovery_by_content == False:
                recovered_by_content = len(self.binary_df[self.binary_df['dereplication_status'] == 'readded_by_content']['dereplication_status'].to_list())
                LOG.info(f"Total hits recovered by alternative gene cluster composition: {recovered_by_content}")
                if self.no_recovery_by_score == False:
                    recovered_by_score = len(self.binary_df[self.binary_df['dereplication_status'] == 'readded_by_score']['dereplication_status'].to_list())
                    LOG.info(f"Total hits recovered by outlier cblaster score: {recovered_by_score}")
                    
            self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])

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
                LOG.debug("-> Bypassing")
                continue
            # Now we go to the scaffold level and loop over each scaffold associated with this organism
            for hit_idx, hit in reversed(list(enumerate(org['scaffolds']))):
                # Obtain the prefixed hit for proper exclusion (in local mode, the user provides scaffolds to exclude in the form of <assembly:scaffold> to prevent issues with duplicate scaffold IDs)
                prefixed_scaffold = org_full_name + ':' + hit['accession']
                # If we have to bypass it, jump to the next one in line
                if self.bypass_scaffolds != {''} and prefixed_scaffold.endswith(tuple(self.bypass_scaffolds)):
                    LOG.debug(f"-> Bypassing scaffold {hit['accession']}")
                    continue
                # If it's not in the list of dereplicated scaffolds, remove it.
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
        # Navigate to the output folder:
        os.chdir(self.OUT_DIR)
        
        # Generate the outputs
        # Session file
        LOG.debug("Writing filtered session file.")
        with open("filtered_session.json", "w") as filtered_session_handle:
            self.filtered_session.to_json(fp = filtered_session_handle)
            
        # Binary table
        LOG.debug("Writing filtered binary table.")
        with open("filtered_binary.txt", 'w') as filtered_binary_handle:
            self.filtered_session.format(form = "binary", delimiter = "\t", fp = filtered_binary_handle)
            
        # Summary file
        LOG.debug("Writing filtered summary file.")
        with open("filtered_summary.txt", "w") as filtered_summary_handle:
            self.filtered_session.format(form = "summary", fp = filtered_summary_handle)    
            
        # List of cluster numbers
        LOG.debug("Writing list of retained cluster numbers.")
        filtered_cluster_numbers = [cluster['number'] 
                                    for organism in Session.to_dict(self.filtered_session)['organisms'] 
                                    for scaffold in organism['scaffolds'] 
                                    for cluster in scaffold['clusters']]
        with open("retained_cluster_numbers.txt", "w") as numbers_handle:
            numbers_handle.write(','.join([str(nb) for nb in filtered_cluster_numbers]))
                
        # Cluster sizes:
        LOG.debug("Writing genome cluster sizes.")
        self.binary_df.groupby('representative').size().to_frame(name='cluster_size').to_csv('genome_cluster_sizes.txt', sep='\t')
        
        # Keep temp output
        if self.keep_intermediate or self.keep_downloads:
            LOG.info('Copying downloaded genomes to output folder.')
            shutil.copytree(self.TEMP_DIR / 'genomes', 'genomes', dirs_exist_ok = True, ignore_dangling_symlinks = True)
        if self.keep_intermediate or self.keep_dereplication:
            LOG.info("Copying dereplication results to output folder.")
            shutil.copytree(self.TEMP_DIR / 'derep_out', 'derep_out', dirs_exist_ok = True, ignore_dangling_symlinks = True)
            
        # Extended binary:
        LOG.debug("Writing extended binary.")
        self.binary_df.to_csv('extended_binary.txt', sep='\t', index = False)
        
        LOG.info(f"Finished! Output files can be found in {self.OUT_DIR}.")
                
        return None
    
    @abstractmethod
    def run():
        pass
    
