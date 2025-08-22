
# Internal imports:
from cagecleaner import util

# External libraries:
import pandas as pd
import os
import sys
import subprocess
import gzip
import re
import shutil
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from copy import deepcopy
from scipy.stats import zscore
from random import choice
from itertools import batched
from Bio import SeqIO
from cblaster.classes import Session
from importlib import resources


class Run(ABC):
    """
    This is an abstract class representing a typical CAGEcleaner Run.
    A Run is either a LocalRun or a RemoteRun, depending on the mode in which cblaster was executed.
    It contains an initializer, shared by both LocalRun and RemoteRun, that internally parses and stores the arguments given through the command line.
    The static method fromArgs() is used to initialize a LocalRun or RemoteRun based on the mode of the provided session file.
        
    One can visualize the workflow as follows:
        
       -- Remote -> fetchAssemblyIDs() -> downloadGenomes() -> mapAssembliesToBinary() -> EXTENDED_BINARY -> |
      |                                                                                                      |
Run --|                                                                                                       -> dereplicateGenomes() -> mapSkderOutToBinary() -> recoverHits() -> filterSession() -> generateOutput()
      |                                                                                                      |
       -- Local/HMM ------------------------------------> prepareGenomes() --------------------------------> |
        
        With EXTENDED_BINARY containing links to the corresponding genome file for each row.
        (In local mode this link essentially already exists through the 'Organism' column)
        
    """                                                         
    def __init__(self, args):
        
        ## Some defensive checks: ##
        assert args is not None, "No arguments were given. ArgParse object is None."
        assert args.session_file.exists() and args.session_file.is_file(), "Provided session file does not exist or is not a file."
        assert args.output_dir.exists() and args.output_dir.is_dir(), "Provided output directory does not exist or is not a directory."
        assert args.temp_dir.exists() and args.temp_dir.is_dir(), "Provided temp directory does not exist or is not a directory."
        assert args.genome_dir.exists() and args.genome_dir.is_dir(), "Provided genome directory does not exist or is not a directory."
        assert args.ani >= 82 and args.ani <= 100, "ANI threshold should be a number between 82 and 100 (see skani documentation)."
        assert args.zscore_outlier_threshold > 0, "Z-score threshold for recovery should be greater than zero."
        assert args.minimal_score_difference >= 0, "Minimal score difference for recovery cannot be smaller than zero."
        assert args.cores > 0, "Amount of CPU cores to use should be greater than zero."
        
        ## User-defined variables: ##
        # Parse the general arguments:
        self.cores: int = args.cores   
        self.verbose: bool = args.verbose
        
        # Parse IO arguments:
        self.session: Session = Session.from_file(args.session_file.resolve())  # Stores the session file as a Session object
        self.OUT_DIR: Path = args.output_dir.resolve()  # Output directory
        self.TEMP_DIR: Path = Path(tempfile.TemporaryDirectory(dir = args.temp_dir).name).resolve()  # Temporary directory (within the given temporary directory)
        self.TEMP_DIR.mkdir(exist_ok=True)
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
        self.ani: float = args.ani
        
        # Hit recovery arguments:
        self.no_recovery_by_content: bool = args.no_recovery_by_content
        self.no_recovery_by_score: bool = args.no_recovery_by_score
        self.zscore_outlier_threshold: float = args.zscore_outlier_threshold
        self.minimal_score_difference: float = args.minimal_score_difference
            
        ## Now some non-user defined variables follow: ##
        # Make a binary table file from the session object:
        with (self.TEMP_DIR / 'binary.txt').open('w') as handle:
           self.session.format("binary", fp = handle)

        # Store it internally as a dataframe:     
        self.binary_df = pd.read_table(self.TEMP_DIR / 'binary.txt', 
                                       sep = "\\s{2,}", 
                                       engine = 'python',
                                       converters= {'Organism': util.removeSuffixes})  # removeSuffixes only relevant in local mode. 
                  
        # This variable will store the filtered session file, the end result.
        self.filtered_session: Session = None
        
        # Set directory for skDER output and Path to dereplication script
        self.SKDER_OUT_DIR: Path = self.TEMP_DIR / 'skder_out'  # Folder where skder will place its output. The directory will be made by the dereplication script when it is called.
        self.GENOME_DIR: Path = self.TEMP_DIR / 'genomes'  # Default path where genomes will be stored. In local mode this can change to USER_GENOME_DIR.
        self.DEREPLICATE_SCRIPT: Path = Path(resources.files(__name__)) / 'dereplicate_assemblies.sh'  # Path to the dereplication script

    @staticmethod
    def fromArgs(args):
        """
        This is the first entry point for initializing a run.
        We read the session file, determine its mode, and return a corresponding LocalRun or RemoteRun
        """
        print("\n--- Loading session file. ---")
        mode = Session.from_file(args.session_file).params['mode'] 
        
        if mode == 'local' or mode == 'hmm':
            print(f"Detected {mode} mode.")
            return LocalRun(args)
        
        elif mode == 'remote':
            print("Detected remote mode.")
            return RemoteRun(args)
        
        else:
            print(f"CAGEcleaner does not support cblaster {mode} mode for the moment. Exiting the program.")
            sys.exit()
            
    def dereplicateGenomes(self):
        """
        This method takes the path to a genome folder and dereplicates the genomes using skDER.
        skDER output is stored in TEMP_DIR/skder_out.
        """    
        # Define current working directory:
        home = os.getcwd()
        # Navigate to the temp directory:
        os.chdir(self.TEMP_DIR)
        # Initiate the dereplication script:
        self.VERBOSE("Calling dereplication script.")
        subprocess.run(['bash', str(self.DEREPLICATE_SCRIPT), str(self.ani), str(self.cores), str(self.GENOME_DIR)], check = True)
        # Go back home
        os.chdir(home)
        
        return None
    
    @abstractmethod
    def mapSkderOutToBinary(self):
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
            print("Skipping hit recovery.")
        
        else:   
            if self.no_recovery_by_score == True:
                print("Skipping hit rcovery by score.")
            # Loop over each representative cluster:
            grouped_by_rep = self.binary_df.groupby('representative')  # Define the grouping
            for representative, cluster in grouped_by_rep:
                self.VERBOSE(f"Recovering hits in the group of {representative}.")
                self.VERBOSE(f"-> {len(cluster)} hits in this group")
                # Now we have to create subgroups within each group based on the amount of genes in each cluster:        
                # Loop over the grouping:
                grouped_by_content = cluster.groupby(self.session.queries)
                self.VERBOSE(f"-> {len(grouped_by_content)} subgroups based on gene cluster composition.")
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
                        pass
        
            if self.no_recovery_by_content == False:
                recovered_by_content = len(self.binary_df[self.binary_df['dereplication_status'] == 'readded_by_content']['dereplication_status'].to_list())
                print(f"Total hits recovered by alternative gene cluster composition: {recovered_by_content}")
                if self.no_recovery_by_score == False:
                    recovered_by_score = len(self.binary_df[self.binary_df['dereplication_status'] == 'readded_by_score']['dereplication_status'].to_list())
                    print(f"Total hits recovered by outlier cblaster score: {recovered_by_score}")

        return None
    
    def filterSession(self) -> None:
        """
        Filter the original session file, keeping only those entries that are not marked as 'redundant' after dereplication.
        
        Mutates:
            self.filtered_session: Session: The filtered Session object.
        """
        def bypassScaffold(scaffold) -> bool:
            # Returns true if the given scaffold is contained within any of the scaffolds to be excluded.
            return scaffold.endswith(tuple(self.bypass_scaffolds))

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
            self.VERBOSE(f"Carving out organism {org_full_name}")
            # First we check whether we need to bypass this assembly. If yes, don't pop this one and proceed with the next one
            if self.bypass_organisms != {''} and org_full_name in self.bypass_organisms:
                self.VERBOSE("-> Bypassing")
                continue
            # Now we go to the scaffold level and loop over each scaffold associated with this organism
            for hit_idx, hit in reversed(list(enumerate(org['scaffolds']))):
                # Obtain the prefixed hit for proper exclusion (in local mode, the user provides scaffolds to exclude in the form of <assembly:scaffold> to prevent issues with duplicate scaffold IDs)
                prefixed_scaffold = org_full_name + ':' + hit['accession']
                # If we have to bypass it, jump to the next one in line
                if self.bypass_scaffolds != {''} and prefixed_scaffold.endswith(tuple(self.bypass_scaffolds)):
                    self.VERBOSE(f"-> Bypassing scaffold {hit['accession']}")
                    continue
                # If it's not in the list of dereplicated scaffolds, remove it.
                elif hit['accession'] not in dereplicated_scaffolds:
                    filtered_session_dict['organisms'][org_idx]['scaffolds'].pop(hit_idx)
                    scaffolds_removed += 1  # Update counter
            # Tell the user how many scaffolds were removed in this organism:
            self.VERBOSE(f"-> Removed {scaffolds_removed} scaffolds for {org_full_name}")
            # If the scaffold list for this organism is now empty, remove the entire organism:
            if not filtered_session_dict['organisms'][org_idx]['scaffolds']:
                filtered_session_dict['organisms'].pop(org_idx)
            # Update the total counter
            scaffolds_removed_total += scaffolds_removed  
        
        # Store the filtered session internally:
        self.filtered_session = Session.from_dict(filtered_session_dict)
        
        print(f"Filtering done. Removed {scaffolds_removed_total} redundant scaffolds.")

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
        self.VERBOSE("Writing filtered session file.")
        with open("filtered_session.json", "w") as filtered_session_handle:
            self.filtered_session.to_json(fp = filtered_session_handle)
            
        # Binary table
        self.VERBOSE("Writing filtered binary table.")
        with open("filtered_binary.txt", 'w') as filtered_binary_handle:
            self.filtered_session.format(form = "binary", fp = filtered_binary_handle)
            
        # Summary file
        self.VERBOSE("Writing filtered summary file.")
        with open("filtered_summary.txt", "w") as filtered_summary_handle:
            self.filtered_session.format(form = "summary", fp = filtered_summary_handle)    
            
        # List of cluster numbers
        self.VERBOSE("Writing list of retained cluster numbers.")
        filtered_cluster_numbers = [cluster['number'] 
                                    for organism in Session.to_dict(self.filtered_session)['organisms'] 
                                    for scaffold in organism['scaffolds'] 
                                    for cluster in scaffold['clusters']]
        with open("retained_cluster_numbers.txt", "w") as numbers_handle:
            numbers_handle.write(','.join([str(nb) for nb in filtered_cluster_numbers]))
                
        # Cluster sizes:
        self.VERBOSE("Writing genome cluster sizes.")
        self.binary_df.groupby('representative').size().to_frame(name='cluster_size').to_csv('genome_cluster_sizes.txt', sep='\t')
        
        # Keep temp output
        if self.keep_intermediate or self.keep_downloads:
            print('Copying downloaded genomes to output folder.')
            shutil.copytree(self.TEMP_DIR / 'genomes', 'genomes', dirs_exist_ok = True)
        if self.keep_intermediate or self.keep_dereplication:
            print("Copying skDER results to output folder.")
            shutil.copytree(self.TEMP_DIR / 'skder_out', 'skder_out', dirs_exist_ok = True)
            
        # Extended binary:
        self.VERBOSE("Writing extended binary.")
        self.binary_df.to_csv('extended_binary.txt', sep='\t', index = False)
        
        print(f"Finished! Output files can be found in {self.OUT_DIR}.")
                
        return None
    
    @abstractmethod
    def run():
        pass
    
    def VERBOSE(self, message: str) -> None:
        """
        Function to print verbose statements. Prevents having to write an if statement every time.
        """
        if self.verbose:
            print('[VERBOSE] ' + message)
        else:
            pass
    
 #####################################################################################################
   
class LocalRun(Run):
    
    def __init__(self, args):
                
        # Call the parent class initiator
        super().__init__(args)
        
        # Can't keep downloads in local mode:
        assert self.keep_downloads == False, "Can't keep downloads in local mode."
        # Also make sure keep_intermediate doesn't fail:
        if self.keep_intermediate:
            self.keep_intermediate = False
            self.keep_dereplication = True
        
        # Make sure there is no exotic stuff in the provided genome folder
        for file in self.USER_GENOME_DIR.iterdir():
            assert file.is_file(), f"The following object in the genome folder is not a file: {file}. Please move or remove it."
            assert (util.isFasta(str(file)) or util.isGff(str(file)) or util.isGenbank(str(file))), f"The following file is not in the correct format: {file}.\nWe only accept the following suffices: [.fasta, .fna, .fa, .gff3, .gff, .gbff, .gbk, .gb]. The '.gz' extension is allowed."
        
        # Remove organisms that the user wants to be excluded:
        if self.excluded_organisms != {''}:
            self.VERBOSE(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]

        # Remove scaffold IDs specified by the user:
        if self.excluded_scaffolds != {''}:
            self.VERBOSE(f"Excluding the following scaffolds: {', '.join(self.excluded_scaffolds)}")
            # Here the approach slightly differs as users might have provided prefixed scaffold IDs:
            # Add a column with a prefixed scaffold based on the Organism column
            self.binary_df['prefixed_scaffold'] = self.binary_df['Organism'] + ':' + self.binary_df['Scaffold']
            # If the prefixed scaffold ends with any of the strings in the set of scaffolds to exclude, remove it:
            self.binary_df = self.binary_df[self.binary_df['prefixed_scaffold'].str.endswith(tuple(self.excluded_scaffolds)) == False]
            # Clean up:
            self.binary_df = self.binary_df.drop(columns=['prefixed_scaffold'])
            
    def prepareGenomes(self) -> None:
        """
        This function is called to inspect the provided genome folder and convert GenBank files to FASTA files if necessary.
        If a FASTA file is found, it is assumed that all genomes are stored in the provided folder in FASTA format.
        If no FASTA is found, but instead a GenBank is found, it is assumed that all genomes are in GenBank format.
        Temporary FASTA-converted copies of these GenBank files will be placed in the temporary folders.
        Only FASTA files can be passed on to skDER.
        If neither GenBank nor FASTA is found, the program exits.
    
        Mutates:
            self.GENOME_DIR: Path: Path to the genome directory
        """
        
        # Assert that Organism and filenames correspond:
        files_in_genomes_dir: set = {util.removeSuffixes(file) for file in os.listdir(self.USER_GENOME_DIR)}
        organisms_in_session: set = {util.removeSuffixes(organism.name) for organism in self.session.organisms}
        
        assert files_in_genomes_dir == organisms_in_session, "Organisms in session and files in genome directory do not correspond. CAGEcleaner leverages the fact that the Organism column of a cblaster binary table contains the file names of the genomes. Please make sure you havent changed the genome file names between a cblaster run and a CAGEcleaner run. Alternatively, make sure that there are no other files or folders in the provided genome directory."
        
        # Check if there are FASTA files in the genome folder:
        fasta_in_folder = [util.isFasta(str(file)) for file in self.USER_GENOME_DIR.iterdir()]
        genbank_in_folder = [util.isGenbank(str(file)) for file in self.USER_GENOME_DIR.iterdir()]

        if any(fasta_in_folder):
            # In this case the genome folder path should remain the same
            print(f"Detected {sum(fasta_in_folder)} FASTA files in {self.USER_GENOME_DIR}. These will be used for dereplication.")
            # Redirect the genome dir to the user-provided folder:
            self.GENOME_DIR = self.USER_GENOME_DIR
            
        
        elif any(genbank_in_folder):
            # In this case we convert to FASTA and redirect to genome folder, which is in the temp folder by default.
            print(f"Detected {sum(genbank_in_folder)} GenBank files in {self.USER_GENOME_DIR}. Now converting to FASTA for dereplication.")
            # Convert to FASTA files:
            self.GENOME_DIR.mkdir(exist_ok=True)  # Make the output folder if it does not exist already.
            util.convertGenbankToFasta(self.USER_GENOME_DIR, self.GENOME_DIR)
            print(f"Migrated genomes in FASTA format to {self.GENOME_DIR}")
            
        else:
            # If there are no FASTA or GenBank files, the program cannot proceed:
            print("No FASTA files or GenBank files were detected in the provided genome folder. Exiting the program.")
            sys.exit()
        
        return None  
    
    def mapSkderOutToBinary(self) -> None:
        """
        This function maps the skDER clustering table to our binary table.
        Each row in our binary table is coupled to a representative genome and its status (redundant or representative).
        This is done by leveraging the fact that the entries in the 'Organism' column are derived from the file names of the genomes from which they come (given that the user has not changed these names between a cblaster and CAGEcleaner run).
        If an 'Organism' and assembly file name are stripped from their suffices (.gz, .gbk, .fasta...), they should be identical. Thus permitting a simple join operation.
        
        Input:
            self: Current instance of the class.   
     
        Mutates:
            self.binary_df: pd.DataFrame: The binary table derived from a cblaster Session object.
        """
        def extractAssembly(file_path: str) -> str:
            return util.removeSuffixes(os.path.basename(file_path))
        
        def renameLabel(label: str) -> str:
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant'}
            return mapping[label]
        
        self.VERBOSE("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.SKDER_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        skder_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 converters = {'assembly': extractAssembly,
                                               'representative': extractAssembly,
                                               'dereplication_status': renameLabel},
                                 names = ['assembly', 'representative', 'dereplication_status'],
                                 usecols = [0,1,4], header = 0, index_col = 'assembly'
                                 )
        # Join with binary df on Organism column. 
        # Every Organism row is retained (left join).
        # If there is a match between binary_df['Organism'] and skder_df['assembly'] (index), the representative and status is added.
        self.VERBOSE("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(skder_df, on='Organism')
        
        # Extract scaffolds that could not be linked to an assembly:
        scaffolds_with_na = self.binary_df[self.binary_df['representative'].isna()]['Scaffold'].to_list()
        
        if scaffolds_with_na:
            print(f"The following {len(scaffolds_with_na)} scaffolds could not be linked to a representative genome: {', '.join(scaffolds_with_na)}. Omitting for further analysis.")    
            # Drop the NA values:
            self.binary_df = self.binary_df.dropna()
            
        print("Mapping done!")
        
        return None
    
    def run(self):
        """
        Run the entire LocalRun workflow.
        """
        print("\n--- STEP 1: Staging genomes for dereplication. ---")
        self.prepareGenomes()
        
        print("\n--- STEP 2: Dereplicating genomes. ---")
        self.dereplicateGenomes()
        
        print("\n--- STEP 3: Mapping skDER output to binary table. ---")
        self.mapSkderOutToBinary()
        
        print("\n--- STEP 4: Recovering hit diversity. ---")
        self.recoverHits()
        
        print("\n--- STEP 5: Filtering session file. ---")
        self.filterSession()
        
        print("\n--- STEP 6: Generating output files")
        self.generateOutput()
        
        # Remove the temporary directory:
        print("Cleaning up temporary directory.")
        shutil.rmtree(self.TEMP_DIR)

 #####################################################################################################

class RemoteRun(Run):
    
    def __init__(self, args):
        # Call the parent constructor:
        super().__init__(args)
        
        # Defensive check:
        assert args.download_batch > 0, "Download batch should be larger than 0."
        
        # Set path to accessions and download script:
        self.ACCESSIONS_SCRIPT: Path = Path(resources.files(__name__)) / 'get_accessions.sh'  # Path to the script that maps scaffold IDs to assembly IDs
        self.DOWNLOAD_SCRIPT: Path = Path(resources.files(__name__)) / 'download_assemblies.sh'  # Path to the script that download the genomes from given assembly IDs
        
        # Variable to store assembly accessions:
        self.assembly_accessions: list = []
        
        # Dictionary to store the scaffold:assembly mapping:
        self.scaffold_assembly_pairs: dict = {}
        
        # Remove Organisms specified by the user:
        if self.excluded_organisms != {''}:
            # Replace colons with spaces to get the correct matching in the binary table.
            self.excluded_organisms = {org.replace(':', ' ', regex=False) for org in self.excluded_organisms}
            self.VERBOSE(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            # Exclude them:
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]
        
        # Remove scaffolds that the user wants excluded:
        if self.excluded_scaffolds != {''}:
            self.VERBOSE(f"Excluding the following scaffolds: {', '.join(self.excluded_scaffolds)}")
            self.binary_df = self.binary_df[~self.binary_df['Scaffold'].isin(self.excluded_scaffolds)]
        
        # Replace colons in the bypass assemblies as wel:
        if self.bypass_organisms != {''}:
            self.bypass_organisms = {org.replace(':', ' ') for org in self.bypass_organisms}
        
    def fetchAssemblyIDs(self) -> None:
        """
        This function writes the scaffold IDs from the binary table to a file.
        It then calls a bash script that fetches the assembly ID for each scaffold ID using NCBI entrez-direct utilities.
        The results are then written to a file by the bash script, read by this python file, and stored internally as a list of assembly IDs.
        
        Mutates:
            self.assembly_accessions: list: A list of assembly IDs to be downloaded later on.
        """
        # First we extract the scaffold IDs out of the binary table:
        scaffolds = self.binary_df['Scaffold'].to_list()
        
        # Write these to a file with every scaffold ID on a new line:
        scaffolds_file: Path = self.TEMP_DIR / 'scaffolds.txt'
        with scaffolds_file.open('w') as file:
            file.writelines('\n'.join(scaffolds))

        # Now we use the helper bash script to map the scaffolds to the NCBI assembly IDs that host them.
        home = os.getcwd()
        os.chdir(self.TEMP_DIR)
        subprocess.run(['bash', str(self.ACCESSIONS_SCRIPT), 'scaffolds.txt'], check = True)
        os.chdir(home)
        
        # Read the result file
        result_file: Path = self.TEMP_DIR / 'assembly_accessions.txt'
        with result_file.open('r') as file:
            # Store the assembly accessions internally:
            self.assembly_accessions = [line.rstrip() for line in file.readlines()]

        return None
    
    def downloadGenomes(self) -> None:
        """
        This function writes the assembly IDs found by fetchAssemblyIDs to a file in batches (default size of 300). Each line in the file is a batch.
        It then calls a bash script that downloads the genomes for each line of assembly IDs.
        Genome files are placed in the TEMP_DIR/genomes folder (in gzipped format)
        """
        # First we cut off the version digits of our assembly IDs, and rely on NCBI datatsets to fetch the latest version:
        versionless_assemblies = [acc.split('.')[0] for acc in self.assembly_accessions]
        
        self.VERBOSE("Writing assembly ID batches for download script.")
        # Now we write the assembly accesions to a file in batches:
        download_batches_file: Path = self.TEMP_DIR / 'download_batches.txt'
        with download_batches_file.open('w') as file:
            download_batches = list(batched(versionless_assemblies, self.download_batch))
            for batch in download_batches:
                file.write(' '.join(batch) + '\n')
                
        # Run the bash script to download genomes:
        home = os.getcwd()
        os.chdir(self.TEMP_DIR)
        self.VERBOSE("Calling download script.")
        subprocess.run(["bash", str(self.DOWNLOAD_SCRIPT)], check=True)
        os.chdir(home)

        return None
    
    def mapScaffoldsToAssemblies(self) -> None:
        """
        This function maps every scaffold in the binary table to its corresponding assembly file.
        The mapping is stored internally as a dictionary.
        It loops over every file in the genome folder and extracts the set of scaffold using BioPython SeqIO module.
        This set is then compared to the set of scaffolds in the binary table, and the intersection of both sets provied us the correct mapping.
        Scaffolds are stripped off their prefixes because NCBI sometimes omits these when downloading genomes.
        
        Mutates:
            self.scaffold_assembly_pairs: dict: Dictionary mapping of scaffold:assembly_file pairs
        
        """
        def removePrefix(scaffold: str) -> str:
            # ^ matches beginning of sting
            # [^_]+ matches one or more characters that are not an underscore
            # _ matches the underscore itself
            pattern = r'^[^_]+_' 
            # Replace the prefix with an empty string:
            return re.sub(pattern, '', scaffold)
        
        def addPrefix(deprefixed_scaffold: str, scaffolds_in_host_assembly: set) -> str:
            # Map back to the prefixed scaffold by scanning the assembly that hosts the given scaffold.
            return [s for s in scaffolds_in_host_assembly if deprefixed_scaffold in s][0]

        # Define the set of scaffolds (without prefix) that can be found in the binary:
        scaffolds_in_binary: set = {scaffold.strip() for scaffold in self.binary_df['Scaffold'].to_list()}
        # Scaffold set in binary without any prefixes:
        scaffolds_in_binary_no_prefix: set = {removePrefix(scaffold) for scaffold in scaffolds_in_binary}
        
        # Loop over the directory containing all genomes:
        for file in self.GENOME_DIR.iterdir():
            # Only read fasta files:
            if util.isFasta(file.name):
                self.VERBOSE(f"Reading {file.name}")
                # Open the file:
                with gzip.open(file, 'rt') as assembly:
                    # Extract the set of scaffold IDs in the file:
                    scaffolds_in_this_assembly: set = {record.id.strip() for record in SeqIO.parse(assembly, 'fasta')}
                    # Remove the prefixes:
                    scaffolds_in_this_assembly_no_prefix: set = {removePrefix(scaffold) for scaffold in scaffolds_in_this_assembly}
                    # Now we take the intersection of both sets. All the scaffolds in this intersection can be mapped to the current file in the loop:
                    found_scaffolds_no_prefix: set = scaffolds_in_binary_no_prefix.intersection(scaffolds_in_this_assembly_no_prefix)
                    # Now we have to add the prefix again by using the original scaffold list from the host assembly:
                    found_scaffolds: set = {addPrefix(scaffold, scaffolds_in_this_assembly) for scaffold in found_scaffolds_no_prefix}
                    # Tell the user what we found in this assembly file
                    self.VERBOSE(f"Found {','.join(found_scaffolds)}") if len(found_scaffolds) > 0 else self.VERBOSE("No scaffolds found.")
                    # Finalize the mapping in a dictionary:
                    for scaffold in found_scaffolds:
                        self.scaffold_assembly_pairs[scaffold] = file.name
                        
        return None
                
    def mapAssembliesToBinary(self) -> None:
        """
        This function maps each row in the binary table to a corresponding assembly file based on the mapping obtained by mapScaffoldsToAssemblies().
        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
        
        """
        # Read the dictionary mapping as a dataframe with the scaffold IDs as index and the assembly file to which it belongs as 'assembly_file':
        scaffold_assembly_pairs_df: pd.DataFrame = pd.DataFrame.from_dict(self.scaffold_assembly_pairs, orient='index', columns = ['assembly_file'])
        
        # Join the binary table on the 'Scaffold' column:
        self.binary_df = self.binary_df.join(scaffold_assembly_pairs_df, on='Scaffold')
            
        # Extract scaffolds that could not be linked to an assembly:
        scaffolds_with_na = self.binary_df[self.binary_df['assembly_file'].isna()]['Scaffold'].to_list()
        
        if scaffolds_with_na:
            print(f"The following {len(scaffolds_with_na)} scaffolds could not be linked to a genome assembly: {', '.join(scaffolds_with_na)}. Omitting for further analysis.") 
            # Drop the NA values:
            self.binary_df = self.binary_df.dropna()
        
        return None
    
    def mapSkderOutToBinary(self) -> None:
        """
        After dereplicating the genomes, map the skDER clustering table to the binary table.
        skDER clustering table is converted to a df and a join is performed with the binary table based on the assembly_file column.
        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
        """
        def extractFileName(file_path: str) -> str:
            # Extract basename from full file path
            return Path(file_path).name
        
        def renameLabel(label: str) -> str:
            # Rename some of the skDER labels
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant'}
            return mapping[label]
        
        self.VERBOSE("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.SKDER_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        skder_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 converters = {'assembly': extractFileName,
                                               'representative': extractFileName,
                                               'dereplication_status': renameLabel},
                                 names = ['assembly', 'representative', 'dereplication_status'],
                                 usecols = [0,1,4], header = 0, index_col = 'assembly'
                                 )
        # Join with binary df on assembly_file column. 
        # Every assembly_file row is retained (left join).
        # If there is a match between binary_df['assembly_file'] and skder_df['assembly'] (its index column), the representative and status is added.
        self.VERBOSE("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(skder_df, on='assembly_file')

        return None
    
    def run(self) -> None:
        """
        Run the entire remote workflow.
        """
        print("\n--- STEP 1: Fetching assembly IDs from NCBI for each scaffold ID in the cblaster binary table. ---")
        self.fetchAssemblyIDs()  # Stores a list of NCBI assembly IDs 
        
        print("\n--- STEP 2: Downloading genomes for each assembly ID. ---")
        self.downloadGenomes()  # Downloads genome for each assembly ID
        
        print("\n--- STEP 3: Mapping scaffold IDs to assembly IDs ---")
        self.mapScaffoldsToAssemblies()  # Results in a dictionary of scaffold:assembly_file pairs
        self.mapAssembliesToBinary()  # Each row in the binary table is now mapped to its assembly file
        
        print("\n--- STEP 4: Dereplicating genomes ---")
        self.dereplicateGenomes()  # Dereplicate using skDER. Output is in self.SKDER_OUT
        
        print("\n--- STEP 5: Mapping skDER output to binary table ---")
        self.mapSkderOutToBinary()  # Map each row in the binary with its representative genome
        
        print("\n--- STEP 6: Recovering hit diversity ---")
        self.recoverHits()  # Recover hits by content and score depending on user input
        
        print("\n--- STEP 7: Filtering original session file. ---")
        self.filterSession()  # Filter the original session file, retaining only the dereplicated hits.
        
        print("\n--- STEP 8: Generating output files ---")
        self.generateOutput()  # Generate all output files.
        
        # Remove the temporary directory:
        print("Cleaning up temporary directory.")
        shutil.rmtree(self.TEMP_DIR)
        
        return None
        
