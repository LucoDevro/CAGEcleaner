#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.remote_run import RemoteRun
from cagecleaner.genome_run import GenomeRun
from cagecleaner.file_utils import is_fasta
from cagecleaner.communication import get_assembly_accessions
from cagecleaner.utils import run_command

import logging
import os
import re
import gzip
import shutil
import pandas as pd
from pathlib import Path
from itertools import batched
from Bio import SeqIO
from tqdm import tqdm


LOG = logging.getLogger()


class RemoteGenomeRun(RemoteRun, GenomeRun):
    """
    Subclass orchestrating the workflow for dereplication by whole-genome similary using genomes downloaded from NCBI.
    
    This class combines remote genome downloading functionality with local genome dereplication.
    It links NCBI assembly accession IDs to the scaffold IDs in the binary table, downloads the
    corresponding genome assemblies, maps the scaffold IDs back to the downloaded assemblies, performs
    whole-genome dereplication using skDER, and integrates dereplication results back into
    the binary table. Unmapped scaffolds (scaffolds for which no assemby file was found) are tracked
    and reported separately.
    
    Inherits from:
        RemoteRun: Intermediary class providing common downloading configurations.
        GenomeRun: Intermediary class providing genome dereplication utilities.
    """
    
    def __init__(self, args):
        """
        Initialise a RemoteGenomeRun instance.
        
        Creates temporary directories for genome storage, and initialises data structures 
        for tracking assembly accessions and scaffold-to-assembly mappings.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
        
        Raises:
            ValueError: If download_batch is not greater than 0.
        """
        
        super().__init__(args)
        
        # Defensive check:
        if not(args.download_batch > 0): 
            raise ValueError("Download batch should be larger than 0.")
        
        # Variable to store assembly accessions:
        self.assembly_accessions: list = []
        
        # Dictionary to store the scaffold:assembly mapping:
        self.scaffold_assembly_pairs: dict = {}
        
        # Make folder to download genomes to
        self.TEMP_GENOME_DIR = self.TEMP_DIR / "genomes"
        self.TEMP_GENOME_DIR.mkdir(parents = True)
        self.DEREP_IN_DIR = self.TEMP_GENOME_DIR
        
        return None
    
    
    def fetch_assembly_ids(self) -> None:
        """
        Fetch NCBI assembly accession IDs for all scaffold IDs in the binary table.
        
        Extracts scaffold IDs from the binary table and categorises them as RefSeq or GenBank IDs.
        For each category, retrieves corresponding assembly accession IDs using NCBI Entrez Direct utilities.
        Combines and deduplicates the results, storing them in self.assembly_accessions.
        
        Mutates:
            self.assembly_accessions (list): Populated with deduplicated NCBI assembly accession IDs.
            
        Raises:
            RuntimeError: If no assembly IDs have been retrieved.
            
        Returns:
            None
        """
        ## First we extract the scaffold IDs out of the binary table
        scaffolds = self.binary_df['Scaffold'].to_list()
        LOG.info(f"Got {len(scaffolds)} scaffold IDs to crossref")
        
        ## Split up in RefSeq and Genbank IDs
        refseq_scaffolds = [sc for sc in scaffolds if sc.startswith('NC_') or sc.startswith('NZ_')]
        genbank_scaffolds = [sc for sc in scaffolds if not(sc.startswith('NC_') or sc.startswith('NZ_'))]
        
        LOG.info(f'..., of which {len(refseq_scaffolds)} are RefSeq IDs')
        LOG.info(f'..., and {len(genbank_scaffolds)} are Genbank IDs.')
        
        ## Fetch RefSeq assembly IDs
        LOG.info('Fetching RefSeq Assembly IDs')
        refseq_assembly_accessions = get_assembly_accessions(refseq_scaffolds, 'RefSeq', no_progress = self.no_progress)
        
        ## Fetch Genbank assembly IDs
        LOG.info('Fetching Genbank Assembly IDs')
        genbank_assembly_accessions = get_assembly_accessions(genbank_scaffolds, 'Genbank', no_progress = self.no_progress)
        
        ## Gather and deduplicate assembly IDs
        LOG.info('Merging ID sets')
        assembly_accessions = refseq_assembly_accessions + genbank_assembly_accessions
        
        ## Stop if empty
        if len(assembly_accessions) == 0:
            msg = 'No assembly IDs retrieved!'
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        LOG.info(f'Got {len(assembly_accessions)} accession IDs in total')
        self.assembly_accessions = assembly_accessions

        return None
    
    
    def fetch_genomes(self) -> None:
        """
        Fetch genome assemblies from NCBI in batches using the datasets command-line tool.
        
        Downloads assemblies in batches (default 300 per batch) using the NCBI Datasets CLI. For each batch,
        fetches genome data using NCBI datasets, rehydrates the dehydrated files with gzip compression,
        and moves all resulting genome files to the temporary genome directory. Removes version digits
        from assembly IDs to ensure the latest versions are downloaded.
        
        Mutates:
            Populates self.TEMP_GENOME_DIR with downloaded gzipped genome files.
        
        Returns:
            None
        """
        LOG.info('Preparing to download assembies...')
        
        # First we cut off the version digits of our assembly IDs, and rely on NCBI datatsets to fetch the latest version:
        versionless_assemblies = [acc.split('.')[0] for acc in self.assembly_accessions]
        
        # Make download batches
        download_batches = list(batched(versionless_assemblies, self.download_batch))
        
        # Then download all batches
        nb_assemblies = len(versionless_assemblies)
        nb_batches = len(download_batches)
        
        LOG.info(f'Got {nb_assemblies} assembly IDs. Downloading genomes in {nb_batches} batches.')
        
        for idx, batch in tqdm(list(enumerate(download_batches)), leave = False, disable = self.no_progress):
            LOG.info(f"Downloading batch {idx+1} out of {nb_batches}")
            DOWNLOADS = self.TEMP_DIR / 'downloads'
            DOWNLOADS.mkdir(parents = True, exist_ok = True)
            datasets_zip = DOWNLOADS / "ncbi_dataset.zip"
            
            # download command
            fetch_cmd = ['datasets', 'download', 'genome', 'accession',
                         ','.join(batch),
                         '--filename', str(datasets_zip.resolve()),
                         '--dehydrated',
                         '--no-progressbar']
            
            try:
                run_command(fetch_cmd)
            except RuntimeError:
                msg = 'Fetching download URLs from NCBI failed!'
                LOG.critical(msg)
                raise RuntimeError(msg)
            
            # Unzip the package
            shutil.unpack_archive(datasets_zip, extract_dir = DOWNLOADS, format = "zip")
            os.remove(datasets_zip)
                
            # Do the actual download
            download_cmd = ['datasets', 'rehydrate',
                            '--directory', str(DOWNLOADS.resolve()),
                            '--gzip',
                            '--no-progressbar',
                            '--max-workers', str(self.download_workers)]
            
            try:
                run_command(download_cmd)
            except RuntimeError:
                msg = 'Downloading genome assemblies from NCBI failed!'
                LOG.critical(msg)
                raise RuntimeError(msg)
            
            # Move all the genome files to the genome folder
            for file in DOWNLOADS.glob('ncbi_dataset/data/GC*/*'):
                shutil.move(file, self.TEMP_GENOME_DIR)
            
            # Remove the downloads folder
            shutil.rmtree(DOWNLOADS)
            
        LOG.info(f'{len(list(self.TEMP_GENOME_DIR.iterdir()))} genomes in {str(self.TEMP_GENOME_DIR)}')
            
        return None
    
    
    def map_scaffolds_to_assemblies(self) -> None:
        """
        Map each scaffold ID from the binary table to its corresponding downloaded assembly file.
        
        Iterates through all genome files in the temporary genome directory and extracts scaffold IDs
        using BioPython's SeqIO. For each assembly file, compares its scaffolds (with and without prefixes)
        to those in the binary table. When matches are found, stores the mapping in self.scaffold_assembly_pairs.
        Prefixes are stripped during comparison because NCBI sometimes omits them in downloaded genomes.
        
        Mutates:
            self.scaffold_assembly_pairs (dict): Populated with (scaffold_id: assembly_filename) pairs.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If the dereplication input folder does not contain any fasta file.        
        """
        def remove_prefix(scaffold: str) -> str: 
            """
            Remove the prefix from a scaffold ID by removing all characters up to and including the first underscore.
            
            Args:
                scaffold (str): A scaffold ID string.
            
            Returns:
                str: The scaffold ID with its prefix removed.
            """
            return re.sub(r'^[^_]+_', '', scaffold)
        
        def find_prefix(deprefixed_scaffold: str, scaffolds_in_host_assembly: set) -> str:
            """
            Find and return the original prefixed scaffold ID from the set of scaffolds in the local assembly file.
            
            Args:
                deprefixed_scaffold (str): A scaffold ID without its prefix.
                scaffolds_in_host_assembly (set): Set of all scaffold IDs in the assembly, with original prefixes intact.
            
            Returns:
                str: The prefixed scaffold ID that matches the deprefixed version.
            """
            return [s for s in scaffolds_in_host_assembly if deprefixed_scaffold in s][0]

        # Define the set of scaffolds (without prefix) that can be found in the binary:
        scaffolds_in_binary = {scaffold.strip() for scaffold in self.binary_df['Scaffold'].to_list()}
        # Scaffold set in binary without any prefixes:
        scaffolds_in_binary_no_prefix = {remove_prefix(scaffold) for scaffold in scaffolds_in_binary}
        
        # Verify that the dereplication input directory is not empty
        try:
            next(filter(is_fasta, self.DEREP_IN_DIR.iterdir()))
        except StopIteration:
            msg = "No fasta files found in the dereplication input folder!"
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        # Loop over all fasta files in the genome directory:
        for file in filter(is_fasta, self.DEREP_IN_DIR.iterdir()):
            LOG.debug(f"Reading {file.name}")
            # Open the file:
            with gzip.open(file, 'rt') as assembly:
                # Extract the set of scaffold IDs in the file and remove prefixes
                scaffolds_in_this_assembly = {record.id.strip() for record in SeqIO.parse(assembly, 'fasta')}
                # Remove the prefixes
                scaffolds_in_this_assembly_no_prefix = {remove_prefix(record.id.strip()) for record in SeqIO.parse(assembly, 'fasta')}
                # Now we take the intersection of both sets. All the scaffolds in this intersection are found in the the current file
                found_scaffolds_no_prefix = scaffolds_in_binary_no_prefix.intersection(scaffolds_in_this_assembly_no_prefix)
                # Now we have to add the prefix again by using the original scaffold list from the host assembly:
                found_scaffolds = {find_prefix(scaffold, scaffolds_in_this_assembly) for scaffold in found_scaffolds_no_prefix}
                # Tell the user what we found in this assembly file
                LOG.debug(f"Found {','.join(found_scaffolds)}") if len(found_scaffolds) > 0 else LOG.debug("No scaffolds found.")
                # Finalize the mapping in a dictionary:
                for scaffold in found_scaffolds:
                    self.scaffold_assembly_pairs[scaffold] = file.name
                        
        return None
    
                
    def join_assemblies_with_binary(self) -> None:
        """
        This function maps each row in the binary table to a corresponding assembly file based on the mapping obtained by map_scaffolds_to_assemblies().
        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
            
        Join the scaffold-to-assembly mapping with the binary table and remove unmapped scaffolds.
        
        Converts self.scaffold_assembly_pairs dictionary to a DataFrame and joins with self.binary_df 
        on the Scaffold column. Identifies scaffolds that could not be linked to any assembly file, 
        logs a warning, saves them to an output file, and removes them from the binary table.
        
        Mutates:
            self.binary_df (pd.DataFrame): Adds 'assembly_file' column and removes rows with unmatched scaffolds.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If the binary table is empty after joining with the mapping table.
        """
        # Read the dictionary mapping as a dataframe with the scaffold IDs as index and the assembly file to which it belongs as 'assembly_file':
        scaffold_assembly_pairs_df: pd.DataFrame = pd.DataFrame.from_dict(self.scaffold_assembly_pairs, orient='index', columns = ['assembly_file'])
        
        # Join the binary table on the 'Scaffold' column:
        self.binary_df = self.binary_df.join(scaffold_assembly_pairs_df, on='Scaffold')
        
        # Extract scaffolds that could not be linked to an assembly:
        scaffolds_with_na = self.binary_df[self.binary_df['assembly_file'].isna()]['Scaffold']
        
        if not(scaffolds_with_na.empty):
            scaffolds_with_na.to_csv(self.OUT_DIR / "unmapped.scaffolds.txt", index = False, header = False)
            LOG.warning(f"{scaffolds_with_na.size} scaffolds could not be linked to a genome assembly. See unmapped.scaffolds.txt") 
            # Drop the NA values:
            self.binary_df = self.binary_df.dropna()
            
        if self.binary_df.empty:
            msg = "No scaffold could be linked with an assembly!"
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        return None
    
    
    def join_dereplication_with_binary(self) -> None:
        """
        Join the genome dereplication clustering results with the binary table.
        
        Reads the skDER clustering output file and converts it to a DataFrame with assembly filenames
        and dereplication status. Joins this data with self.binary_df on the assembly_file column.
        Sorts the resulting table by representative genome and dereplication status.
        
        Mutates:
            self.binary_df (pd.DataFrame): Adds 'representative' and 'dereplication_status' columns and sorts by these values.
        
        Returns:
            None
            
        Raises:
            FileNotFoundError: If the dereplication table cannot be read
            RuntimeError: If the dereplication table is empty.
            RuntimeError: If the binary table is empty after joining with the dereplication table.
            
        Note:
            This is the workflow-specific implementation of the abstract method inherited from its grandparent class Run.
        """
        
        def extract_filename(file_path: str) -> str: return Path(file_path).name
        
        def rename_skder_label(label: str) -> str:
            # Rename some of the skDER labels
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant',
                       'outside_cutoffs_requested': 'redundant'} # edge case that clusters by skani dist, but fails the within-cluster cutoffs
            return mapping[label]
        
        LOG.debug("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        try:
            derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                     converters = {'assembly': extract_filename,
                                                   'representative': extract_filename,
                                                   'dereplication_status': rename_skder_label},
                                     names = ['assembly', 'representative', 'dereplication_status'],
                                     usecols = [0,1,4], header = 0, index_col = 'assembly'
                                     )
        except FileNotFoundError as err:
            LOG.critical(f'{err}')
            raise err
        if derep_df.empty:
            msg = "Dereplication table is empty!"
            LOG.error(msg)
            raise RuntimeError(msg)
            
        # Join with binary df on assembly_file column. 
        # Every assembly_file row is retained (left join).
        # If there is a match between binary_df['assembly_file'] and derep_df['assembly'] (its index column), the representative and status is added.
        LOG.debug("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(derep_df, on='assembly_file')
        
        # Verify binary table is not empty
        if self.binary_df.empty:
            msg = "Binary table is empty after joining with the dereplication table!"
            LOG.error(msg)
            raise RuntimeError(msg)
        
        # Sort by representative ID and then by dereplication status
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])
        LOG.info("Mapping done!")
            
        return None
    
    
    def run(self):
        """
        Execute the complete remote genome dereplication pipeline.
        
        Orchestrates all processing steps in sequence: fetches NCBI assembly IDs for scaffold IDs,
        downloads genomes, maps scaffold IDs to assembly IDs, performs whole-genome dereplication,
        joins dereplication results with the binary table, recovers hit diversity according to
        user parameters, filters the original session file, and generates output files.
        Cleans up the temporary directory upon completion.
        
        Returns:
            None
        """
        
        LOG.info("--- STEP 1: Fetching NCBI assembly IDs for each scaffold ID. ---")
        self.fetch_assembly_ids()  # Stores a list of NCBI assembly IDs 
    
        LOG.info("--- STEP 2: Downloading genomes for each assembly ID. ---")
        self.fetch_genomes()  # Downloads genome for each assembly ID
        
        LOG.info("--- STEP 3: Mapping scaffold IDs to assembly IDs ---")
        self.map_scaffolds_to_assemblies()  # Results in a dictionary of scaffold:assembly_file pairs
        self.join_assemblies_with_binary()  # Each row in the binary table is now mapped to its assembly file
        
        LOG.info("--- STEP 4: Dereplicating ---")
        self.dereplicate_genomes()  # Dereplicate.
        
        LOG.info("--- STEP 5: Mapping dereplication clustering to binary table ---")
        self.join_dereplication_with_binary()  # Map each row in the binary table with its representative
        
        LOG.info("--- STEP 6: Recovering hit diversity ---")
        self.recover_hits()  # Recover hits by content and score depending on user input
        
        LOG.info("--- STEP 7: Filtering original session file. ---")
        self.filter_session()  # Filter the original session file, retaining only the dereplicated hits.
        
        LOG.info("--- STEP 8: Generating output files ---")
        self.generate_output()  # Generate all output files.
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
    