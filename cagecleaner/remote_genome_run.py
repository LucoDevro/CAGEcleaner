#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.remote_run import RemoteRun
from cagecleaner.genome_run import GenomeRun
from cagecleaner import util

import logging
import os
import subprocess
import re
import gzip
import sys
import pandas as pd
from pathlib import Path
from itertools import batched
from Bio import SeqIO
from importlib import resources


LOG = logging.getLogger()


class RemoteGenomeRun(RemoteRun, GenomeRun):
    
    def __init__(self, args):
        
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
        
        # Make folder to download genomes to
        self.TEMP_GENOME_DIR = self.TEMP_DIR / "genomes"
        self.TEMP_GENOME_DIR.mkdir(parents = True)
        
        return None
    
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
        
        LOG.debug("Writing assembly ID batches for download script.")
        # Now we write the assembly accesions to a file in batches:
        download_batches_file: Path = self.TEMP_DIR / 'download_batches.txt'
        with download_batches_file.open('w') as file:
            download_batches = list(batched(versionless_assemblies, self.download_batch))
            for batch in download_batches:
                file.write(' '.join(batch) + '\n')
                
        # Run the bash script to download genomes:
        home = os.getcwd()
        os.chdir(self.TEMP_DIR)
        LOG.debug("Calling download script.")
        subprocess.run(["bash", str(self.DOWNLOAD_SCRIPT), str(self.download_workers)], check=True)
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
        for file in self.DEREP_IN_DIR.iterdir():
            # Only read fasta files:
            if util.isFasta(file.name):
                LOG.debug(f"Reading {file.name}")
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
                    LOG.debug(f"Found {','.join(found_scaffolds)}") if len(found_scaffolds) > 0 else LOG.debug("No scaffolds found.")
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
            with open(self.OUT_DIR / "unmapped.scaffolds.txt", "w") as handle:
                handle.writelines(scaffolds_with_na)
            LOG.warning(f"{len(scaffolds_with_na)} scaffolds could not be linked to a genome assembly. See unmapped.scaffolds.txt") 
            # Drop the NA values:
            self.binary_df = self.binary_df.dropna()
        
        return None
    
    def mapDereplicationToBinary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).        
        Mutates:
            self.binary_df: pd.DataFrame: Internal representation of the binary table.
        """
        # Full genome dereplication using skDER
        def extractFileName(file_path: str) -> str:
            # Extract basename from full file path
            return Path(file_path).name
        
        def renameLabel(label: str) -> str:
            # Rename some of the skDER labels
            mapping = {'representative_to_self': 'dereplication_representative',
                       'within_cutoffs_requested': 'redundant',
                       'outside_cutoffs_requested': 'redundant'} # edge case that clusters by skani dist, but fails the clustering cutoffs
            return mapping[label]
        
        LOG.debug("Reading skDER clustering table.")
        # Read the skder out clustering table:
        path_to_cluster_file: Path = self.DEREP_OUT_DIR / 'skDER_Clustering.txt'
        # Convert to dataframe:
        derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                 converters = {'assembly': extractFileName,
                                               'representative': extractFileName,
                                               'dereplication_status': renameLabel},
                                 names = ['assembly', 'representative', 'dereplication_status'],
                                 usecols = [0,1,4], header = 0, index_col = 'assembly'
                                 )
        # Join with binary df on assembly_file column. 
        # Every assembly_file row is retained (left join).
        # If there is a match between binary_df['assembly_file'] and derep_df['assembly'] (its index column), the representative and status is added.
        LOG.debug("Joining skDER clustering table and cblaster binary table.")
        self.binary_df = self.binary_df.join(derep_df, on='assembly_file')
        
        # Sort by representative ID and then by dereplication status
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])
        LOG.info("Mapping done!")
            
        return None
    
    def run(self):
        
        LOG.info("--- STEP 1: Fetching NCBI assembly IDs for each scaffold ID. ---")
        self.fetchAssemblyIDs()  # Stores a list of NCBI assembly IDs 
    
        LOG.info("--- STEP 2: Downloading genomes for each assembly ID. ---")
        self.downloadGenomes()  # Downloads genome for each assembly ID
        
        LOG.info("--- STEP 3: Mapping scaffold IDs to assembly IDs ---")
        self.mapScaffoldsToAssemblies()  # Results in a dictionary of scaffold:assembly_file pairs
        self.mapAssembliesToBinary()  # Each row in the binary table is now mapped to its assembly file
        
        LOG.info("--- STEP 4: Dereplicating ---")
        self.dereplicateGenomes()  # Dereplicate.
        
        LOG.info("--- STEP 5: Mapping dereplication clustering to binary table ---")
        self.mapDereplicationToBinary()  # Map each row in the binary table with its representative
        
        LOG.info("--- STEP 6: Recovering hit diversity ---")
        self.recoverHits()  # Recover hits by content and score depending on user input
        
        LOG.info("--- STEP 7: Filtering original session file. ---")
        self.filterSession()  # Filter the original session file, retaining only the dereplicated hits.
        
        LOG.info("--- STEP 8: Generating output files ---")
        self.generateOutput()  # Generate all output files.
        
        # Remove the temporary directory:
        LOG.info("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        return None
    
    