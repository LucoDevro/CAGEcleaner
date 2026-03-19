#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Internal imports:
from cagecleaner import util
from cagecleaner.classes import Run

# External libraries:
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

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
        
        assert files_in_genomes_dir >= organisms_in_session, "The genomes of all organisms in the cblaster session have not been found in the genome directory. Check the paths and please make sure you have not changed the genome filenames between a cblaster run and a CAGEcleaner run."
        
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
            
        assembly_files = [acc + ".fasta.gz" for acc in self.binary_df['Organism']]
        self.binary_df['assembly_file'] = assembly_files
        
        return None
    
    def mapDereplicationToBinary(self) -> None:
        """
        After dereplication, map the dereplication clustering table to the binary table.
        The dereplication clustering table is converted to a dataframe and joined with the binary table based on
        assembly ID (full genome dereplication) or scaffold ID (region dereplication).  
        
        Mutates:
            self.binary_df: pd.DataFrame: The binary table derived from a cblaster Session object.
        """
        # Full genome dereplication using skDER
        if not(self.regions):
            def extractAssembly(file_path: str) -> str:
                return util.removeSuffixes(os.path.basename(file_path))
            
            def renameLabel(label: str) -> str:
                mapping = {'representative_to_self': 'dereplication_representative',
                           'within_cutoffs_requested': 'redundant'}
                return mapping[label]
            
            self.VERBOSE("Reading skDER clustering table.")
            # Read the skder out clustering table:
            path_to_cluster_file: Path = self.DEREP_OUT_DIR / 'skDER_Clustering.txt'
            # Convert to dataframe:
            derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                     converters = {'assembly': extractAssembly,
                                                   'representative': extractAssembly,
                                                   'dereplication_status': renameLabel},
                                     names = ['assembly', 'representative', 'dereplication_status'],
                                     usecols = [0,1,4], header = 0, index_col = 'assembly'
                                     )
            # Join with binary df on Organism column. 
            # Every Organism row is retained (left join).
            # If there is a match between binary_df['Organism'] and derep_df['assembly'] (index), the representative and status is added.
            self.VERBOSE("Joining skDER clustering table and cblaster binary table.")
            self.binary_df = self.binary_df.join(derep_df, on='Organism')
        
        # Region dereplication using MMseqs2
        else:
            # Read the MMseqs clustering table
            path_to_cluster_file: Path = self.TEMP_DIR / "derep_out_cluster.tsv"
            # Convert to a dataframe
            derep_df: pd.DataFrame = pd.read_table(path_to_cluster_file,
                                     names = ['representative', 'Scaffold'],
                                     header = None, index_col = 'Scaffold'
                                     )
            # Add dereplication status column based on whether a scaffold's ID is the same as the representative's one
            # Every assembly_file row is retained (left join).
            # If there is a match between binary_df['assembly_file'] and derep_df['assembly'] (its index column), the representative and status is added.
            derep_df['dereplication_status'] = derep_df.index == derep_df['representative']
            derep_df['dereplication_status'] = np.where(derep_df['dereplication_status'], 'dereplication_representative', 'redundant')
            self.VERBOSE("Joining MMseqs2 clustering table and cblaster binary table.")
            self.binary_df = self.binary_df.join(derep_df, on = "Scaffold")
            
        self.binary_df = self.binary_df.sort_values(['representative', 'dereplication_status'])           
        print("Mapping done!")
        
        return None
    
    def run(self):
        """
        Run the entire LocalRun workflow.
        """
        print("\n--- STEP 1: Staging genomes for dereplication. ---")
        self.prepareGenomes()
        
        print("\n--- STEP 2: Dereplicating. ---")
        self.dereplicate()
        
        print("\n--- STEP 3: Mapping dereplication output to binary table. ---")
        self.mapDereplicationToBinary()
        
        print("\n--- STEP 4: Recovering hit diversity. ---")
        self.recoverHits()
        
        print("\n--- STEP 5: Filtering session file. ---")
        self.filterSession()
        
        print("\n--- STEP 6: Generating output files")
        self.generateOutput()
        
        # Remove the temporary directory:
        print("Cleaning up temporary directory.")
        self.TEMP_DIR_CONTEXT.cleanup()
        
        