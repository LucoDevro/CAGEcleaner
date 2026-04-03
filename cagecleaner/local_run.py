#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.file_utils import isFasta, isGff, isGenbank, removeSuffixes, convertGenbankToFasta

import logging
import os
import sys
from abc import abstractmethod

LOG = logging.getLogger()


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
            assert (isFasta(str(file)) or isGff(str(file)) or isGenbank(str(file))), f"The following file is not in the correct format: {file}.\nWe only accept the following suffices: [.fasta, .fna, .fa, .gff3, .gff, .gbff, .gbk, .gb]. The '.gz' extension is allowed."
        
        # Remove organisms that the user wants to be excluded:
        if self.excluded_organisms != {''}:
            LOG.debug(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]

        # Remove scaffold IDs specified by the user:
        if self.excluded_scaffolds != {''}:
            LOG.debug(f"Excluding the following scaffolds: {', '.join(self.excluded_scaffolds)}")
            # Here the approach slightly differs as users might have provided prefixed scaffold IDs:
            # Add a column with a prefixed scaffold based on the Organism column
            self.binary_df['prefixed_scaffold'] = self.binary_df['Organism'] + ':' + self.binary_df['Scaffold']
            # If the prefixed scaffold ends with any of the strings in the set of scaffolds to exclude, remove it:
            self.binary_df = self.binary_df[self.binary_df['prefixed_scaffold'].str.endswith(tuple(self.excluded_scaffolds)) == False]
            # Clean up:
            self.binary_df = self.binary_df.drop(columns=['prefixed_scaffold'])
            
        return None
    
    @abstractmethod
    def mapDereplicationToBinary(self):
        pass
            
    def prepareGenomes(self) -> None:
        """
        This function is called to inspect the provided genome folder and convert GenBank files to FASTA files if necessary.
        If a FASTA file is found, it is assumed that all genomes are stored in the provided folder in FASTA format.
        If no FASTA is found, but instead a GenBank is found, it is assumed that all genomes are in GenBank format.
        Temporary FASTA-converted copies of these GenBank files will be placed in the temporary folders.
        Only FASTA files can be passed on to skDER.
        If neither GenBank nor FASTA is found, the program exits.
    
        Mutates:
            self.DEREP_IN_DIR: Path: Path to the genome directory
        """
        
        # Assert that Organism and filenames correspond:
        files_in_genomes_dir: set = {removeSuffixes(file) for file in os.listdir(self.USER_GENOME_DIR)}
        organisms_in_session: set = {removeSuffixes(organism.name) for organism in self.session.organisms}
        
        assert files_in_genomes_dir >= organisms_in_session, "The genomes of all organisms in the cblaster session have not been found in the genome directory. Check the paths and please make sure you have not changed the genome filenames between a cblaster run and a CAGEcleaner run."
        
        # Check if there are FASTA files in the genome folder:
        fasta_in_folder = [isFasta(str(file)) for file in self.USER_GENOME_DIR.iterdir()]
        genbank_in_folder = [isGenbank(str(file)) for file in self.USER_GENOME_DIR.iterdir()]

        if any(fasta_in_folder):
            # In this case the genome folder path should remain the same
            LOG.info(f"Detected {sum(fasta_in_folder)} FASTA files in {self.USER_GENOME_DIR}. These will be used for dereplication.")
            # Redirect the genome dir to the user-provided folder:
            self.TEMP_GENOME_DIR = self.USER_GENOME_DIR
            
        elif any(genbank_in_folder):
            # In this case we convert to FASTA and redirect to genome folder, which is in the temp folder by default.
            LOG.info(f"Detected {sum(genbank_in_folder)} GenBank files in {self.USER_GENOME_DIR}.")
            # Convert to FASTA files:
            self.TEMP_GENOME_DIR = self.TEMP_DIR / 'genomes'
            self.TEMP_GENOME_DIR.mkdir(exist_ok=True)  # Make the temporary genome folder if it does not exist already.
            convertGenbankToFasta(self.USER_GENOME_DIR, self.TEMP_GENOME_DIR, workers = self.cores)
            LOG.info(f"Saved genomes in FASTA format to {self.TEMP_GENOME_DIR}")
            
        else:
            # If there are no FASTA or GenBank files, the program cannot proceed:
            LOG.critical("No FASTA files or GenBank files were detected in the provided genome folder. Exiting the program.")
            sys.exit()
            
        assembly_files = [acc + ".fasta.gz" for acc in self.binary_df['Organism']]
        self.binary_df['assembly_file'] = assembly_files
        
        return None
    
    