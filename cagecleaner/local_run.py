#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.file_utils import is_fasta, is_genbank, remove_suffixes, convert_genbanks_to_fastas

import logging
from abc import abstractmethod

LOG = logging.getLogger()


class LocalRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving local sequence files.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         LocalRegionRun: Region-based dereplication for hits in local sequences.
         LocalGenomeRun: Whole-genome dereplication for hits in local sequences.
    """
    
    def __init__(self, args):
        """
        Initialise a LocalRun instance.
        
        Runs the base class init and checks for valid keep flags.
        Checks whether genome folder only contains valid genome sequence files.
        Excludes scaffolds from the analysis as specified by the user.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
            
        Returns:
            None
            
        Raises:
            ValueError: If the download flag has been enabled in local mode.
            ValueError: If the user-supplied local genome directory does not contain any fasta or genbank file.
            ValueError: If there are no hits left in the binary table after excluding scaffolds or organisms.
        """
                
        # Call the parent class initiator
        super().__init__(args)
        
        # Can't keep downloads in local mode. Also make sure keep_intermediate doesn't fail.
        if self.keep_downloads == True:
            raise ValueError("Can't keep downloads in local mode.")
        if self.keep_intermediate:
            self.keep_intermediate = False
            self.keep_dereplication = True
        
        # Make sure there is no exotic stuff in the provided genome folder
        try:
            next(filter(lambda x: is_fasta(x) or is_genbank(x), self.USER_GENOME_DIR.iterdir()))
        except StopIteration:
            msg = "The user-supplied genome directory does not contain any fasta or genbank file!"
            LOG.critical(msg)
            raise ValueError(msg)
        
        # Remove organisms that the user wants to be excluded:
        if self.excluded_organisms != {''}:
            LOG.debug(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]
            
            if self.binary_df.empty:
                msg = "No hits left after excluding organisms!"
                LOG.error(msg)
                raise ValueError(msg)

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
            
            if self.binary_df.empty:
                msg = "No hits left after excluding scaffolds!"
                LOG.error(msg)
                raise ValueError(msg)
            
        return None
    
    
    @abstractmethod
    def join_dereplication_with_binary(self):
        """
        Join dereplication results with the binary table.
        
        Mutates:
            self.binary_df: Adds 'representative' and 'dereplication_status' columns.
            
         Expected Result:
             self.binary_df should now have columns:
             - representative: Genome ID of the dereplication representative
             - dereplication_status: 'dereplication_representative' | 'redundant'
             
         Returns:
             None
             
        Notes:
            This is the abstract method inherited from the Run parent class.
            It is not meant to be implemented at this level. Only the child classes inheriting this method
            are expected to provide a workflow-dependent specific implementation.
        """
        pass
    
            
    def prepare_genomes(self) -> None:
        """
        Prepare the genome sequence files in the specified genome directory for dereplication.
        
        Checks whether The filenames of the genome sequence files are among the names of the organisms in the Session object,
        ignoring file extensions. Checks whether there are fasta and genbank files in the user-specified genome folder,
        converting genbank files to fasta files on-the-fly.
        
        Adds a column assembly_file to binary table specifying the filepath of each scaffold's associated genome assembly.
        In Genbank mode, this will point to converted files in the temporary directories.
        
        Mutates:
            self.binary_df (pd.DataFrame): Updated in-place with an additional column for
                'assembly_file' and 'dereplication_status'.
                
        Returns:
            None
            
        Raises:
            ValueError: If an organism is found of which the genome is not present in the user-supplied genome directory.
            RuntimeError: If no fasta or genbank files have been found in the supplied genome directory.
            
        Notes:
            The sequence files in the user genome folder should be either all fasta files or all genbank files. There is
            no mix case support.
        """
        
        # Assert that Organism and filenames correspond:
        files_in_genomes_dir = {remove_suffixes(file.name) for file in self.USER_GENOME_DIR.iterdir()}
        organisms_in_session = {remove_suffixes(organism.name) for organism in self.session.organisms}
        if not(files_in_genomes_dir >= organisms_in_session): 
            raise ValueError("Not all genomes of the organisms in the session have been found in the genome directory. Make sure you have not changed the genome filenames between a cblaster run and a CAGEcleaner run.")
        
        # Check if there are valid sequence files in the genome folder:
        try:
            next(filter(lambda x: is_fasta(x) or is_genbank(x), self.USER_GENOME_DIR.iterdir()))
        except StopIteration:
            msg = "No fasta files or Genbank files were found in the provided genome folder!"
            LOG.critical(msg)
            raise RuntimeError(msg)
            
        fasta_in_folder = [is_fasta(str(file)) for file in self.USER_GENOME_DIR.iterdir()]
        genbank_in_folder = [is_genbank(str(file)) for file in self.USER_GENOME_DIR.iterdir()]

        if any(filter(is_fasta, self.USER_GENOME_DIR.iterdir())):
            # In this case the genome folder path should remain the same
            LOG.info(f"Detected {sum(fasta_in_folder)} FASTA files in {self.USER_GENOME_DIR}. These will be used for dereplication.")
            # Redirect the genome dir to the user-provided folder:
            self.TEMP_GENOME_DIR = self.USER_GENOME_DIR
            
        elif any(filter(is_genbank, self.USER_GENOME_DIR.iterdir())):
            # In this case we convert to FASTA and redirect to genome folder, which is in the temp folder by default.
            LOG.info(f"Detected {sum(genbank_in_folder)} GenBank files in {self.USER_GENOME_DIR}.")
            # Convert to FASTA files:
            self.TEMP_GENOME_DIR = self.TEMP_DIR / 'genomes'
            self.TEMP_GENOME_DIR.mkdir(exist_ok = True)  # Make the temporary genome folder if it does not exist already.
            convert_genbanks_to_fastas(self.USER_GENOME_DIR, self.TEMP_GENOME_DIR, workers = self.cores)
            LOG.info(f"Saved genomes in FASTA format to {self.TEMP_GENOME_DIR}")
            
        assembly_files = [acc + ".fasta.gz" for acc in self.binary_df['Organism']]
        self.binary_df['assembly_file'] = assembly_files
        
        return None
    
    