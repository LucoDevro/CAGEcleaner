#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.file_utils import is_fasta, is_genbank, remove_suffixes, convert_genbanks_to_fastas

import pandas as pd
import logging
import shutil
from abc import abstractmethod
from cblaster.extract_clusters import get_sorted_cluster_hierarchies

LOG = logging.getLogger(__name__)


class LocalRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving local sequence files.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         LocalRegionRun: Region-based dereplication for hits in local sequences.
         LocalGenomeRun: Whole-genome dereplication for hits in local sequences.
    """
    
    def __init__(self, parsed_args):
        """
        Initialise a LocalRun instance.
        
        Runs the base class init and excludes scaffolds from the analysis as specified by the user.
        
        Args:
            parsed_args (dict): Parsed and validated command-line arguments
            
        Returns:
            None
            
        Raises:
            RuntimeError: If there are no hits left in the binary table after excluding scaffolds or organisms.
        """
                
        # Call the parent class initiator
        super().__init__(parsed_args)
        
        # Automatically toggle off keeping both downloads and dereplication in local mode,
        # since there won't be downloaded anything.
        if self.keep_intermediate:
            self.keep_intermediate = False
            self.keep_downloads = False
            self.keep_dereplication = True
        
        # Remove organisms that the user wants to be excluded:
        if self.excluded_organisms != {''}:
            LOG.debug(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]
            
            if self.binary_df.empty:
                msg = "No hits left after excluding organisms!"
                LOG.error(msg)
                raise RuntimeError(msg)

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
                raise RuntimeError(msg)
            
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
        Prepare the genome sequence files in the specified genome directory for dereplication, ignoring genomes without hit.
        
        Checks whether The filelabels of the genome sequence files are among the names of the organisms in the Session object,
        ignoring file extensions. Selects only genome files with a hit, and checks whether these are fasta and genbank files.
        Converts any genbank file to fasta on-the-fly.
        
        Adds a column assembly_file to binary table specifying the filepath of each scaffold's associated genome assembly.
        For Genbank files, this will be the filename of the converted fasta file.
        
        Mutates:
            self.binary_df (pd.DataFrame): Updated in-place with an additional column for
                'assembly_file'.
                
        Returns:
            None
            
        Raises:
            ValueError: If an organism is found of which the genome is not present in the user-supplied genome directory.
            RuntimeError: If no fasta or genbank files have been found in the supplied genome directory.
        """
        ## Get an overview of the genome files in the genome directory, and the hit-containing genomes included in the session
        # Make it a mapping from filename to filelabel
        files_in_genomes_dir = [file.name for file in self.USER_GENOME_DIR.iterdir()]
        assemblies_in_session = [hit[2] for hit in get_sorted_cluster_hierarchies(self.session, max_clusters = None)]

        ## Assert that all hit-coding assemblies included in the session have an assembly file
        files_in_genomes_dir_no_suffix = [remove_suffixes(file) for file in files_in_genomes_dir]
        assemblies_in_session_no_suffix = [remove_suffixes(label) for label in assemblies_in_session]
        if not(set(files_in_genomes_dir_no_suffix) >= set(assemblies_in_session_no_suffix)): 
            raise ValueError("Not all relevant genome files have been found in the genome directory. Make sure you have not changed the genome filenames between a cblaster run and a CAGEcleaner run.")
        
        ## Make a mapping from filelabel to actual path, discarding assembly files without a hit
        relevant_paths = set([[path for path in files_in_genomes_dir if label in str(path)][0]
                              for label in assemblies_in_session])
        
        ## Gather all genomes file, checking whether the files have valid file extensions
        # Make the temporary genome folder if it does not exist already
        self.TEMP_GENOME_DIR = self.TEMP_DIR / 'genomes'
        self.TEMP_GENOME_DIR.mkdir(exist_ok = True)
        # Same for the Genbank temporary subfolder
        genbanks_temp_subfolder = self.TEMP_GENOME_DIR / 'genbanks'
        genbanks_temp_subfolder.mkdir()
        
        fastas_found = 0
        genbanks_found = 0
        for path in relevant_paths:
            if is_fasta(path):
                fastas_found += 1
                link_path = self.TEMP_GENOME_DIR / path
            elif is_genbank(path):
                genbanks_found += 1
                link_path = genbanks_temp_subfolder / path
            try:
                link_path.symlink_to(self.USER_GENOME_DIR / path)
            except FileExistsError:
                link_path.unlink()
                link_path.symlink_to(self.USER_GENOME_DIR / path)
                
        if not(fastas_found + genbanks_found):
            msg = "No fasta files or Genbank files were found in the provided genome folder!"
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        LOG.info(f"Detected {fastas_found} relevant FASTA files.")
        LOG.info(f"Detected {genbanks_found} relevant GenBank files.")
        
        # Convert Genbank files to fasta format
        if genbanks_found:
            convert_genbanks_to_fastas(genbanks_temp_subfolder, self.TEMP_GENOME_DIR, workers = self.cores)
        
        # Delete temporary Genbank subfolder
        shutil.rmtree(genbanks_temp_subfolder)
            
        LOG.info(f"Generated {fastas_found + genbanks_found} FASTA files in {self.TEMP_GENOME_DIR}")
            
        ## Add the assembly file column to the extended binary table
        LOG.info('Matching hit accessions with genome file paths')
        
        # First strip the filenames and the accessions from any suffix
        accessions_no_suffices = self.binary_df['Organism'].apply(remove_suffixes)
        accessions_no_suffices.name = 'filename'
        filenames_no_suffices = [{'path': str(f.resolve()),
                                  'filename': remove_suffixes(f.name)}
                                 for f in self.TEMP_GENOME_DIR.iterdir()]
        filenames_no_suffices = pd.DataFrame.from_records(filenames_no_suffices)
        
        # Then match the file paths to these shortened accessions
        accessions_by_filenames = pd.merge(accessions_no_suffices, filenames_no_suffices, on = 'filename', how = 'left')
        broken_links = accessions_by_filenames[accessions_by_filenames['path'].isna()]
        if not broken_links.empty:
            msg = f'No genome file found for accessions {", ".join(broken_links["filename"].to_list())}'
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        # Set if no broken links
        self.binary_df['assembly_file'] = accessions_by_filenames['path']
        
        return None
    
    