#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.utils import run_command

import logging


LOG = logging.getLogger()


class GenomeRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving whole-genome dereplication.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         LocalGenomeRun: Full-genome dereplication for hits in local sequences.
         RemoteGenomeRun: Full-genome dereplication for hits in remote sequences.
    """
    
    def __init__(self, args):
        """
        Initialise a GenomeRun instance.
        
        Runs the base class init and checks for a valid ANI threshold.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
            
        Raises:
            AssertionError: If identity threshold is not between 82 and 100, the threshold region support by skani.
            
        Returns:
            None
        """
        
        super().__init__(args)
        
        assert args.identity <= 100 and args.identity >= 82, "Identity threshold should be between 82 % and 100 % in case of full-genome-based dereplication (see skani documentation)."
        
        return None
    
    def dereplicate_genomes(self):
        """
        Dereplicate the gathered genome files using whole-genome ANI similarity with skDER.
        
        Sets the dereplication input directory to the full genome folder, and runs the skDER dereplication command.
        skDER output is stored in TEMP_DIR/dereplication.
        
        Returns:
            None
        """
        self.DEREP_IN_DIR = self.TEMP_GENOME_DIR
        
        LOG.info(f'Dereplicating genomes in {str(self.DEREP_IN_DIR)} with identity cutoff of {str(self.identity)} % and coverage cutoff of {str(self.coverage)} %')
        
        LOG.info("Starting skDER")
        
        cmd = ['skder',
               '-g', str(self.DEREP_IN_DIR),
               '-o', str(self.DEREP_OUT_DIR),
               '-i', str(self.identity),
               '-f', str(self.coverage),
               '-c', str(self.cores),
               '-d', "low_mem_"*self.low_mem + 'greedy',
               '-n'
               ]
        run_command(cmd)
        
        LOG.info("Dereplication done!")
        
        extensions = {'.fna','.fa','.fasta','.fna.gz','.fa.gz','.fasta.gz'}
        paths = [str(p) for p in self.DEREP_IN_DIR.iterdir() if extensions & set(p.suffixes)]
        before = len(paths)
        after = len(list((self.DEREP_OUT_DIR / 'Dereplicated_Representative_Genomes').iterdir()))
        LOG.info(f'{before} genomes were reduced to {after} genomes.')
        
        return None
    
    