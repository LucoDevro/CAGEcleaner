#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.utils import run_command

import logging


LOG = logging.getLogger()


class GenomeRun(Run):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        assert args.identity <= 100 and args.identity >= 82, "Identity threshold should be between 82 % and 100 % in case of full-genome-based dereplication (see skani documentation)."
        
        return None
    
    def dereplicateGenomes(self):
        """
        This method takes the path to a genome folder and dereplicates the genomes using skDER.
        skDER output is stored in TEMP_DIR/dereplication.
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
    
    