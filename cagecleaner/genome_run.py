#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.communication import _stream_reader

import logging
import shutil
import subprocess
import threading
import sys


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
        skder_executable = shutil.which('skder')
        
        cmd = [skder_executable,
               '-g', str(self.DEREP_IN_DIR),
               '-o', str(self.DEREP_OUT_DIR),
               '-i', str(self.identity),
               '-f', str(self.coverage),
               '-c', str(self.cores),
               '-d', "low_mem_"*self.low_mem + 'greedy',
               '-n'
               ]
        
        LOG.debug(f'Running command: {" ".join(cmd)}')
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Capture stdout and stderr in realtime and wrap it in the logs
        def skder_stdout_log(s): return LOG.debug(s.rstrip())
        def skder_stderr_log(s): return LOG.warning(s.rstrip())
        
        t_out = threading.Thread(target=_stream_reader, args=(proc.stdout, skder_stdout_log))
        t_err = threading.Thread(target=_stream_reader, args=(proc.stderr, skder_stderr_log))
        t_out.daemon = True
        t_err.daemon = True
        t_out.start()
        t_err.start()
    
        returncode = proc.wait()
        t_out.join()
        t_err.join()
    
        # Wrap up
        if returncode != 0:
            LOG.critical(f"skDER exited with code {returncode}")
            sys.exit()
        else:
            LOG.info('skDER finished successfully.')
        
        LOG.info("Dereplication done!")
        
        extensions = {'.fna','.fa','.fasta','.fna.gz','.fa.gz','.fasta.gz'}
        paths = [str(p) for p in self.DEREP_IN_DIR.iterdir() if extensions & set(p.suffixes)]
        before = len(paths)
        after = len(list((self.DEREP_OUT_DIR / 'Dereplicated_Representative_Genomes').iterdir()))
        LOG.info(f'{before} genomes were reduced to {after} genomes.')
        
        return None
    
    