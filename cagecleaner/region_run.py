#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner import util

import logging
import shutil
import subprocess
import threading
import sys
from pathlib import Path


LOG = logging.getLogger()


class RegionRun(Run):
    
    def __init__(self, args):
        
        super().__init__(args)
        
        assert args.identity <= 100 and args.identity >= 0, "Identity threshold should be between 0 % and 100 % in case of region-based dereplication."
        assert args.margin >= 0, "Region margin cannot be negative when dereplicating regions."
        
        self.DEREP_IN_DIR: Path = self.TEMP_DIR / 'regions' # Path where the genomic regions will be saved temporarily for region-based dereplication
        
        # Region directory should not exist yet.
        try:
            self.DEREP_IN_DIR.mkdir(parents = True)
        except FileExistsError:
            if args.force:
                LOG.warning('Region folder already exists, but it will be overwritten.')
            else:
                LOG.error('Region folder already exists! Rerun with -f to overwrite it.')
                sys.exit()
                
        return None
    
    def dereplicateRegions(self):
        """
        This method takes the path to a genomic regions folder and dereplicates them using MMseqs2.
        MMseqs2 output is stored in TEMP_DIR/derep_out.
        """
        mmseqs_executable = shutil.which('mmseqs')
        mmseqs_verbosity = str(min(self.verbosity, 3))
        self.DEREP_OUT_DIR.mkdir(parents = True, exist_ok = True)
        
        cmd = [mmseqs_executable, 'easy-cluster',
               *[str(p) for p in self.DEREP_IN_DIR.iterdir()],
               str(self.DEREP_OUT_DIR / 'derep'),
               str(self.DEREP_OUT_DIR / 'tmp'),
               '--min-seq-id', str(self.identity/100),
               '-c', str(self.coverage/100),
               '--threads', str(self.cores),
               '-v', mmseqs_verbosity
               ]
        
        LOG.debug(f'Running command: {" ".join(cmd)}')
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Capture stdout and stderr in realtime and wrap it in the logs
        def mmseqs_stdout_log(s): return LOG.debug(s.rstrip())
        def mmseqs_stderr_log(s): return LOG.warning(s.rstrip())
        
        t_out = threading.Thread(target=util._stream_reader, args=(proc.stdout, mmseqs_stdout_log))
        t_err = threading.Thread(target=util._stream_reader, args=(proc.stderr, mmseqs_stderr_log))
        t_out.daemon = True
        t_err.daemon = True
        t_out.start()
        t_err.start()
    
        returncode = proc.wait()
        t_out.join()
        t_err.join()
    
        # Wrap up
        if returncode != 0:
            LOG.critical(f"MMseqs2 exited with code {returncode}")
            sys.exit()
        else:
            LOG.info('MMseqs2 finished successfully.')
        
        return None
    
    