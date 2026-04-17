#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.utils import run_command

import logging
from pathlib import Path


LOG = logging.getLogger()


class RegionRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving region-based dereplication.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         LocalRegionRun: Region-based dereplication for hits in local sequences.
         RemoteRegionRun: Region-based dereplication for hits in remote sequences.
    """
    
    def __init__(self, args):
        """
        Initialise a RegionRun instance.
        
        Runs the base class init and checks for a valid identity threshold and sequence margin.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
            
        Raises:
            ValueError: If identity threshold is not a percentage value, or if sequence margin is not a positive number.
            
        Returns:
            None
        """
        
        super().__init__(args)
        
        try:
            if not(args.identity <= 100 and args.identity >= 0):
                raise ValueError("Identity threshold should be between 0 and 100 in case of region-based dereplication.")
            if not(args.margin >= 0):
                raise ValueError("Region margin cannot be negative when dereplicating regions.")
        except ValueError as err:
            LOG.error(f'{err}')
            raise err
        
        self.DEREP_IN_DIR: Path = self.TEMP_DIR / 'regions' # Path where the genomic regions will be saved temporarily for region-based dereplication
        
        return None
    
    
    def dereplicate_regions(self):
        """
        This method takes the path to a genomic regions folder and dereplicates them using MMseqs2.
        MMseqs2 output is stored in TEMP_DIR/derep_out.
        
        Dereplicate the gathered genome files using whole-genome ANI similarity with skDER.
        
        Sets the dereplication input directory to the full genome folder, and runs the skDER dereplication command.
        skDER output is stored in TEMP_DIR/dereplication.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If the MMseqs2 command run fails, or if the input folder is empty or does not exist.
        """
        mmseqs_verbosity = str(min(self.verbosity, 3))
        self.DEREP_OUT_DIR.mkdir(parents = True, exist_ok = True)
        
        if not(self.DEREP_IN_DIR.exists()):
            msg = "The dereplication input folder does not exist!"
            LOG.critical(msg)
            raise RuntimeError(msg)
            
        try:
            next(self.DEREP_IN_DIR.iterdir())
        except StopIteration:
            msg = "The dereplication input folder is empty!"
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        cmd = ['mmseqs', 'easy-cluster',
               *[str(p) for p in self.DEREP_IN_DIR.iterdir()],
               str(self.DEREP_OUT_DIR / 'derep'),
               str(self.DEREP_OUT_DIR / 'tmp'),
               '--min-seq-id', str(self.identity/100),
               '-c', str(self.coverage/100),
               '--threads', str(self.cores),
               '-v', mmseqs_verbosity
               ]
        
        try:
            run_command(cmd)
        except RuntimeError:
            msg = 'Dereplicating regions with MMseqs2 failed!'
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        return None
    
    