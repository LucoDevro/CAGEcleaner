#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run
from cagecleaner.utils import run_command
from cagecleaner.file_utils import parse_cdhit_clustering

import logging
import pandas as pd
from pathlib import Path


LOG = logging.getLogger(__name__)


class RegionRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving region-based dereplication.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         LocalRegionRun: Region-based dereplication for hits in local sequences.
         RemoteRegionRun: Region-based dereplication for hits in remote sequences.
    """
    
    def __init__(self, parsed_args):
        """
        Initialise a RegionRun instance.
        
        Runs the base class init and checks for a valid identity threshold and sequence margin.
        
        Args:
            parsed_args (dict): Parsed and validated command-line arguments
            
        Raises:
            ValueError: If identity threshold is not a percentage value, or if sequence margin is not a positive number.
            
        Returns:
            None
        """
        
        super().__init__(parsed_args)
        
        self.DEREP_IN_DIR: Path = self.TEMP_DIR / 'regions' # Path where the genomic regions will be saved temporarily for region-based dereplication
        self.DEREP_IN_DIR.mkdir(parents = True)
        
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
               '-v', mmseqs_verbosity,
               '--kmer-per-seq', str(80),
               '--max-seqs', str(300),
               '--cluster-reassign',
               '--seq-id-mode', str(1),
               ]
        
        try:
            run_command(cmd)
        except RuntimeError:
            msg = 'Dereplicating regions with MMseqs2 failed!'
            LOG.critical(msg)
            raise RuntimeError(msg)
        
        return None
    
    
    def reassign_clusters(self):
        """
        Reassign cluster representatives.
        
        Reassigns cluster representatives using an a posteriori CD-HIT wide-bandwidth clustering.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If the CD-HIT command run fails.
        """
        # Rename uncorrected files
        for uncorr_file_old in self.DEREP_OUT_DIR.glob('derep_*'):
            uncorr_file = uncorr_file_old.with_name(uncorr_file_old.name.replace('derep_', 'derep_uncorr_'))
            uncorr_file_old.rename(uncorr_file)
        
        # Launch a CD-HIT clustering for the representatives
        cmd = ['cd-hit-est',
               '-i', str(self.DEREP_OUT_DIR / 'derep_uncorr_rep_seq.fasta'),
               '-o', str(self.DEREP_OUT_DIR / 'derep_rep_seq.fasta'),
               '-c', str(self.identity/100),
               '-A', str(self.coverage/100),
               '-b', str(self.cdhit_bandwidth),
               '-T', str(self.cores),
               '-M', str(0),
               '-d', str(0),
               ]
        try:
            run_command(cmd)
        except RuntimeError:
            msg = 'Reassigning clusters with CD-HIT failed!'
            LOG.critical(msg)
            raise RuntimeError(msg)
            
        # Parse CD-HIT clustering file
        cdhit_clustering = parse_cdhit_clustering(self.DEREP_OUT_DIR / 'derep_rep_seq.fasta.clstr')
            
        # Reassign representatives directly to a corrected MMseqs clustering file
        mmseqs_uncorr = pd.read_table(self.DEREP_OUT_DIR / 'derep_uncorr_cluster.tsv',
                                      header = None, names = ['representative', 'region'])
        mmseqs_uncorr['representative'] = mmseqs_uncorr['representative'].replace(cdhit_clustering)
        mmseqs_uncorr.to_csv(self.DEREP_OUT_DIR / 'derep_cluster.tsv', sep = "\t", header = False, index = False)
            
        return None
    
    