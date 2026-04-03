#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.local_genome_run import LocalGenomeRun
from cagecleaner.local_region_run import LocalRegionRun
from cagecleaner.remote_genome_run import RemoteGenomeRun
from cagecleaner.remote_region_run import RemoteRegionRun

import argparse
import tempfile
import sys
import logging
from pathlib import Path
from importlib.metadata import version
from cblaster.classes import Session

__version__ = version("cagecleaner")

LOG = logging.getLogger()

def parseArguments():
    """
    This function parses the arguments given through the command line.
    """
    
    # Initiate an ArgumentParser object
    parser = argparse.ArgumentParser(
        prog = 'cagecleaner',
                epilog = 
                """
                Lucas De Vrieze, Miguel Biltjes
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                CAGEcleaner: A tool to remove redundancy among gene mining hits.
   
                CAGEcleaner reduces redundancy in cblaster hit sets by dereplicating the genome (region)s containing the hits. 
                It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster content
                or an outlier homology score.
                """,
                add_help = False
                )

    args_general = parser.add_argument_group('General')
    args_general.add_argument('--cores', dest = 'cores', default = 1, type = int, help = "Number of cores available to use (default: 1)")    
    args_general.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: False).")
    args_general.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4], help = "Console verbosity level (default: 3 (info))")
    args_general.add_argument('-np', '--no-progress', dest = "no_progress", default = False, action = 'store_true', help = "Don't show progress bar (default: False).")
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      
    
    args_io = parser.add_argument_group('File inputs and outputs')
    args_io.add_argument('-s', '--session', dest = "session_file", type = Path, help = "Path to cblaster session file", required = True)
    args_io.add_argument('-g', '--genomes', dest = "genome_dir", type = Path, default = '.', help = "[Only relevant for local cblaster sessions] Path to local genome folder containing genome files. Accepted formats are FASTA and GenBank [.fasta; .fna; .fa; .gbff; .gbk; .gb]. Files can be gzipped. Folder can contain other files. (default: current working directory)")
    args_io.add_argument('-o', '--output', dest = "output_dir", type = Path, default = '.', help = "Output directory (default: current working directory)")
    args_io.add_argument('-t', '--temp', dest = "temp_dir", type = Path, default = tempfile.gettempdir(), help = "Path to store temporary files (default: your OS's default temporary directory).")
    args_io.add_argument('--keep_downloads', dest = "keep_downloads", default = False, action = "store_true", help = "Keep downloaded genomes")
    args_io.add_argument('--keep_dereplication', dest = "keep_dereplication", default = False, action = "store_true", help = "Keep dereplication intermediary output")
    args_io.add_argument('--keep_intermediate', dest = "keep_intermediate", default = False, action = "store_true", help = "Keep all intermediate data. This overrules other keep flags.")
 
    #! Arguments for bypassing scaffolds or assemblies:
    args_id_io = parser.add_argument_group('Analysis inputs and outputs', description = "For local cblaster sessions, duplicate scaffold IDs can be further specified using the following format: <organism_ID>:<scaffold_ID>. Discard any file extension.")
    args_id_io.add_argument('-bys', '--bypass_scaffolds', dest = "bypass_scaffolds", default = '', help = "Scaffold IDs in the binary table that should bypass dereplication (comma-separated). These will end up in the final output in any case.")
    args_id_io.add_argument('-byo', '--bypass_organisms', dest = "bypass_organisms", default = '', help = "Organisms in the binary table that should bypass dereplication (comma-separated). These will end up in the final output in any case.")
    args_id_io.add_argument('-exs', '--exclude_scaffolds', dest = 'excluded_scaffolds', default = '', help = "Scaffolds IDs in the binary table to be excluded from the hit set (comma-separated). ")
    args_id_io.add_argument('-exo', '--exclude_organisms', dest = 'excluded_organisms', default = '', help = "Organisms in the binary table to be excluded from the hit set (comma-seperated).")
  
    args_download = parser.add_argument_group('Download options')
    args_download.add_argument('-w', '--workers', dest = 'download_workers', default = 1, type = int, help = "Number of download workers (default: 1).")
    args_download.add_argument('--download_batch', dest = 'download_batch', default = 300, type = int, help = "Number of genomes to download in one batch when downloading genomes (default: 300)")
    
    args_dereplication = parser.add_argument_group('Dereplication')
    args_dereplication.add_argument('--method', dest = 'method', default = "genomes", choices = ['genomes', 'regions'], type = str, 
                                    help = "Dereplication method: full genome-based ('genomes') or genomic neighbourhood-based ('regions') (default: genomes)")
    args_dereplication.add_argument('-i', '--identity', dest = 'identity', default = 99.0, type = float, help = "Identity dereplication cutoff (default: 99.0)")
    args_dereplication.add_argument('-c', '--coverage', dest = 'coverage', default = 80.0, type = float, help = "Coverage dereplication cutoff (default: 80.0)")

    args_genome_dereplication = parser.add_argument_group('Genome-based dereplication (applies skani clustering via skDER)')
    args_genome_dereplication.add_argument('--low_mem', dest = "low_mem", default = False, action = 'store_true', help = "Use skDER's low-memory mode. Lowers memory requirements substantially at the cost of a slightly lower representative quality.")
    
    args_region_dereplication = parser.add_argument_group('Region-based dereplication (applies MMseqs2 clustering)')
    args_region_dereplication.add_argument('-m', '--margin', dest = 'margin', default = 0, type = int, help = "Sequence margin to add to both sides of the cluster hit in bp. Required in case of region-based dereplication. (default: 0)")
    args_region_dereplication.add_argument('--strict', dest = 'strict_regions', default = False, action = "store_true", help = "Omit genomic regions that, including margin, are at a contig edge.")

    args_recovery = parser.add_argument_group('Hit recovery')
    args_recovery.add_argument('--no_recovery_content', dest = 'no_recovery_by_content', default = False, action = "store_true", help = "Skip recovering hits by cluster content (default: False)")
    args_recovery.add_argument('--no_recovery_score', dest = 'no_recovery_by_score', default = False, action = "store_true", help = "Skip recovering hits by outlier scores (default: False)")
    args_recovery.add_argument('--min_z_score', dest = 'zscore_outlier_threshold', default = 2.0, type = float, help = "z-score threshold to consider hits outliers (default: 2.0)")
    args_recovery.add_argument('--min_score_diff', dest = 'minimal_score_difference', default = 0.1, type = float, help = "minimum score difference between hits to be considered different. Discards outlier hits with a score difference below this threshold. (default: 0.1)")

    args = parser.parse_args()
    
    # Set up logger
    log_levels = {0: logging.CRITICAL,
                  1: logging.ERROR,
                  2: logging.WARNING,
                  3: logging.INFO,
                  4: logging.DEBUG
                  }
    logging.basicConfig(
        level = log_levels[args.verbosity],
        format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
        datefmt="%H:%M:%S",
        handlers = [logging.StreamHandler(sys.stdout)]
        )
        
    return args

def main():
    
    # First we parse the arguments:
    args = parseArguments()
    
    # Initiate the approriate CAGEcleaner Run object:
    LOG.info("--- Loading session file. ---")
    source = Session.from_file(args.session_file).params['mode']
    method = args.method
    mode = (source, method)
    match mode:
        case ('remote', 'genomes'):
            LOG.info('Entering remote genome mode')
            my_run = RemoteGenomeRun(args)
        case ('remote', 'regions'):
            LOG.info('Entering remote region mode')
            my_run = RemoteRegionRun(args)
        case ('local', 'genomes') | ('hmm', 'genomes'):
            LOG.info('Entering local genome mode')
            my_run = LocalGenomeRun(args)
        case ('local', 'regions') | ('hmm', 'regions'):
            LOG.info('Entering local region mode')
            my_run = LocalRegionRun(args)
        case _:
            LOG.critical(f'Invalid mode detected: {mode}. Exiting.')
            sys.exit()
    
    # Run the initialised workflow:
    my_run.run()
    
     
if __name__ == "__main__":
    main()
