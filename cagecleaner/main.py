#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import tempfile
import sys
import logging
import warnings
from pathlib import Path
from importlib.metadata import version
from cblaster.classes import Session

from cagecleaner.local_genome_run import LocalGenomeRun
from cagecleaner.local_region_run import LocalRegionRun
from cagecleaner.remote_genome_run import RemoteGenomeRun
from cagecleaner.remote_region_run import RemoteRegionRun
from cagecleaner.validators import parse_and_validate_arguments

__version__ = version("cagecleaner")

warnings.filterwarnings(action = 'ignore', module = 'cblaster')


LOG = logging.getLogger(__name__)


def create_parser() -> argparse.Namespace:
    """
    This function creates a parser object that will collect the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        parser (argparse.ArgumentParser): An ArgumentParser object holding the CLI ready to collect the arguments when called
        
    Note:
        Also configures the logger.
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
                CAGEcleaner: A tool to reduce redundancy among gene mining hit sets.
   
                CAGEcleaner reduces redundancy by dereplicating the genome (region)s harbouring the hits. 
                It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster layout
                or an outlier homology score.
                """,
                add_help = False
                )

    args_general = parser.add_argument_group('General')
    args_general.add_argument('--cores', dest = 'cores', default = 1, type = int, help = "Number of cores available to use (default: 1).")    
    args_general.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: False).")
    args_general.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4], help = "Console verbosity level (default: 3 (info)).")
    args_general.add_argument('-np', '--no-progress', dest = "no_progress", default = False, action = 'store_true', help = "Hide most progress bars (default: False).")
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit.")      
    
    args_io = parser.add_argument_group('File inputs and outputs')
    args_io.add_argument('-s', '--session', dest = "session", type = Path, required = True, help = "Path to cblaster session (either obtained from a search run or from cagecleaner-generate-session).")
    args_io.add_argument('-g', '--genomes', dest = "genome_dir", type = Path, default = Path('.'), 
                         help = "[Only relevant for local searches] Path to local genome folder containing genome files. Accepted formats are FASTA and Genbank [.fasta; .fna; .fa; .gbff; .gbk; .gb]. Files can be gzipped. (default: current working directory)")
    args_io.add_argument('-o', '--output', dest = "output", type = Path, default = Path('.'), help = "Output directory (default: current working directory).")
    args_io.add_argument('-t', '--temp', dest = "temp", type = Path, default = Path(tempfile.gettempdir()), help = "Path to store temporary files (default: your OS's default temporary directory).")
    args_io.add_argument('--keep_downloads', dest = "keep_downloads", default = False, action = "store_true", help = "Keep downloaded genomes.")
    args_io.add_argument('--keep_dereplication', dest = "keep_dereplication", default = False, action = "store_true", help = "Keep dereplication intermediary output.")
    args_io.add_argument('--keep_intermediate', dest = "keep_intermediate", default = False, action = "store_true", help = "Keep all intermediate data. This overrules other keep flags.")
 
    args_id_io = parser.add_argument_group('Exclusions and bypasses', description = "For local searches, duplicate scaffold IDs can be further specified using the following format: <organism_ID>:<scaffold_ID>. Discard any file extension.")
    args_id_io.add_argument('-bys', '--bypass_scaffolds', dest = "bypass_scaffolds", default = '', help = "Scaffold IDs in the binary table to bypass dereplication (comma-separated). These will end up in the final output anyway.")
    args_id_io.add_argument('-byo', '--bypass_organisms', dest = "bypass_organisms", default = '', help = "Organisms in the binary table to bypass dereplication (comma-separated). These will end up in the final output anyway.")
    args_id_io.add_argument('-exs', '--exclude_scaffolds', dest = 'excluded_scaffolds', default = '', help = "Scaffold IDs in the binary table to be excluded from the hit set (comma-separated).")
    args_id_io.add_argument('-exo', '--exclude_organisms', dest = 'excluded_organisms', default = '', help = "Organisms in the binary table to be excluded from the hit set (comma-seperated).")
  
    args_download = parser.add_argument_group('Download options')
    args_download.add_argument('-w', '--workers', dest = 'download_workers', default = 2, type = int, help = "Number of download workers (default: 2).")
    args_download.add_argument('--download_batch', dest = 'download_batch', default = 300, type = int, help = "Batch size when downloading genomes (default: 300).")
    
    args_dereplication = parser.add_argument_group('Dereplication')
    args_dereplication.add_argument('--method', dest = 'method', default = "genomes", choices = ['genomes', 'regions'], type = str, 
                                    help = "Dereplication method: full genome-based ('genomes') or genomic neighbourhood-based ('regions') (default: genomes)")
    args_dereplication.add_argument('-i', '--identity', dest = 'identity', default = 99.0, type = float, help = "Identity dereplication cutoff (default: 99.0)")
    args_dereplication.add_argument('-c', '--coverage', dest = 'coverage', default = 80.0, type = float, help = "Coverage dereplication cutoff (default: 80.0)")

    args_genome_dereplication = parser.add_argument_group('Genome-based dereplication (applies skani clustering via skDER)')
    args_genome_dereplication.add_argument('--low_mem', dest = "low_mem", default = False, action = 'store_true', help = "Use skDER's low-memory mode. Lowers memory requirements substantially at the cost of a slightly lower representative quality.")
    
    args_region_dereplication = parser.add_argument_group('Region-based dereplication (applies MMseqs2 clustering)')
    args_region_dereplication.add_argument('-m', '--margin', dest = 'margin', default = 0, type = int, help = "Sequence margin at both sides of the cluster in bp. Required in case of region-based dereplication. (default: 0)")
    args_region_dereplication.add_argument('--strict', dest = 'strict_regions', default = False, action = "store_true", help = "Omit genomic regions that, including margin, are at a contig edge.")

    args_recovery = parser.add_argument_group('Hit recovery')
    args_recovery.add_argument('--no_recovery_content', dest = 'no_recovery_by_content', default = False, action = "store_true", help = "Skip recovering hits by cluster layout (default: False)")
    args_recovery.add_argument('--no_recovery_score', dest = 'no_recovery_by_score', default = False, action = "store_true", help = "Skip recovering hits by outlier homology scores (default: False)")
    args_recovery.add_argument('--min_z_score', dest = 'zscore_outlier_threshold', default = 2.0, type = float, help = "z-score threshold to consider hits outliers (default: 2.0)")
    args_recovery.add_argument('--min_score_diff', dest = 'minimal_score_difference', default = 0.1, type = float, help = "minimum score difference between hits to be considered different. Discards outlier hits with a score difference below this threshold. (default: 0.1)")

    return parser


def setup_logging(verbosity: int) -> None:
    """
    Set up the root logger if it has not been set up yet.
    
    Args:
        verbosity (int): Verbosity level (choices: 0,1,2,3,4).
        
    Returns:
        None
    """
    root_logger = logging.getLogger()
    if root_logger.handlers:
        return None
    
    log_levels = {0: logging.CRITICAL,
                  1: logging.ERROR,
                  2: logging.WARNING,
                  3: logging.INFO,
                  4: logging.DEBUG
                  }
    logging.basicConfig(
        level = log_levels[verbosity],
        format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
        datefmt="%H:%M:%S",
        handlers = [logging.StreamHandler(sys.stdout)]
        )
    
    return None
    

def main():
    # First we parse the arguments:
    parser = create_parser()
    args = parser.parse_args()
    
    # Check whether the session file exist and the verbosity value is valid. We'll definitely need these
    # to find out which workflow we need to initiate, which validation checks to run and to set up the logger.
    if not args.session.exists():
        msg = f"Session not found at {args.session}"
        LOG.error(msg)
        raise FileNotFoundError(msg)
    if args.verbosity not in [0,1,2,3,4]:
        msg = f'Invalid verbosity level: {args.verbosity}. Valid levels are 0, 1, 2, 3, and 4.'
        LOG.error(msg)
        raise ValueError(msg)
        
    # Set up logger
    setup_logging(args.verbosity)
    
    # Validate arguments
    parsed_args = parse_and_validate_arguments(args)
    
    # Initiate the approriate CAGEcleaner Run object
    source = Session.from_file(args.session).params['mode']
    method = parsed_args['method']
    mode = (source, method)
    match mode:
        case ('remote', 'genomes'):
            LOG.info('Entering remote genome mode')
            my_run = RemoteGenomeRun(parsed_args)
        case ('remote', 'regions'):
            LOG.info('Entering remote region mode')
            my_run = RemoteRegionRun(parsed_args)
        case ('local', 'genomes') | ('hmm', 'genomes'):
            LOG.info('Entering local genome mode')
            my_run = LocalGenomeRun(parsed_args)
        case ('local', 'regions') | ('hmm', 'regions'):
            LOG.info('Entering local region mode')
            my_run = LocalRegionRun(parsed_args)
        case _:
            LOG.critical(f'Invalid mode detected: {mode}. Exiting.')
            raise argparse.ArgumentError(f'Invalid search mode detected: {mode}')
    
    # Run the initialised workflow:
    my_run.run()
    
     
if __name__ == "__main__":
    main()
