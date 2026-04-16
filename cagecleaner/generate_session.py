#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.utils import generate_cblaster_session

import argparse
import sys
import logging
from pathlib import Path


LOG = logging.getLogger()
logging.basicConfig(
    level = 3,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)]
    )


def parse_arguments() -> argparse.Namespace:
    """
    This function parses the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        A Namespace object holding the parsed arguments
        
    Note:
        Also configures the logger.
    """
    
    # Initiate an ArgumentParser object
    parser = argparse.ArgumentParser(
        prog = 'cagecleaner-generate-session',
                epilog = 
                """
                Lucas De Vrieze, Miguel Biltjes
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                Helper tool to generate a cblaster session from TSV tables.
                """,
                add_help = False
                )

    parser.add_argument('--clusters', dest = "clusters", type = Path, required = True, help = "Path to clusters TSV file.")
    parser.add_argument('--hits', dest = 'hits', type = Path, required = True, help = 'Path to hits TSV file.')
    parser.add_argument('--queries', dest = "queries", type = Path, required = True, help = "Path to queries TSV file.")
    parser.add_argument('-m', '--mode', dest = "mode", type = str, required = True, choices = ['local', 'remote'],
                        help = "Dereplication mode ('remote' for dereplication using NCBI IDs and server queries; 'local' for dereplication using local sequence DBs).")
    parser.add_argument('-s', '--session', dest = "session", type = Path, default = Path('session.json'), help = "Output path (default: session.json)")
    parser.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting session (default: False).")
    parser.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit.")
    
    # Parse arguments
    args = parser.parse_args()
    
    return args

def main():
    
    # Parse arguments
    args = parse_arguments()
    
    # Check output path
    if args.session.exists():
        if args.force:
            LOG.warning('Output path already exists, but it will be overwritten.')
        else:
            LOG.error('Output path already exists! Rerun with -f to overwrite it.')
            sys.exit()
    
    # Generate session
    LOG.info("Generating cblaster session.")
    session = generate_cblaster_session(args.hits, args.clusters, args.queries, args.mode)
    
    # Save session on harddisk
    LOG.info(f"Writing session to disk at {args.session.resolve()}")
    with open(args.session, "w") as handle:
        session.to_json(handle)
    
    LOG.info("DONE!")
        

if __name__ == "__main__":
    main()
    
    