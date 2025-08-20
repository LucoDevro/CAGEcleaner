#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Import util from the directory above:

import argparse
import tempfile
from pathlib import Path
from importlib.metadata import version
from .classes import Run

__version__ = version("cagecleaner")

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
                (c) 2025 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                CAGEcleaner: A tool to remove redundancy from cblaster hits.
   
                CAGEcleaner reduces redundancy in cblaster hit sets by dereplicating the genomes containing the hits. 
                It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster content
                or an outlier cblaster score.
                """,
                add_help = False
                )

    args_general = parser.add_argument_group('General')
    args_general.add_argument('-c', '--cores', dest = 'cores', default = 1, type = int, help = "Number of cores to use (default: 1)")    
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      
    args_general.add_argument('--verbose', dest = 'verbose', default = False, action = 'store_true', help = "Enable verbose logging")
    
    args_io = parser.add_argument_group('File inputs and outputs')
    args_io.add_argument('-s', '--session', dest = "session_file", type = Path, help = "Path to cblaster session file", required = True)
    args_io.add_argument('-g', '--genomes', dest = "genome_dir", type = Path, default = '.', help = "[Only relevant for local cblaster sessions] Path to local genome folder containing genome files. Accepted formats are FASTA and GenBank [.fasta; .fna; .fa; .gbff; .gbk; .gb]. Files can be gzipped. Folder can contain other files. (default: current working directory)")
    args_io.add_argument('-o', '--output', dest = "output_dir", type = Path, default = '.', help = "Output directory (default: current working directory)")
    args_io.add_argument('-t', '--temp', dest = "temp_dir", type = Path, default = tempfile.gettempdir(), help = "Path to store temporary files (default: your system's default temporary directory).")
    args_io.add_argument('--keep_downloads', dest = "keep_downloads", default = False, action = "store_true", help = "Keep downloaded genomes")
    args_io.add_argument('--keep_dereplication', dest = "keep_dereplication", default = False, action = "store_true", help = "Keep skDER output")
    args_io.add_argument('--keep_intermediate', dest = "keep_intermediate", default = False, action = "store_true", help = "Keep all intermediate data. This overrules other keep flags.")
 
    #! Arguments for by-pass scaffolds or assemblies:
    args_id_io = parser.add_argument_group('Analysis inputs and outputs', description = "For local cblaster sessions, duplicate scaffold IDs can be further specified using the following format: <organism_ID>:<scaffold_ID>. Discard any file extension.")
    args_id_io.add_argument('-bys', '--bypass_scaffolds', dest = "bypass_scaffolds", default = '', help = "Scaffold IDs in the binary table that should bypass dereplication (comma-separated). These will end up in the final output in any case.")
    args_id_io.add_argument('-byo', '--bypass_organisms', dest = "bypass_organisms", default = '', help = "Organisms in the binary table that should bypass dereplication (comma-separated). These will end up in the final output in any case.")
    args_id_io.add_argument('-exs', '--exclude_scaffolds', dest = 'excluded_scaffolds', default = '', help = "Scaffolds IDs in the binary table to be excluded from the hit set (comma-separated). ")
    args_id_io.add_argument('-exo', '--exclude_organisms', dest = 'excluded_organisms', default = '', help = "Organisms in the binary table to be excluded from the hit set (comma-seperated).")
  
    args_download = parser.add_argument_group('Download')
    args_download.add_argument('--download_batch', dest = 'download_batch', default = 300, type = int, help = "Number of genomes to download in one batch (default: 300)")
    
    args_dereplication = parser.add_argument_group('Dereplication')
    args_dereplication.add_argument('-a', '--ani', dest = 'ani', default = 99.0, type = float, help = "ANI dereplication threshold (default: 99.0)")
    
    args_recovery = parser.add_argument_group('Hit recovery')
    args_recovery.add_argument('--no_recovery_content', dest = 'no_recovery_by_content', default = False, action = "store_true", help = "Skip recovering hits by cluster content (default: False)")
    args_recovery.add_argument('--no_recovery_score', dest = 'no_recovery_by_score', default = False, action = "store_true", help = "Skip recovering hits by outlier scores (default: False)")
    args_recovery.add_argument('--min_z_score', dest = 'zscore_outlier_threshold', default = 2.0, type = float, help = "z-score threshold to consider hits outliers (default: 2.0)")
    args_recovery.add_argument('--min_score_diff', dest = 'minimal_score_difference', default = 0.1, type = float, help = "minimum cblaster score difference between hits to be considered different. Discards outlier hits with a score difference below this threshold. (default: 0.1)")

    args = parser.parse_args()
        
    return args

def main():
    
    # First we parse the arguments:
    args = parseArguments()
    
    # Initiate a CAGECLEANER Run object:
    my_run = Run.fromArgs(args)  # This is now a LocalRun or RemoteRun depending on the mode of the session in args
    
    # Run the entire Local or Remote workflow:
    my_run.run()
    
     
if __name__ == "__main__":
    main()
