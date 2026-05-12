#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import argparse
import tempfile
from pathlib import Path
from cblaster.classes import Session

from cagecleaner.file_utils import is_fasta, is_genbank


LOG = logging.getLogger(__name__)


def parse_and_validate_arguments(args: argparse.Namespace, bypass_source: None | str = None) -> dict:
    """
    This function validates the parsed arguments given through the command line.
    
    Args:
        parser (argparse.NameSpace): A NameSpace object with parsed CLI arguments
        bypass_source (None | str): Override CCL's assessment of the sequence source. Defaults to None,
            which sticks to CCL's assessment. Necessary for compatibility with csuite's validation checkers.
    
    Returns:
        parsed_args (dict): A dictionary holding the parsed and validated argument values.
        
    Raises:
        ValueError: if an invalid argument value was given.
    """
    # Determine sequence source
    match bypass_source:
        case None:
            LOG.info("--- Loading session file. ---")
            source = Session.from_file(args.session).params['mode']
        case 'local' | 'remote':
            source = bypass_source
        case _:
            raise ValueError('Overrided source is not valid! Possible choices: "local", "remote".')
    
    # Run appropriate validation checks
    method = args.method
    mode = (source, method)
    match mode:
        case ('remote', 'genomes'):
            validate_remote_genome_run_args(args)
        case ('remote', 'regions'):
            validate_remote_region_run_args(args)
        case ('local', 'genomes') | ('hmm', 'genomes'):
            validate_local_genome_run_args(args)
        case ('local', 'regions') | ('hmm', 'regions'):
            validate_local_region_run_args(args)
        case _:
            LOG.critical(f'Invalid mode detected: {mode}. Exiting.')
            raise argparse.ArgumentError(f'Invalid search mode detected: {mode}')
    
    # Parse arguments
    parsed_args = vars(args)
    
    return parsed_args


def validate_run_args(args: argparse.Namespace):
    """
    Validates the arguments processed at the level of the base Run class.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
        FileExistsError: If the output directory already exists and the force flag was not on.
    """
    try:
        if not(args.genome_dir.exists()):
            raise ValueError(f"Provided genome directory does not exits: {args.genome_dir}.")
        if not (args.coverage >= 0 and args.coverage <= 100):
            raise ValueError("Coverage threshold should be a number between 0 and 100.")
        if not(args.zscore_outlier_threshold > 0):
            raise ValueError("Z-score threshold for recovery should be greater than zero.")
        if not(args.minimal_score_difference >= 0):
            raise ValueError("Minimal score difference for recovery cannot be smaller than zero.")
        if not(args.cores > 0):
            raise ValueError("Amount of CPU cores to use should be greater than zero.")
    except ValueError as err:
        LOG.error(f'{err}')
        raise err
        
    # Output directory should not exist yet, unless flagged.
    try:
        args.output.mkdir(parents = True)
    except FileExistsError as err:
        if args.force:
            LOG.warning('Output folder already exists, but it will be overwritten.')
        else:
            LOG.error('Output folder already exists! Rerun with -f to overwrite it.')
            raise err
            
    # Temporary directory should not exist yet if not default, unless flagged.
    if args.temp != Path(tempfile.gettempdir()):
        try:
            args.temp.mkdir(parents = True)
        except FileExistsError as err:
            if args.force:
                LOG.warning('Temporary folder already exists, but it will be overwritten.')
            else:
                msg = 'Temporary folder already exists! Rerun with -f to overwrite it.'
                LOG.error(msg)
                raise err
    args.temp = Path(tempfile.mkdtemp(dir = args.temp))
            
    return None


def validate_local_run_args(args: argparse.Namespace, skip_base: bool = False):
    """
    Validates the arguments processed at the level of the intermediary Local class.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        skip_base (bool): Skip validating the base class argument values. This eliminates redundant checks when
            validating the arguments of the parent classes in this multiple inheritance class system.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    if not skip_base:
        validate_run_args(args)
    
    if args.keep_downloads == True:
        raise ValueError("Can't keep downloads in local mode.")
        
    # Make sure there is no exotic stuff in the provided genome folder
    try:
        next(filter(lambda x: is_fasta(x) or is_genbank(x), args.genome_dir.iterdir()))
    except StopIteration:
        print(args.genome_dir.resolve())
        msg = "The user-supplied genome directory does not contain any fasta or genbank file!"
        LOG.critical(msg)
        raise ValueError(msg)
        
    return None


def validate_remote_run_args(args: argparse.Namespace, skip_base: bool = False):
    """
    Validates the arguments processed at the level of the intermediary Remote class.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        skip_base (bool): Skip validating the base class argument values. This eliminates redundant checks when
            validating the arguments of the parent classes in this multiple inheritance class system.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    if not skip_base:
        validate_run_args(args)
    
    return None


def validate_genome_run_args(args: argparse.Namespace, skip_base: bool = False):
    """
    Validates the arguments processed at the level of the intermediary Genome class.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        skip_base (bool): Skip validating the base class argument values. This eliminates redundant checks when
            validating the arguments of the parent classes in this multiple inheritance class system.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid identity threshold value was passed
    """
    
    if not skip_base:
        validate_run_args(args)
    
    if not(args.identity <= 100 and args.identity >= 82):
        raise ValueError("Identity threshold should be between 82 % and 100 % in case of full-genome-based dereplication (see skani documentation).")
        
    return None


def validate_region_run_args(args: argparse.Namespace, skip_base: bool = False):
    """
    Validates the arguments processed at the level of the intermediary Region class.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        skip_base (bool): Skip validating the base class argument values. This eliminates redundant checks when
            validating the arguments of the parent classes in this multiple inheritance class system.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    if not skip_base:
        validate_run_args(args)
    
    try:
        if not(args.identity <= 100 and args.identity >= 0):
            raise ValueError("Identity threshold should be between 0 and 100 in case of region-based dereplication.")
        if not(args.margin >= 0):
            raise ValueError("Region margin cannot be negative when dereplicating regions.")
    except ValueError as err:
        LOG.error(f'{err}')
        raise err
        
    return None


def validate_local_genome_run_args(args: argparse.Namespace) -> None:
    """
    Validates the arguments for a local genome-based dereplication run.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    validate_local_run_args(args)
    validate_genome_run_args(args, skip_base = True)
    
    return None


def validate_local_region_run_args(args: argparse.Namespace) -> None:
    """
    Validates the arguments for a local region-based dereplication run.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    validate_local_run_args(args)
    validate_region_run_args(args, skip_base = True)
    
    return None


def validate_remote_genome_run_args(args: argparse.Namespace) -> None:
    """
    Validates the arguments for a remote genome-based dereplication run.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    validate_remote_run_args(args)
    validate_genome_run_args(args, skip_base = True)
    
    if not(args.download_batch > 0): 
        raise ValueError("Download batch should be larger than 0.")
        
    return None


def validate_remote_region_run_args(args: argparse.Namespace) -> None:
    """
    Validates the arguments for a remote region-based dereplication run.
    
    Args:
        args (argparse.Namespace): Namespace object holding the arguments parsed through the CLI.
        
    Returns:
        None
        
    Raises:
        ValueError: If an invalid argument value was passed
    """
    
    validate_remote_run_args(args)
    validate_region_run_args(args, skip_base = True)
    
    return None

