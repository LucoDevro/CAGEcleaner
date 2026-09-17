#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import logging
import gzip
import gb_io
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from tqdm.contrib.logging import logging_redirect_tqdm
from Bio import SeqIO


LOG = logging.getLogger(__name__)


SEQ_FILE_EXTENSIONS = r'\.(fasta|fna|fa|gb|gbk|gbff)(\.gz)?$' # valid extensions for fasta and genbank files, potentially gzipped
FASTA_FILE_EXTENSIONS = r'\.(fasta|fna|fa)(\.gz)?$'
GENBANK_FILE_EXTENSIONS = r'\.(gbff|gbk|gb)(\.gz)?$'


def remove_suffixes(string: str | Path) -> str:
    """
    Split off any valid sequence file suffix either in the form of .<suffix> or .<suffix>.gz.
    
    Args:
        string (str | Path): string from which the suffix should be split off, if any.
        
    Returns:
        str: the string without the sequence file suffix
    """
    
    return re.sub(SEQ_FILE_EXTENSIONS, '', str(string))


def is_fasta(file: str | Path) -> bool:
    """
    Check whether the file ends in any of the accepted fasta suffices. Gzipped allowed.
    
    Args:
        file (str | Path): filename to check
        
    Returns:
        bool: Boolean result of the check
    """
    if re.search(FASTA_FILE_EXTENSIONS, str(file)) is None:
        return False
    else:
        return True
    

def is_genbank(file: str | Path) -> bool:
    """
    Check whether the file ends in any of the accepted genbank suffices. Gzipped allowed.
    
    Args:
        file (str | Path): filename to check
        
    Returns:
        bool: Boolean result of the check
    """
    if re.search(GENBANK_FILE_EXTENSIONS, str(file)) is None:
        return False
    else:
        return True
    
    
def read_file(file: str | Path, mode: str = "t"):
    """
    Open the appropriate file handle for a genome file.
    
    Automatically distinguishes between compressed and uncompressed files based on the file extension.
    
    Args:
        file (str | Path): genome file to open
        mode (str): file mode (text ('t'; default) or binary ('b'))
        
    Returns:
        handle: A file handle to open the genome file
    """
    if '.gz' in Path(file).suffixes:
        open_file = gzip.open(file, mode = f'r{mode}')
    else:
        open_file = open(file, mode = 'r')
        
    return open_file
    
    
def _extract_one_region(row: dict, margin: int, in_dir: Path, out_dir: Path, strict: bool) -> bool:
    """
    Extract a genomic region from a fasta file.
    
    A genomic region is extracted from a fasta file parsed using BioPython's SeqIO. The original cluster
    genomic coordinates are used in the filename and the sequence ID of the extracted region to make them
    appear in MMseqs2 dereplication table downstream in the pipeline. Extracted region's sequence files are gzipped.
    
    Regions that extend beyond contig boundaries are treated as specified by the user (strict_regions flag).
    When strict_regions is enabled, regions at contig edges are skipped.When disabled (permissive mode),
    such regions are retained but clipped to the contig boundaries.
    
    Args:
        row (dict): dictionary form of a row from the binary table
        margin (int): size of the sequence margin to add to both sides of the cluster before extraction
        in_dir (Path): path of the input directory containing the genome assembly files
        out_dir (Path): path of the output directory to which the region files will be extracted.
        strict (bool): flag for contig edge behaviour (strict discarding or permissive clipping)
    
    Returns:
        bool: Whether the extracted region was at a contig edge
        
    Raises:
        FileNotFoundError: If the input file cannot be found.
        FileNotFoundError: If the output file cannot be opened.
    
    """
    contig_end = False
    
    assembly_file = Path(row['assembly_file'])
    scaffold = row['Scaffold']
    begin_cluster = row['Start']
    begin = begin_cluster - margin
    end_cluster = row['End']
    end = end_cluster + margin
    
    in_file = in_dir / assembly_file
    try:
        with read_file(in_file) as handle:
            # Parse sequences
            seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    except FileNotFoundError as err:
        LOG.error(err)
        raise err
            
    # Select the contig to extract the region from
    scaffold_to_extract_from = seqs[scaffold]
    
    # Check for contig edges
    length = len(scaffold_to_extract_from)
    if end >= length or begin <= 0:
        contig_end = True
        # If strict, skip regions that are at a contig edge
        if strict:
            return contig_end
        end = min(end, length)
        begin = max(0, begin)
        
    # Extract genomic region from assembly
    region = scaffold_to_extract_from[begin:end]
    
    # Write in a new compressed fasta file, using the original cluster coordinates as sequence ID and filename
    region.id = '§'.join([scaffold_to_extract_from.id, str(begin_cluster), str(end_cluster)])
    # Flexible search pattern to remove sequence file extensions
    out_file = Path(re.sub(SEQ_FILE_EXTENSIONS, '', str(out_dir / region.id)) + '.fasta.gz')
    try:
        with gzip.open(out_file, "wt") as out_handle:
            SeqIO.write(region, out_handle, "fasta")
    except FileNotFoundError as err:
        LOG.error(err)
        raise err
    
    return contig_end


def _convert_one_genbank_to_fasta(input_output_paths: tuple[Path]) -> None:
    """
    Convert a genbank file to a fasta file.
    
    Converts a genbank file to a fasta file using a gb-io iterator. Genbank file may be compressed,
    which is handled by read_file().
    
    Args:
        input_output_paths (tuple[Path]): Tuple of the paths of the input genbank and the output fasta file
        
    Returns:
        None
        
    Raises:
        FileNotFoundError: If any of the paths cannot be opened.
    """
    in_file, out_file = input_output_paths
    
    try:
        with open(out_file.with_suffix('.fasta'), "wb") as fasta_handle:
            with read_file(in_file, mode = "b") as genbank_handle:
                for record in gb_io.iter(genbank_handle):
                    fasta_handle.write((f">{record.version}\n".encode()))
                    fasta_handle.write(record.sequence + b'\n')
            
    except FileNotFoundError as err:
        LOG.error(err)
        raise err
    
    return None
    

def convert_genbanks_to_fastas(in_dir: Path, out_dir: Path, workers: int = 1, no_progress: bool = False) -> None:
    """
    Convert all genbank files in an input directory to compressed fasta files.
    
    All genbank files in the input directory are converted to compressed fasta
    files in parallel using any2fasta and gzip.
    
    Args:
        in_dir (Path): input directory containing the genbank files
        out_dir (Path): output directory where the new fasta files will be written
        workers (int): number of threads for parallellisation
        no_progress (bool): flag to disable showing a progress bar while converting
        
    Returns:
        None
    """
    input_output_paths = [(i, Path(re.sub(SEQ_FILE_EXTENSIONS, '', str(out_dir / i.name)) + '.fasta'))
                          for i in in_dir.iterdir()]
    LOG.info(f'Converting {len(input_output_paths)} Genbank genomes to Fasta format.')
    with logging_redirect_tqdm(loggers = [LOG]):
        thread_map(_convert_one_genbank_to_fasta, input_output_paths,
                   max_workers = workers,
                   leave = False,
                   disable = no_progress)
    
    return None


def parse_cdhit_clustering(path: Path) -> dict:
    """
    Parse a CD-HIT clustering file.
    
    Parses a CD-HIT clustering file into a dictionary with key the ID of the representative region
    and value the ID of a cluster member region.
    
    Args:
        path (Path): path of the CD-HIT clustering file
        
    Returns:
        clustering (dict): parsed clustering
    """
    # Initialise containers
    clustering = {}
    new_clustering_lines = {}
    representative = ""
    region_ids = []
    
    # Iterate over all lines in the file
    with open(path, "r") as handle:
        for line in handle:
            # If a new cluster starts, generate and save the new clustering lines, and reset
            if line.startswith('>'):
                # Save new cluster
                new_clustering_lines = {rid: representative for rid in region_ids}
                clustering |= new_clustering_lines
                
                # Reset containers and continue iterating
                region_ids = []
                representative = ""
                new_clustering_lines = {}
                continue
            
            # Collect region IDs
            region_id = re.findall('>.+\.\.\.', line)[0]
            region_id = region_id[1:-3]
            region_ids.append(region_id)
            
            # Record the representative when found
            if line.endswith("*\n"):
                representative = region_id
        
        # Save final cluster
        new_clustering_lines = {rid: representative for rid in region_ids}
        clustering |= new_clustering_lines
        
    return clustering

