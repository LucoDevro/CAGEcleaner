import re
import logging
import subprocess
import gzip
import os
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from Bio import SeqIO


LOG = logging.getLogger()


def remove_suffixes(string: str) -> str:
    """
    Split off any valid sequence file suffix either in the form of .<suffix> or .<suffix>.gz.
    
    Args:
        string (str): string from which the suffix should be split off, if any.
        
    Returns:
        str: the string without the sequence file suffix
    """
    
    pattern = r'\.(fasta|fna|fa|gb|gbk|gbff)(\.gz)?$' # valid extensions for fasta and genbank files, potentially gzipped
    
    return re.sub(pattern, '', string)

def is_fasta(file: str) -> bool:
    """
    Check whether the file ends in any of the accepted fasta suffices. Gzipped allowed.
    
    Args:
        file (str): filename to check
        
    Returns:
        bool: Boolean result of the check
    """
    # ? = might occur but does not have to
    # $ = this should be the end of the string
    pattern = r'\.(fasta|fna|fa)(\.gz)?$'
    if re.search(pattern, file) is None:
        return False
    else:
        return True

def is_genbank(file: str) -> bool:
    """
    Check whether the file ends in any of the accepted genbank suffices. Gzipped allowed.
    
    Args:
        file (str): filename to check
        
    Returns:
        bool: Boolean result of the check
    """
    pattern = r'\.(gbff|gbk|gb)(\.gz)?$'
    if re.search(pattern, file)==None:
        return False
    else:
        return True
    
    
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
    
    """
    contig_end = False
    
    assembly_file = Path(row['assembly_file'])
    scaffold = row['Scaffold']
    begin_cluster = row['Start']
    begin = begin_cluster - margin
    end_cluster = row['End']
    end = end_cluster + margin
    
    in_file = in_dir / assembly_file
    with gzip.open(in_file, "rt") as handle:
        # Parse sequences
        seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        
        # Select the contig to extract the region from
        scaffold_to_extract_from = seqs[scaffold]
        
        # Check for contig edges
        length = len(scaffold_to_extract_from)
        if end >= length or begin < 0:
            contig_end = True
            # If strict, skip regions that are at a contig edge
            if strict:
                return contig_end
            end = min(end, length)
            begin = max(0, begin)
            
        # Extract genomic region from assembly
        region = scaffold_to_extract_from[begin:end]
        
        # Write in a new compressed fasta file, using the original cluster coordinates as sequence ID and filename
        region.id = '|'.join([scaffold_to_extract_from.id, str(begin_cluster), str(end_cluster)])
        out_file = str(Path(out_dir / region.id)) + '.fasta.gz'
        with gzip.open(out_file, "wt") as out_handle:
            SeqIO.write(region, out_handle, "fasta")
            
    return contig_end


def _convert_one_genbank_to_fasta(input_output_paths: tuple) -> None:
    """
    Convert a genbank file to a fasta file.
    
    Converts a genbank file to a fasta file using any2fasta. Gzips the new file.
    
    Args:
        input_output_paths (tuple): Tuple of the paths of the input genbank and the output fasta file
        
    Returns:
        None
    """
    in_file, out_file = input_output_paths
    
    # Open the output file and redirect the output of any2fasta to it.
    with open(out_file, "w") as handle:
        # use -q for quiet mode, text=True because output is not in byte form.
        subprocess.run(['any2fasta', '-q', '-g', str(in_file)], stdout=handle, check=True, text=True)
    with open(out_file, 'r') as handle:
        with gzip.open(out_file.with_suffix('.fasta.gz'), "wt") as compressed_handle:
            compressed_handle.writelines(handle)
    os.remove(out_file)
    
    return None
    

def convert_genbanks_to_fastas(in_dir: Path, out_dir: Path, workers: int = 1, no_progress: bool = False) -> None:
    """
    Convert all genbank files in an input directory to fasta files.
    
    All genbank files in the input directory are converted to fasta files in parallel using any2fasta.
    
    Args:
        in_dir (Path): input directory containing the genbank files
        out_dir (Path): output directory where the new fasta files will be written
        workers (int): number of threads for parallellisation
        no_progress (bool): flag to disable showing a progress bar while converting
        
    Returns:
        None
    """
    input_output_paths = [(i, out_dir / i.with_suffix('.fasta').name) for i in in_dir.iterdir()]
    LOG.info(f'Converting {len(input_output_paths)} Genbank genomes to Fasta format.')
    thread_map(_convert_one_genbank_to_fasta, input_output_paths,
               max_workers = workers,
               leave = False,
               disable = no_progress)
    
    return None

