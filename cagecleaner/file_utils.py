import re
import logging
import subprocess
import gzip
import os
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from Bio import SeqIO


LOG = logging.getLogger(__name__)


def removeSuffixes(string: str) -> str:
    """
    Splits off any suffix either in the form of .<suffix> or .<suffix>.gz.
    Is used in local mode to ensure that the 'Organism' column in the binary table matches exactly with the name of the genome file.
    """
    pattern = r'\.(fasta|fna|fa|gff|gff3|gb|gbk|gbff)(\.gz)?$'
    
    # Substitute the match with an empty string and return:
    return re.sub(pattern, '', string)

def isFasta(file: str) -> bool:
    """
    Returns true if the file ends in any of the accepted fasta suffices. Gzipped allowed.
    """
    # ? = might occur but does not have to
    # $ = this should be the end of the string
    pattern = r'\.(fasta|fna|fa)(\.gz)?$'
    if re.search(pattern, file) is None:
        return False
    else:
        return True

def isGff(file: str) -> bool:
    """
    Returns true if the file ends in any of the accepted gff suffices. Gzipped allowed.
    """
    pattern = r'\.(gff|gff3)(\.gz)?$'
    if re.search(pattern, file) is None:
        return False
    else:
        return True

def isGenbank(file: str) -> bool:
    """
    Returns true if the file ends in any of the accepted genbank suffices. Gzipped allowed.
    """
    pattern = r'\.(gbff|gbk|gb)(\.gz)?$'
    if re.search(pattern, file)==None:
        return False
    else:
        return True
    
    
def _extractOneRegion(row: dict, margin: int, in_dir: Path, out_dir: Path, strict: bool) -> bool:
    assembly = row['assembly_file']
    scaffold = row['Scaffold']
    begin = row['Start'] - margin
    end = row['End'] + margin
    contig_end = False
    
    with gzip.open(in_dir / assembly, "rt") as handle:
        seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        # If strict, skip regions that are at a contig edge
        scaffold_to_extract_from = seqs[scaffold]
        length = len(scaffold_to_extract_from)
        if end >= length or begin < 0:
            contig_end = True
            if strict:
                return contig_end
            end = min(end, length)
            begin = max(0, begin)
        # Extract genomic region from assembly
        region = scaffold_to_extract_from[begin:end]
        # Write in a new compressed fasta file
        with gzip.open(out_dir / assembly, "wt") as out_handle:
            SeqIO.write(region, out_handle, "fasta")
            
    return contig_end


def _convertOneGenbankToFasta(input_output_paths: tuple) -> None:
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
    

def convertGenbankToFasta(genome_dir: Path, out_dir: Path, workers: int = 1, no_progress: bool = False) -> None:
    """
    This function takes the path to a genome folder containing genbank files.
    It then uses any2fasta in a subprocess to convert them to FASTA format and store them in the out folder.
    """
    
    assert genome_dir.exists() and genome_dir.is_dir(), "Provided genome folder for GenBank conversion does not exist or is not a directory."
    assert out_dir.exists() and out_dir.is_dir(), "Provided output folder for GenBank conversion does not exist or is not a directory."
    
    # Convert the GenBank files in parallel
    ios = [(i, out_dir / i.with_suffix('.fasta').name) for i in genome_dir.iterdir()]
    LOG.info(f'Converting {len(ios)} Genbank genomes to Fasta format.')
    thread_map(_convertOneGenbankToFasta, ios,
               max_workers = workers,
               leave = False,
               disable = no_progress)
    
    return None

