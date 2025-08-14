import re
import subprocess
from pathlib import Path

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
    

def convertGenbankToFasta(genome_dir: Path, out_dir: Path) -> None:
    """
    This function takes the path to a genome folder containing genbank files.
    It then uses any2fasta in a subprocess to convert them to FASTA format and store them in the out folder.
    """
    
    assert genome_dir.exists() and genome_dir.is_dir(), "Provided genome folder for GenBank conversion does not exist or is not a directory."
    assert out_dir.exists() and out_dir.is_dir(), "Provided output folder for GenBank conversion does not exist or is not a directory."
    
    # Loop over the files in the directory:
    for file in genome_dir.iterdir():
        # Define the output file
        out_file = out_dir / removeSuffixes(file.name)
        out_file = out_file.with_suffix('.fasta')
        # Open the output file and redirect the output of any2fasta to it.
        with out_file.open('w') as out_file:
            # use -q for quiet mode, text=True because output is not in byte form.
            subprocess.run(['any2fasta', '-q', str(file)], stdout=out_file, check=True, text=True)
    
    return None