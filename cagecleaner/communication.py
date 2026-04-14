import re
import logging
import subprocess
import gzip
import os
import shutil
import warnings
import sys
import pandas as pd
from itertools import batched
from pathlib import Path
from Bio import Entrez
from tqdm import tqdm
from tqdm.contrib.concurrent import thread_map
from tqdm.contrib.logging import logging_redirect_tqdm
from subprocess import CalledProcessError


LOG = logging.getLogger()
warnings.filterwarnings(action = "ignore", module = "Bio")


def _download_one_region(region: tuple, out_dir: Path, max_attempts = 3) -> None:
    """
    Download a genomic region sequence from NCBI in FASTA format with automatic compression.
    
    Uses ncbi-acc-download to retrieve a specific nucleotide region from an NCBI accession.
    The downloaded FASTA file is gzipped afterwards.
    
    Args:
        region (tuple): A two-element tuple containing:
            - accession (str): NCBI nucleotide accession identifier
            - nucl_range (str): Nucleotide range in format "start:end" (e.g., "1000:2000")
        out_dir (Path): Directory where the compressed FASTA file will be saved.
        max_attempts (int): Number of retries in case of a failure. Defaults to 3.
    
    Returns:
        None
    
    Notes:
        - Output filename is constructed as "{out_dir}/{accession}|{start}|{end}.fasta.gz"
        - Retries up to 3 times on download failure before giving up.
        - The uncompressed intermediate FASTA file is automatically deleted after compression.
        - Requires ncbi-acc-download to be installed and available in the system PATH.
        - Errors are logged but do not raise exceptions; failures return None silently.
    """
    accession, nucl_range = region
    out_file = Path('|'.join([str(out_dir / accession)] + nucl_range.split(':')) + '.fasta')
    
    acc_downloader_executable = shutil.which('ncbi-acc-download')
    cmd = [acc_downloader_executable,
           '-m', 'nucleotide',
           '-F', 'fasta',
           '-o', '/dev/stdout',
           '-g', nucl_range,
           accession]
    
    for attempt in range(max_attempts):
        with open(out_file, 'w') as handle:
            try:
                subprocess.run(cmd, stdout = handle, check = True, text = True)
            except CalledProcessError:
                if attempt < max_attempts:
                    LOG.warning(f"Error downloading {accession}. Retrying...")
                else:
                    LOG.error(f"Error downloading {accession}.")
                    return None
        
    compressed_file = out_file.with_suffix('.fasta.gz')
    with open(out_file, "r") as handle:
        with gzip.open(compressed_file, "wt") as compressed_handle:
            compressed_handle.writelines(handle)
    os.remove(out_file)
    
    return None


def download_regions(regions: pd.DataFrame, directory: Path, download_workers: int, no_progress: bool = False) -> None:
    """
    Download multiple genomic sequence regions in parallel from NCBI.
    
    Accepts a DataFrame of genomic regions and uses multi-threaded downloads to retrieve
    each region from NCBI. Automatically limits concurrent downloads to a maximum of 2.
    
    Args:
        regions (pd.DataFrame): A DataFrame with at minimum the following columns:
            - 'Scaffold': NCBI nucleotide accession identifiers
            - 'Start': Start position of the region (integer or convertible to int)
            - 'End': End position of the region (integer or convertible to int)
        directory (Path): Directory where downloaded files will be saved.
        download_workers (int): Number of concurrent download threads requested. If greater
            than 2, will be automatically reduced to 2.
        no_progress (bool, optional): If True, suppress progress bar display. Defaults to False.
    
    Returns:
        None
    """
    accessions = regions['Scaffold'].to_list()
    ranges = (regions['Start'].astype(str) + ':' + regions['End'].astype(str)).to_list()
    regions = list(zip(accessions, ranges))
    
    if download_workers > 2:
        max_workers = 2
        LOG.warning("You are requesting too many download workers by NCBI policy. Limiting the number to 2 for compliance.")
    else:
        max_workers = download_workers
    
    thread_map(lambda x: _download_one_region(x, directory), regions, 
               max_workers = max_workers,
               leave = False,
               disable = no_progress)
    
    return None


def get_assembly_accessions(scaffolds: list, source: str, no_progress: bool = False, 
                            max_attempts: int = 3, batch_size: int = 100) -> list:
    """
    Retrieve Assembly accession IDs from NCBI given a list of Nucleotide scaffold IDs.
    
    Converts Nucleotide accession IDs to their associated Assembly IDs using BioPython's
    Entrez elink and esummary services. Automatically redirects WGS (Whole Genome Shotgun)
    records to their master records when requesting Genbank IDs.
    
    Args:
        scaffolds (list): A list of NCBI Nucleotide accession identifiers (e.g., GenBank or
            RefSeq accessions). WGS Genbank records (pattern: XXXX12345678.1) are automatically
            converted to their master records.
        source (str): The database to extract IDs from. Only possible values include
            'Genbank' or 'RefSeq'. This determines which synonym field is extracted from
            the Assembly summary.
        no_progress (bool, optional): If True, suppress progress bar display. Defaults to False.
        max_attempts (int): Number of retries in case of a failure. Defaults to 3.
        batch_size (int): Number of accession IDs to link in one batch. Defaults to 100.
    
    Returns:
        list: A list of Assembly accession IDs associated with the input scaffolds. Returns
            an empty list if the retrieval fails after all retry attempts.
    
    Notes:
        - Processes scaffolds in batches of 100.
        - Retries up to 3 times on network/server errors.
    
    Raises:
        AssertationError: If the specified source is not a valid NCBI Nucleotide database.
    """
    assert source in {'Genban', 'Refseq'}, 'Invalid NCBI Nucleotide database!'
    
    # Only in case of requesting Genbank IDs, redirect WGS records to their master record
    if source == 'Genbank':
        r = re.compile('([A-Z]{4,}[0-9]{8,})\.[0-9]+')
        non_wgs_scaffolds = [sc for sc in scaffolds if not r.findall(sc)]
        wgs_scaffolds = [sc for sc in scaffolds if r.findall(sc)]
        wgs_scaffolds = [wgs_sc[:-11] + '0'*9 for wgs_sc in wgs_scaffolds]
        scaffolds = wgs_scaffolds + non_wgs_scaffolds
        LOG.info(f'Redirected {len(wgs_scaffolds)} Nucleotide WGS IDs to their master records')
    
    # Get Assembly IDs in batches of 100 in max. 3 attempts
    n_batches = int((len(scaffolds) + 0.1) // batch_size) + 1 # +0.1 to make the hundreds get right
    LOG.info(f'Getting Assembly accession IDs in {n_batches} batches of {batch_size}')
    
    # Get associated Assembly objects using Elink
    with logging_redirect_tqdm(loggers = [LOG]):
        LOG.debug('Elinking')
        records = []
        for b_idx, batch in tqdm(list(enumerate(batched(scaffolds, batch_size))), leave = False, disable = no_progress):
            for attempt in range(max_attempts):
                try:
                    with Entrez.elink(dbfrom = "nucleotide", db = 'assembly', id = batch) as handle:
                        batch_records = Entrez.read(handle)
                    records += batch_records
                    break
                except:
                    if attempt+1 < max_attempts:
                        LOG.warning(f'Error Elinking batch {b_idx+1} in attempt {attempt+1}. Max. attempts: {max_attempts}')
                    else:
                        LOG.error(f'Error Elinking batch {b_idx+1} after {max_attempts} attempts. Skipping...')
        
    # Extract UIDs from Elink records
    uids = [str(l['Id']) for record in records for linksetdb in record['LinkSetDb'] for l in linksetdb['Link']]
    
    # Get Esummary objects for every UID in one go in max. 3 attempts
    LOG.debug('Retrieving Assembly IDs from Elink responses')
    for attempt in range(max_attempts):
        try:
            with Entrez.esummary(db = 'assembly', id = uids) as handle:
                summary_set = Entrez.read(handle)
            break
        except:
            if attempt+1 < max_attempts:
                LOG.warning(f'Error getting assembly IDs in attempt {attempt+1}. Max. attempts: {max_attempts}')
            else:
                LOG.error(f'Error getting assembly IDs after {max_attempts} attempts. Skipping...')
                return []
    
    # Extract Assembly accession ID from Esummary object
    summaries = summary_set['DocumentSummarySet']['DocumentSummary']
    accessions = [str(s['Synonym'][source]) for s in summaries]
    
    LOG.info(f'Got {len(accessions)} accession IDs from {source}')
    
    return accessions


def fetch_contig_lengths(contig_ids: list, max_attempts = 3):
    """
    Retrieve sequence lengths for given NCBI contig accessions using BioPython's Entrez API.
    
    Fetches summary information for a list of nucleotide contig IDs from NCBI using Efetch and extracts
    their sequence lengths. Returns results as a deduplicated DataFrame.
    
    Args:
        contig_ids (list): A list of NCBI nucleotide accession identifiers (contig IDs)
            for which to retrieve sequence lengths.
        max_attempts (int, optional): Maximum number of retry attempts on network/server errors.
            Defaults to 3.
    
    Returns:
        pd.DataFrame: A DataFrame with columns:
            - 'Scaffold': NCBI accession version identifiers
            - 'Contig_length': Sequence length in base pairs (integer)
            Duplicates are removed and the index is reset.
    """
    LOG.debug('Retrieving contig lengths via Entrez')
    for attempt in range(max_attempts):
        try:
            with Entrez.efetch(db = 'nucleotide', id = contig_ids, rettype = 'docsum') as handle:
                records = Entrez.read(handle)
        except:
            if attempt+1 < max_attempts:
                LOG.warning(f'Error getting contig lengths in attempt {attempt+1}. Max. attempts: {max_attempts}')
            else:
                LOG.error(f'Error getting contig lengths after {max_attempts} attempts. Exiting...')
                sys.exit()
                
    # Parse the response
    lengths = [int(s['Length']) for s in records]
    accessions_with_lengths = [str(s['AccessionVersion']) for s in records]
    lengths_df = pd.DataFrame({'Scaffold': accessions_with_lengths, 'Contig_length': lengths})
    lengths_df.drop_duplicates(inplace = True, ignore_index = True)
    
    return lengths_df

