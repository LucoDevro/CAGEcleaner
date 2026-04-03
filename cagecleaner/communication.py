import re
import logging
import subprocess
import gzip
import os
import shutil
import warnings
from itertools import batched
from pathlib import Path
from Bio import Entrez
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from subprocess import CalledProcessError


LOG = logging.getLogger(__name__)


def _stream_reader(pipe, write_func):
    try:
        with pipe:
            for chunk in iter(lambda: pipe.readline(), b''):
                if not chunk:
                    break
                try:
                    text = chunk.decode('utf-8', 'replace')
                except Exception:
                    text = chunk.decode('latin-1', 'replace')
                write_func(text)
    except Exception:
        LOG.exception("stream reader error")
        

def _downloadOneRegion(region: tuple, out_dir: Path) -> None:
    accession, nucl_range = region
    out_file = (out_dir / accession).with_suffix('.fasta')
    
    acc_downloader_executable = shutil.which('ncbi-acc-download')
    cmd = [acc_downloader_executable,
           '-m', 'nucleotide',
           '-F', 'fasta',
           '-o', '/dev/stdout',
           '-g', nucl_range,
           accession]
    
    with open(out_file, 'w') as handle:
        try:
            subprocess.run(cmd, stdout = handle, check = True, text = True)
        except CalledProcessError:
            LOG.warning(f"Error downloading {accession}")
            return None
        
    with open(out_file, "r") as handle:
        with gzip.open(out_file.with_suffix('.fasta.gz'), "wt") as compressed_handle:
            compressed_handle.writelines(handle)
    os.remove(out_file)
    
    return None


def get_assembly_accessions(scaffolds: list, source: str, no_progress: bool = False) -> list:
    
    # In case of requesting Genbank IDs, redirect WGS records to their master record
    if source == 'Genbank':
        r = re.compile('([A-Z]{4,}[0-9]{8,})\.[0-9]+')
        non_wgs_scaffolds = [sc for sc in scaffolds if not r.findall(sc)]
        wgs_scaffolds = [sc for sc in scaffolds if r.findall(sc)]
        wgs_scaffolds = [wgs_sc[:-11] + '0'*9 for wgs_sc in wgs_scaffolds]
        scaffolds = wgs_scaffolds + non_wgs_scaffolds
        LOG.info(f'Redirected {len(wgs_scaffolds)} Nucleotide WGS IDs to their master records')
    
    # Get Assembly IDs in batches of 100 in max. 3 attempts
    max_attempts = 3
    batch_size = 100
    n_batches = int((len(scaffolds) + 0.1) // batch_size) + 1 # +0.1 to make the hundreds get right
    LOG.info(f'Getting Assembly accession IDs in {n_batches} batches of {batch_size}')
    
    with logging_redirect_tqdm(loggers = [LOG]):
        records = []
        for b_idx, batch in tqdm(list(enumerate(batched(scaffolds, batch_size))), leave = False, disable = no_progress):
            for attempt in range(max_attempts):
                try:
                    with warnings.catch_warnings(action = "ignore"):
                        with Entrez.elink(dbfrom = "nucleotide", db = 'assembly', id = batch) as handle:
                            batch_records = Entrez.read(handle)
                        records += batch_records
                        break
                except:
                    if attempt+1 < max_attempts:
                        LOG.warning(f'Error processing batch {b_idx+1} in attempt {attempt+1}. Max. attempts: {max_attempts}')
                    else:
                        LOG.error(f'Error processing batch {b_idx+1} after {max_attempts} attempts. Skipping...')
        
        # Get UIDs from Elink records
        uids = [str(l['Id']) for record in records for linksetdb in record['LinkSetDb'] for l in linksetdb['Link']]
        
        # Get Esummary objects for every UID in max. 3 attempts
        for attempt in tqdm(range(max_attempts), leave = False, disable = no_progress):
            try:
                with warnings.catch_warnings(action = "ignore"):
                    with Entrez.esummary(db = 'assembly', id = uids) as handle:
                        summary_set = Entrez.read(handle)
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

      