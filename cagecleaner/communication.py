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


LOG = logging.getLogger()
warnings.filterwarnings(action = "ignore", module = "Bio")


def _downloadOneRegion(region: tuple, out_dir: Path) -> None:
    accession, nucl_range = region
    out_file = Path('|'.join([str(out_dir / accession)] + nucl_range.split(':')) + '.fasta')
    max_attempts = 3
    
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


def get_assembly_accessions(scaffolds: list, source: str, no_progress: bool = False) -> list:
    
    # Only in case of requesting Genbank IDs, redirect WGS records to their master record
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
    
    # Get Esummary objects for every UID in max. 3 attempts
    LOG.debug('Retrieving Assembly IDs from Elink responses')
    for attempt in range(max_attempts):
        try:
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

      