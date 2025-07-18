#!/bin/env python

### Program: cagecleaner
### Author: Lucas De Vrieze
### (c) Masschelein lab
### VIB-KU Leuven Centre of Microbiology
### KU Leuven, Department of Biology

"""
Usage: `cagecleaner -s cblaster_session.json <options>`

cageclear removes redundant hits from the cblaster tool by dereplicating the hosting genomes.

It only requires the cblaster session file of the run you would like to dereplicate. 
You can set the ANI cutoff that skDER uses to dereplicate the genomes.
In addition, it recovers some of the gene cluster diversity lost by the dereplication by assessing gene cluster content and hit score outliers.

This tool will produce seven final output files
    - filtered_session.json: a filtered cblaster session file
    - filtered_binary.txt: a cblaster binary presence/absence table, containing only the retained hits.
    - filtered_summary.txt: a cblaster summary file, containing only the retained hits.
    - clusters.txt: the corresponding cluster IDs from the cblaster summary file for each retained hit.
    - genome_cluster_sizes.txt: the number of genomes in a dereplication genome cluster, referred to by the dereplication representative genome.
    - genome_cluster_status.txt: a table with scaffold IDs, their representative genome assembly and their dereplication status.
    - scaffold_assembly_pairs.txt: a table with scaffold IDs and the IDs of the genome assemblies of which they are part.
    
There are four possible dereplication statuses:
    - 'dereplication_representative': this scaffold is part of the genome assembly that has been selected as the representative of a genome cluster.
    - 'readded_by_content': this scaffold has been kept as it contains a hit that is different in content from the one of the dereplication representative.
    - 'readded_by_score': this scaffold has been kept as it contains a hit that has an outlier cblaster score.
    - 'redundant': this scaffold has not been retained and is therefore removed from the final output.
"""

import pandas as pd
import sys
import subprocess
import re
import os
from itertools import batched
import gzip
import shutil
import argparse
import tempfile
from Bio import SeqIO
from scipy.stats import zscore
from random import choice
from importlib.metadata import version
from copy import deepcopy
from cblaster.classes import Session

__version__ = version("cagecleaner")

ACCESSIONS_SCRIPT = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'get_accessions.sh')
DOWNLOAD_SCRIPT = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'download_assemblies.sh')
DEREPLICATE_SCRIPT = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'dereplicate_assemblies.sh')

global TEMP_DIR
global GENOMES
global SKDER_OUT
global VERBOSE

def parse_arguments():
    """
    Argument parsing function
    """
    
    # Auxiliary function to check the value range of the ANI threshold
    def check_percentage(value):
        if 0 <= float(value) <= 100:
            return value
        else:
            raise argparse.ArgumentTypeError("%s should be a percentage value between 0 and 100" % value)
    
    # Auxiliary function to check whether a path exists
    def check_exists(path):
        if os.path.exists(path):
            return path
        else:
            raise argparse.ArgumentTypeError("%s is not a valid path. Please check the file path." % path)
    
    parser = argparse.ArgumentParser(prog = 'cagecleaner',
                                     epilog = """
    Lucas De Vrieze, 2025
    (c) Masschelein lab, VIB
                                     """,
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description = """
    cagecleaner: A tool to remove redundancy from cblaster hits.
    
    cagecleaner reduces redundancy in cblaster hit sets by dereplicating the genomes containing the hits. 
    It can also recover hits that would have been omitted by this dereplication if they have a different gene cluster content
    or an outlier cblaster score.
    
    cagecleaner first retrieves the assembly accession IDs of each cblaster hit via NCBI Entrez-Direct utilities, 
    then downloads these assemblies using NCBI Datasets CLI, and then dereplicates these assemblies using skDER.
    By default, cblaster hits that have an alternative gene cluster content or an outlier cblaster score 
    (calculated via z-scores) are recovered as well.
                                     """,
                                     add_help = False
                                     )

    args_general = parser.add_argument_group('General')
    args_general.add_argument('-c', '--cores', dest = 'cores', default = 1, type = int, help = "Number of cores to use (default: 1)")    
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('--verbose', dest = 'verbose', default = False, action = 'store_true', help = "Enable verbose logging")
    
    args_io = parser.add_argument_group('Input / Output')
    args_io.add_argument('-s', '--session', dest = "session_file", type = check_exists, help = "Path to cblaster session file")
    args_io.add_argument('-o', '--output', dest = "output_dir", default = '.', help = "Output directory (default: current working directory)")
    args_io.add_argument('-t', '--temp', dest = "temp_dir", default = tempfile.gettempdir(), help = "Path to store temporary files (default: your system's default temporary directory)")
    args_io.add_argument('-exs', '--exclude_scaffolds', dest = 'excluded_scaffolds', default = '', help = "Scaffolds IDs to be excluded from the hit set (comma-separated)")
    args_io.add_argument('-exa', '--exclude_assemblies', dest = 'excluded_assemblies', default = '', help = "Assembly IDs to be excluded from the hit set (comma-seperated)")
    args_io.add_argument('--keep_downloads', dest = "keep_downloads", default = False, action = "store_true", help = "Keep downloaded genomes")
    args_io.add_argument('--keep_dereplication', dest = "keep_dereplication", default = False, action = "store_true", help = "Keep skDER output")
    args_io.add_argument('--keep_intermediate', dest = "keep_intermediate", default = False, action = "store_true", help = "Keep all intermediate data. This overrules other keep flags.")
    
    args_download = parser.add_argument_group('Download')
    args_download.add_argument('--download_batch', dest = 'download_batch', default = 300, type = int, help = "Number of genomes to download in one batch (default: 300)")
    
    args_dereplication = parser.add_argument_group('Dereplication')
    args_dereplication.add_argument('-a', '--ani', dest = 'ani', default = 99.0, type = check_percentage, help = "ANI dereplication threshold (default: 99.0)")
    
    args_recovery = parser.add_argument_group('Hit recovery')
    args_recovery.add_argument('--no_recovery_content', dest = 'no_recovery_by_content', default = False, action = "store_true", help = "Skip recovering hits by cluster content")
    args_recovery.add_argument('--no_recovery_score', dest = 'no_recovery_by_score', default = False, action = "store_true", help = "Skip recovering hits by outlier scores")
    args_recovery.add_argument('--min_z_score', dest = 'zscore_outlier_threshold', default = 2.0, type = float, help = "z-score threshold to consider hits outliers (default: 2.0)")
    args_recovery.add_argument('--min_score_diff', dest = 'minimal_score_difference', default = 0.1, type = float, help = "minimum cblaster score difference between hits to be considered different. Discards outlier hits with a score difference below this threshold. (default: 0.1)")

    args = parser.parse_args()
        
    return args

    
def load_session(path_to_session: str) -> Session:
    """
    This function loads the cblaster session object.
    
    :param str path_to_session: Path to the cblaster session json file
    :rtype dict: A cblaster Session object
    """
    print("--- STEP 0: Loading session file. ---")
    
    session = Session.from_file(path_to_session)
    
    if VERBOSE:
        print("Session loaded!")
    
    return session


def get_scaffolds(session: Session, excluded: list = []) -> list:
    """
    This function extracts the scaffold IDs from a temporary copy of the cblaster binary output file.
    
    :param Session session: cblaster session object
    :rtype list: A list containing Genbank and RefSeq Nucleotide IDs.
    """
    print("\n--- STEP 1: Extracting scaffold IDs. ---")
    
    # Make a temporary copy of the binary table
    if VERBOSE:
        print("Generating binary table from session")
    with open(_temp("binary.txt"), "w") as handle:
        session.format("binary", fp = handle)

    # Read the file using a variable series of spaces as separator and extract the second column (containing the scaffold IDs)
    if VERBOSE:
        print("Reading binary table")
    scaffolds = pd.read_table(_temp("binary.txt"), sep = "\\s{2,}", engine = 'python', usecols = [1])['Scaffold'].to_list()
    
    # Exclude requested scaffolds
    if excluded != ['']:
        if VERBOSE:
            print(f"Excluding scaffolds {', '.join(excluded)}")
        scaffolds = list(set(scaffolds) - set(excluded))

    print(f"Extracted {len(scaffolds)} scaffold IDs")
    
    return scaffolds
        

def get_assemblies(scaffolds: list, excluded: list = []) -> list:
    """
    This function obtains the genome assembly ID for each scaffold ID obtained by get_scaffolds().
    It uses the NCBI Entrez-Direct utilities via a bash subprocess.

    :param list scaffolds: A list containing Genbank and RefSeq Nucleotide IDs.
    :rtype list: A list of Genbank and RefSeq Assembly IDs.
    """        

    print("\n--- STEP 2: Retrieving genome assembly IDs. ---")
    
    # save the scaffold list in a temporary file
    if VERBOSE:
        print("Writing list of scaffolds")
    with open(_temp('scaffolds.txt'), 'w') as handle:
        handle.writelines('\n'.join(scaffolds))
    
    # map to assembly IDs using E-utilities
    if VERBOSE:
        print("Calling NCBI Entrez Utilities")
    home = os.getcwd()
    os.chdir(ASSIGNED_TEMP_DIR)
    subprocess.run(['bash', ACCESSIONS_SCRIPT, 'scaffolds.txt'], check = True)
    os.chdir(home)
    
    # read the result file
    if VERBOSE:
        print("Reading results")
    with open(_temp('assembly_accessions'), 'r') as handle:
        assemblies = [l.rstrip() for l in handle.readlines()]
        
    # Excluded requested assemblies
    if excluded != ['']:
        if VERBOSE:
            print(f"Excluding assemblies {', '.join(excluded)}")
        assemblies = list(set(assemblies) - set(excluded))
    
    return assemblies

    
def download_genomes(assemblies: list, batch_size: int = 300) -> None:
    """
    This function downloads the full nucleotide fasta files of all found assemblies using the NCBI Datasets CLI via a bash subprocess.
    It automatically selects for the most recent accession version by omitting the version digit and relying on the NCBI Datasets CLI defaults 
    to download the most recent version.
    The assemblies are downloaded in batches of 300 by default and saved in the temporary folder data/genomes in the working directory.
    
    :param list assemblies: A list containing Genbank and RefSeq assembly IDs.
    :param int batch_size: The number of assemblies to download per batch. [default: 300]
    """
    
    print("\n--- STEP 3: Downloading genomes. ---")
    
    # Cut off the version digits
    versionless_assemblies = [acc.split('.')[0] for acc in assemblies]
    
    # Prepare the batches and save them in a temporary file
    if VERBOSE:
        print('Preparing download batches')
    with open(_temp('download_batches.txt'), 'w') as file:
        download_batches = list(batched(versionless_assemblies, batch_size))
        for batch in download_batches:
            file.write(' '.join(batch) + '\n')

    # Run the bash script to download genomes:
    if VERBOSE:
        print('Calling NCBI Datasets CLI')
    home = os.getcwd()
    os.chdir(ASSIGNED_TEMP_DIR)
    subprocess.run(["bash", DOWNLOAD_SCRIPT], check=True)
    os.chdir(home)
    
    return None


def map_scaffolds_to_assemblies(scaffolds: list, assemblies: list) -> dict:
    """
    This function maps the scaffolds in the list of cblaster hits to the genome assembly they are part of. To get the exact scaffold IDs of each assembly,
    scaffolds are retrieved from the headers of the assembly fasta files. Prefixes are split off from both scaffold ID sets (the one from cblaster,
    and the one from the NCBI download) during mapping as scaffold IDs mapped by the E-utilities sometimes do not correspond with the ones
    in the downloaded fasta files.
    It generates a report text file with the scaffold-assembly pairs.
        
    :param list scaffolds: A list containing Genbank and RefSeq Nucleotide IDs.
    :param assemblies: A list containing Genbank and RefSeq Assembly IDs.
    :rtype dict: A dictionary with scaffold Nucleotide IDs as keys and a containing assembly ID as value
    """
    
    print("\n--- STEP 4: Mapping local scaffold IDs to assembly IDs. ---")
    
    # Auxiliary function to split off prefixes from Nucleotide accession IDs
    def split_off_prefix(txt: str) -> str:
        if '_' not in txt:
            return txt
        else:
            return '_'.join(txt.split('_')[1:])
        
    # Auxiliary function to link the Nucleotide accession ID without prefix back to the original ID
    def map_back(no_prefix: str, scaffolds: list) -> str:
        return [s for s in scaffolds if no_prefix in s][0]
    
    scaffolds_set = set([split_off_prefix(s) for s in scaffolds]) # scaffolds from cblaster, deprefixed
    mappings = {}
    
    for assmbl in assemblies:
        if VERBOSE:
            print(f"Mapping {assmbl}")
        
        # Find the path to the downloaded fasta file for this assembly ID, ignoring version digits
        try:
            assmbl_file = [a for a in os.listdir(GENOMES) if assmbl.split('.')[0] in a][0]
        except IndexError:
            print(f'No assembly file found for {assmbl}!')
            continue
        
        # All genomes were gzip-compressed
        with gzip.open(os.path.join(GENOMES, assmbl_file), "rt") as genome:
            try:
                # get all scaffold Nucleotide accession IDs in this assembly
                scaffolds_in_this_assembly = [record.id for record in SeqIO.parse(genome, 'fasta')]
                
                # split off any prefix
                scaffolds_in_this_assembly = [split_off_prefix(i) for i in scaffolds_in_this_assembly]
                
                # find the ones we have a cblaster hit for
                scaffolds_in_this_assembly = set(scaffolds_in_this_assembly)
                found_scaffolds_no_prefix = list(scaffolds_set.intersection(scaffolds_in_this_assembly))
                
                # map the deprefixed IDs back to the original ones
                found_scaffolds = [map_back(s, scaffolds) for s in found_scaffolds_no_prefix]
                
                if VERBOSE:
                    print(f"Found scaffolds {', '.join(found_scaffolds)}")
                
                # add a mapping item for all scaffolds with a cblaster hit
                for scaff in found_scaffolds:
                    mappings[scaff] = assmbl
                    
            except IndexError:
                print(f'No corresponding scaffold accession ID could be found for {assmbl}!')
    
    # Write scaffold-assembly pairs to report file
    if VERBOSE:
        print("Writing mapping pairs to file")
    with open('scaffold_assembly_pairs.txt', 'w') as handle:
        handle.write('\n'.join([s + '\t' + a for s,a in mappings.items()]))
    
    print(f"Found {len(mappings)} scaffold-assembly links")
    
    return mappings


def dereplicate_genomes(ani_threshold: float = 99.0, nb_cores: int = 1) -> None:
    """
    This function takes a list of assembly IDs and calls a helper bash script that dereplicates the genomes using skDER.
    The default ANI cutoff for dereplication is 99 %.
    
    A file called "dereplicated_assemblies.txt" is then generated by the helper script, in addition to skDER's output which can be found
    in the temporary folder data/skder_out in the working directory.'

    :param float ani_threshold: The percent identity cutoff for dereplicating genomes (see skDER docs) [default: 99.0]
    :param int nb_cores: The number of cores to be available for dereplication [default: 1]
    """
    print("\n--- STEP 5: Dereplicating genomes. ---")
    
    if VERBOSE:
        print("Calling skDER")
    home = os.getcwd()
    os.chdir(ASSIGNED_TEMP_DIR)
    subprocess.run(['bash', DEREPLICATE_SCRIPT, str(ani_threshold), str(nb_cores)], check = True)
    os.chdir(home)
    
    return None


def parse_dereplication_clusters(scaffold_assembly_pairs: dict) -> pd.DataFrame:
    """
    This function parses the secondary clustering result file from skDER to characterise the size and the members of the genome clusters.
    It lists which genomes were selected as representative and to which genome cluster all assemblies belong.
    In addition, it produces a report file with the size of the genome clusters.
    
    :param dict scaffold_assembly_pairs: A dictionary mapping scaffold Nucleotide IDs to their overarching Assembly ID.
    :rtype Pandas dataframe: A dataframe with Nucleotide IDs in the index and columns 'representative' and 'dereplication_status',
                             which, resp., refer to the dereplication representative of that assembly,
                             and the dereplication status of that assembly ('dereplication_representative' or 'redundant').
    """
    print("\n--- STEP 6: Parsing dereplication genome clusters. ---")
    
    # Auxiliary function to extract assembly accession IDs from the file paths listed in the skDER output table
    def extract_assembly(path: str) -> str:
        pattern = re.compile("GC[A|F]_[0-9]{9}\\.[1-9]+")
        assembly = pattern.findall(path)[0]
        return assembly
    
    # Auxiliary function to rename the column names parsed from the skDER output
    def remap_type_label(label: str) -> str:
        mapping = {'representative_to_self': 'dereplication_representative',
                   'within_cutoffs_requested': 'redundant'}
        return mapping[label]
    
    # Auxiliary function to retrieve all scaffold Nucleotide IDs with a cblaster hit for a certain assembly ID
    def map_scaffolds(assembly: str, scaffold_assembly_pairs: dict) -> list:
        return [s for s,a in scaffold_assembly_pairs.items() if a == assembly]
    
    if VERBOSE:
        print("Reading skDER table")
    # Parse the skDER output table. Rename the columns and extract Assembly IDs on-the-fly.
    genome_clusters_df = pd.read_table(os.path.join(SKDER_OUT, "skDER_Clustering.txt"), 
                                       converters={'assembly': extract_assembly,
                                                   'representative': extract_assembly,
                                                   'dereplication_status': remap_type_label},
                                       names = ['assembly', 'representative', 'dereplication_status'],
                                       usecols = [0,1,4], header = 0, index_col = 'assembly'
                                       ).sort_values(by = ['representative','dereplication_status'])
    
    ## Replace Assembly accession IDs by Nucleotide IDs.
    ## If there are multiple Nucleotide IDs with a cblaster hit, the records will be replicated.
    
    # Deconstruct the initial dataframe from the parsing into a dictionary with Assembly accessions as keys and representative assemblies as values
    genome_clusters_records = genome_clusters_df.to_dict(orient = 'index')
    
    # Build a similar dictionary using Nucleotide accessions as keys, replicating representative assemblies if needed
    if VERBOSE:
        print("Building scaffold-representative assembly table")
    genome_clusters = {}
    for assembly, representative_status in genome_clusters_records.items():
        mapped_scaffolds = map_scaffolds(assembly, scaffold_assembly_pairs) # map assembly to a list of scaffold(s) with a cblaster hit
        for scaffold in mapped_scaffolds:
            if VERBOSE:
                print(f'Mapping scaffold {scaffold} from assembly {assembly} to representative {representative_status["representative"]}')
            genome_clusters[scaffold] = representative_status
    
    # Construct a new dataframe from this Nucleotide-keyed dictionary
    genome_clusters = pd.DataFrame.from_dict(genome_clusters, orient = 'index', columns = ['representative', 'dereplication_status'])
    genome_clusters = genome_clusters.sort_values(by = ['representative','dereplication_status'])
    
    # Determine the size of each genome cluster and write report file
    clust_size = pd.DataFrame(genome_clusters.groupby(by = "representative")['representative'].count()).rename(columns = {'representative': 'size'})
    clust_size.to_csv('genome_cluster_sizes.txt', sep = "\t")
    print('Genome cluster sizes written in genome_cluster_sizes.txt')
    
    return genome_clusters


def recover_hits(session: Session, genome_clusters_mapping: pd.DataFrame, not_by_content: bool = False, not_by_score: bool = False, outlier_z: float = 2, min_score_diff: float = 0.1) -> pd.DataFrame:
    """
    This function recovers some variation in gene clusters that was lost due to dereplicating the hosting genomes. It offers two approaches to
    flag interesting gene clusters that will be kept in the output.
    
    1) Different gene cluster content
    A scaffold that is part of a dereplication cluster may have a different gene cluster content than its representative, i.e. a different number of
    identified homologs. The cblaster binary table also lists the number of homologs were found in a gene cluster for each query gene.
    If these numbers are different from the ones of the representative assembly, flag to keep this scaffold and its hit.
    
    2) Outlier cblaster score
    cblaster adds a 'Score' column in the binary output table that captures the numbers of homologs as well as aggregates the level of homology
    of the entire cluster. If this score is significantly different from the other scores, this may indicate multiple gene cluster lineages within
    this genome cluster, or a remarkably fast evolving gene cluster. In any case, it is interesting to retain this hit.
    Signficance is currently determined using z-scores for each cluster content group.
    
    :param Session session: cblaster session object
    :param pandas Dataframe genome_clusters_mapping: dataframe returned by the parse_dereplication_clusters() function, containing the dereplication
                                                     status and representative of each assembly
    :param bool not_by_content: flag to disable recovering gene clusters by gene cluster content, also disables recovery by outlier score [default: False]
    :param bool not_by_score: flag to disable recovering gene clusters by outlier cblaster score [default: False]
    :param float outlier_z: minimum absolute value of the z-score to consider a hit as an outlier
    :param float min_score_diff: minimum difference in cblaster score between a hit and the mode score when determining outlier scores.
                                 Outliers with a score difference below this value are discarded. [default: 0.1]
    :rtype Pandas dataframe: updated dereplication status table, now also containing flags to keep scaffolds and their gene clusters in the final output
    """
    print("\n--- STEP 7: Recovering gene cluster diversity. ---")
    
    # Make a copy of the original genome cluster table
    updated_mapping = genome_clusters_mapping.copy()
    
    # Count the number of recovered scaffolds
    recovered = 0
    
    # Skip this if not recovering by gene cluster content
    if not(not_by_content):
        
        # Read the relevant columns from a temporary copy of the cblaster binary table
        if VERBOSE:
            print("Generating binary table from session file")
        with open(_temp("binary_recovery.txt"), "w") as handle:
            session.format("binary", fp = handle)
        hits = pd.read_table(_temp("binary_recovery.txt"), sep = "\\s{2,}", engine = 'python', 
                             usecols = lambda x: x not in ['Organism', 'Start', 'End'],
                             index_col = "Scaffold")
        
        ## Hits are recovered within dereplication clusters so we will check the hits of each dereplication cluster
        # Get the IDs of all assemblies grouped by their dereplication representative
        if VERBOSE:
            print("Grouping hits by representative assembly")
        grouped_mapping = genome_clusters_mapping.reset_index(names = 'scaffold').groupby('representative').agg(list)
        genome_groups = dict(zip(grouped_mapping.index, grouped_mapping['scaffold']))
        
        for representative, group in genome_groups.items():

            # Get the cblaster hit data for scaffolds of this dereplication cluster
            hits_this_group = hits.loc[group]
            
            if VERBOSE:
                print(f"Grouping hits for representative {representative} by content")
            # Split the dereplication grouping further into gene cluster subgroups as reflected by the number of homologs of each query gene
            hits_this_group_by_content_group = [l.to_list() for l in hits_this_group.groupby(
                                                list(set(hits_this_group.columns).difference({'Score'})) # to get all query gene columns
                                                ).groups.values()]
            
            # Loop over all structural subgroups
            if VERBOSE:
                print("Checking all structural subgroups")
            for content_group in hits_this_group_by_content_group:
                scores_this_content_group = hits_this_group.loc[content_group, 'Score'] # cblaster scores
                mode_score_this_content_group = float(scores_this_content_group.mode().iloc[0]) # modal cblaster score
                zscores_this_content_group = scores_this_content_group.apply(zscore, by_row = False) # zscores
                
                # If the result of the z-score calculation yields all NaNs, all cblaster scores were identical, implying there is no score outlier.
                no_outlier_screening_needed = zscores_this_content_group.isna().all()
                if VERBOSE and no_outlier_screening_needed:
                    print("No outliers in this subgroup!")
                
                # In case there are not score outliers, we can continue flagging different gene cluster contents, if any
                # If the user is not interested in flagging score outlier hits, we end up in this case anyway.
                if not_by_score or no_outlier_screening_needed:
                    
                    # If the dereplication representative is one of the assemblies that has a scaffold in this structural subgroup, 
                    # then we already have a hit from this subgroup, so we can skip this one.
                    dereplication_status_this_content_group = [genome_clusters_mapping.loc[m, 'dereplication_status'] for m in content_group]
                    if "dereplication_representative" in dereplication_status_this_content_group:
                        if VERBOSE:
                            print(f"Skipped representative group {representative}")
                        continue
                    
                    # If the dereplication representative is not in this structural subgroup, randomly pick a member as a representative hit to keep
                    else:
                        content_group_representative = choice(content_group)
                        if VERBOSE:
                            print(f"Picked {content_group_representative} as representative for content subgroup")
                        updated_mapping.at[content_group_representative, 'dereplication_status'] = 'readded_by_cluster_content'
                        recovered += 1
                
                # The case there were different cblaster scores and the user is interested in score outliers
                else:
                    
                    # Select hits that have an absolute z-score above the threshold,
                    # and of which the cblaster score is sufficiently different from the modal score,
                    # to avoid 'false' outliers that are just a little bit more different than most hits
                    if VERBOSE:
                        print(f'Calculating z scores for representative group {representative}')
                    outliers_this_content_group = zscores_this_content_group.loc[
                        (zscores_this_content_group.abs() >= outlier_z) &
                        ((scores_this_content_group - mode_score_this_content_group).abs() >= min_score_diff)].index
                    
                    # All score outliers will be added
                    if VERBOSE:
                        print("Checking all outliers")
                    for outlier in outliers_this_content_group.to_list():
                        # If the overarching assembly of this outlier is the dereplication representative, skip flagging it
                        dereplication_status_outlier = genome_clusters_mapping.loc[outlier, 'dereplication_status']
                        if dereplication_status_outlier == "dereplication_representative":
                            if VERBOSE:
                                print(f"Skipped representative {representative}")
                            continue
                        else:
                            if VERBOSE:
                                print(f'Recovered {outlier} as score outlier')
                            updated_mapping.at[outlier, 'dereplication_status'] = 'readded_by_outlier_score'
                            recovered += 1
                            
                    # We still have to flag a non-outlier representative of this structural subgroup
                    non_outliers_this_content_group = zscores_this_content_group.drop(index = outliers_this_content_group).index.to_list()
                    
                    # If this subgroup contains a scaffold of the dereplication representative, skip this subgroup
                    dereplication_status_non_outliers_this_content_group = [genome_clusters_mapping.loc[m, 'dereplication_status'] 
                                                                              for m in non_outliers_this_content_group]
                    if "dereplication_representative" in dereplication_status_non_outliers_this_content_group:
                        if VERBOSE:
                            print(f"Skipped representative group {representative}")
                        continue
                    
                    # If not, pick a random representative from the non-outlier hits
                    else:
                        content_group_representative = choice(non_outliers_this_content_group)
                        if VERBOSE:
                            print(f"Picked {content_group_representative} as representative for content subgroup")
                        updated_mapping.at[content_group_representative, 'dereplication_status'] = 'readded_by_cluster_content'
                        recovered += 1
    
    # There is nothing to recover in this case
    else:
        print('All revisiting options have been disabled. Not updating the genome grouping labels.')
    
    # Tidy up the updated mapping
    updated_mapping = updated_mapping.sort_values(by = ['representative', 'dereplication_status'])
    
    # Write report file
    if VERBOSE:
        print("Writing dereplication status table")
    updated_mapping.reset_index(names = 'scaffold').to_csv('genome_cluster_status.txt', sep = "\t", index = False)
    
    print(f"Recovered {recovered} scaffolds")
    
    return updated_mapping

                            
def get_dereplicated_scaffolds(genome_clusters: pd.DataFrame) -> list:
    """
    This function retrieves the final retained scaffold IDs, after dereplication and hit recovery.

    :param Pandas dataframe; A dataframe containing the final status of each scaffold and its dereplication representative
    :rtype list: A list containing the retained scaffolds
    """  
    print("\n--- STEP 8: Gathering retained scaffold IDs. ---")
    
    dereplicated_scaffolds = genome_clusters[genome_clusters['dereplication_status'] != 'redundant'].index.to_list()

    print(f"Got {len(dereplicated_scaffolds)} representative scaffolds")
    
    return dereplicated_scaffolds


def generate_output(dereplicated_scaffolds:list, session:Session) -> None:
    """
    This function takes a list of retained scaffold IDs and the cblaster session object.
    From this object, it generates the filtered session, binary and summary file,
    as well as a list of cluster number of the retained gene clusters.

    :param list dereplicated_scaffolds: A list containing the retained scaffold IDs.
    :param Session session: cblaster session object of the unfiltered session
    """
    print("\n--- STEP 9: Generating output files. ---")
    
    # Get a dictionary export of the session object
    session_dict = session.to_dict()
    
    # Make a deep copy to start carving out the new session
    if VERBOSE:
        print("Parsing original session file")
    filtered_session_dict = deepcopy(session_dict)
    
    ## Filtering the json session file
    # Remove all cluster hits that do not link with a dereplicated scaffold
    # Remove strains that have no hits left
    if VERBOSE:
        print("Filtering session file")
    for org_idx, org in reversed(list(enumerate(session_dict['organisms']))):
        for hit_idx, hit in reversed(list(enumerate(org['scaffolds']))):
            if hit['accession'] not in dereplicated_scaffolds:
                filtered_session_dict['organisms'][org_idx]['scaffolds'].pop(hit_idx)
                if not filtered_session_dict['organisms'][org_idx]['scaffolds']:
                    filtered_session_dict['organisms'].pop(org_idx)
                
    # Create the filtered cblaster session object
    if VERBOSE:
        print("Creating filtered cblaster session")
    filtered_session = Session.from_dict(filtered_session_dict)
    
    # Generate the outputs
    # session file
    if VERBOSE:
        print("Writing filtered session file")
    with open("filtered_session.json", "w") as session_handle:
        filtered_session.to_json(fp = session_handle)
        
    # binary table
    if VERBOSE:
        print("Writing filtered binary table")
    with open("filtered_binary.txt", 'w') as filtered_binary_handle:
        filtered_session.format(form = "binary", fp = filtered_binary_handle)
        
    # summary file
    if VERBOSE:
        print('Writing filtered summary file')
    with open("filtered_summary.txt", "w") as filtered_summary_handle:
        filtered_session.format(form = "summary", fp = filtered_summary_handle)    
        
    # list of cluster numbers
    if VERBOSE:
        print("Saving list of retained cluster numbers")
    filtered_cluster_numbers = [cluster['number'] 
                                for organism in filtered_session_dict['organisms'] 
                                for scaffold in organism['scaffolds'] 
                                for cluster in scaffold['clusters']]
    with open("clusters.txt", "w") as numbers_handle:
        numbers_handle.write(','.join([str(nb) for nb in filtered_cluster_numbers]))
            
    print("Output generated!")
            
    return None

def _temp(folder) -> os.path:
    """
    Internal auxiliary function to map a filename to the temporary workfolder.
    """
    if isinstance(folder, str):
        path = [folder]
    else:
        path = folder
    return os.path.join(ASSIGNED_TEMP_DIR, *path)
    

def main():
    """
    Main function
    """
    global ASSIGNED_TEMP_DIR
    global GENOMES
    global SKDER_OUT
    global VERBOSE
    
    ## Parse arguments
    args = parse_arguments()

    work_dir = args.output_dir.rstrip('/')
    temp_dir = args.temp_dir.rstrip('/')
    path_to_session = os.path.join(os.getcwd(), args.session_file)
    nb_cores = int(args.cores)
    VERBOSE = args.verbose
    ani_threshold = float(args.ani)
    no_recovery_by_content = args.no_recovery_by_content
    no_recovery_by_score = args.no_recovery_by_score
    outlier_z = float(args.zscore_outlier_threshold)
    min_score_diff = float(args.minimal_score_difference)
    download_batch_size = int(args.download_batch)
    keep_downloads = args.keep_downloads
    keep_intermediate = args.keep_intermediate
    keep_dereplication = args.keep_dereplication
    excluded_scaffolds = [i.strip() for i in args.excluded_scaffolds.split(',')]
    excluded_assemblies = [i.strip() for i in args.excluded_assemblies.split(',')]

    ## Set up working and temporary directory
    os.makedirs(work_dir, exist_ok = True)
    with tempfile.TemporaryDirectory(dir = temp_dir) as ASSIGNED_TEMP_DIR:
        os.chdir(work_dir)
        GENOMES = _temp(['data', 'genomes'])
        SKDER_OUT = _temp(['data', 'skder_out'])
        
        ## Execute the workflow
        # load session file and generate cluster tables
        session = load_session(path_to_session)
        # extract scaffold IDs from the cblaster output
        scaffolds = get_scaffolds(session, excluded_scaffolds)
        # link to assembly IDs via the NCBI E-utilities
        assemblies = get_assemblies(scaffolds, excluded_assemblies)
        # download assemblies using the NCBI Datasets CLI
        download_genomes(assemblies, download_batch_size)
        # map the downloaded scaffold IDs to assembly IDs
        scaffold_assembly_pairs = map_scaffolds_to_assemblies(scaffolds, assemblies)
        # dereplicate the genomes using skDER
        dereplicate_genomes(ani_threshold, nb_cores)
        # parse the secondary clustering from the skDER output and construct a dereplication status table
        genome_cluster_composition = parse_dereplication_clusters(scaffold_assembly_pairs)
        # recover gene cluster hits and update status table
        updated_status = recover_hits(session, genome_cluster_composition, no_recovery_by_content, no_recovery_by_score, outlier_z, min_score_diff)
        # retrieve the finally retained scaffold IDs from the updated status table
        dereplicated_scaffolds = get_dereplicated_scaffolds(updated_status)
        # generate final output files
        generate_output(dereplicated_scaffolds, session)
        
        ## Finish by copying over intermediate results if flagged, and removing the temporary folder
        if keep_intermediate or keep_downloads:
            if VERBOSE:
                print('Copying downloaded genomes to output folder')
            shutil.copytree(GENOMES, os.path.join('data', 'genomes'), dirs_exist_ok = True)
        if keep_intermediate or keep_dereplication:
            if VERBOSE:
                print("Copying skDER results to output folder")
            shutil.copytree(SKDER_OUT, os.path.join('data', 'skder_out'), dirs_exist_ok = True)
            
    print(f"\nAll done! Results can be found in {work_dir}")


if __name__ == "__main__":
    main()
    