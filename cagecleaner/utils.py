#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import re
import shutil
import threading
import subprocess
import pandas as pd
import networkx as nx
from copy import deepcopy
from cblaster.classes import Session
from pathlib import Path


LOG = logging.getLogger()


def correctLayouts(binary_df: pd.DataFrame) -> pd.DataFrame:
    # Calculate the complementary layout for each cluster
    original = list(zip(binary_df['Strand'],
                        binary_df['Layout_group']))
    complementary = list(zip(binary_df['Strand'].apply(lambda x: tuple([-i for i in x])),
                             binary_df['Layout_group'].apply(lambda x: tuple(reversed(x)))))
    orig_compl_pairs = list(zip(original, complementary))
    
    # Break backloops using a directed graph
    G = nx.DiGraph()
    G.add_edges_from(list(set(orig_compl_pairs)))
    G_pruned = G.copy()
    for u,v in G.edges():
        if G_pruned.has_edge(u,v) and G_pruned.has_edge(v,u):
            G_pruned.remove_edge(u,v)
    backloop_edges = list(G_pruned.edges())
    
    # Now flip complementary layouts
    correction_mapping = dict(backloop_edges)
    corrected = deepcopy(original)
    for orig_idx, orig in enumerate(corrected):
        if orig in correction_mapping.keys():
            corrected[orig_idx] = correction_mapping[orig]
    corrected_layouts = [ly for _,ly in corrected]
    binary_df['Layout_group'] = corrected_layouts
    
    return binary_df


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
        

def run_command(cmd_list: list, max_attempts: int = 3) -> None:
    executable = Path(shutil.which(cmd_list[0]))
    cmd = [str(executable)] + cmd_list[1:]
    def command_stdout_log(s): return LOG.debug(s.rstrip())
    def command_stderr_log(s): return LOG.warning(s.rstrip())
    
    LOG.debug(f'Running command: {" ".join(cmd)}')
    
    for attempt in range(max_attempts):
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        t_out = threading.Thread(target=_stream_reader, args=(proc.stdout, command_stdout_log))
        t_err = threading.Thread(target=_stream_reader, args=(proc.stderr, command_stderr_log))
        t_out.daemon = True
        t_err.daemon = True
        t_out.start()
        t_err.start()
    
        returncode = proc.wait()
        t_out.join()
        t_err.join()
    
        if returncode != 0:
            LOG.error(f"{executable.name} failed with code {returncode}.")
            if attempt < max_attempts:
                LOG.error(f'Retrying. Attempt {attempt+1} out {max_attempts}.')
        else:
            LOG.info(f'{executable.name} finished successfully')
            break
            
    return None


def generate_cblaster_session(tables_folder: Path, mode: str) -> Session:
    
    hits_df = pd.read_table(tables_folder / 'hits.tsv', sep = "\t",
                            usecols = ['db_id', 'query', 'scaff', 'strand', 'coords', 'evalue', 'score', 'seqid', 'tcov'])
    clusters_df = pd.read_table(tables_folder / 'clusters.tsv', sep = "\t",
                                usecols = ['number', 'hits', 'start', 'end', 'length', 'score', 'scaff', 'taxon_name', 'taxon_id'])
    queries_df = pd.read_table(tables_folder / 'queries.tsv', sep = "\t")
    
    session_dict = {}

    ### Queries field
    session_dict['queries'] = queries_df['id'].to_list()

    ### Query field
    cblaster_query = {}
    cblaster_query['indices'] = []
    ## Query Subjects field
    cblaster_query_subjects = []
    previous_end = -500
    for row in queries_df.itertuples(index = False):
        this_subject = {}
        this_subject['id'] = None
        this_subject['hits'] = []
        this_subject['name'] = row.id
        this_subject['ipg'] = None
        
        # Sequence length is determined from the number of modelled residues in the CIF file
        # cblaster session files take a margin of 500 aa between the hits
        length = row.end - row.start
        this_subject['start'] = previous_end + 500
        this_subject['end'] = this_subject['start'] + length
        previous_end = this_subject['end']
        
        this_subject['strand'] = 1
        this_subject['sequence'] = ""
        
        cblaster_query_subjects.append(this_subject)
        
    cblaster_query['subjects'] = cblaster_query_subjects
    cblaster_query['intermediate_genes'] = []
    cblaster_query['score'] = 0
    cblaster_query['start'] = 0
    cblaster_query['end'] = cblaster_query_subjects[-1]['end']
    cblaster_query['number'] = 0

    session_dict['query'] = cblaster_query

    ### Params field
    cblaster_params = {}
    cblaster_params['mode'] = mode
    cblaster_params['database'] = []
    cblaster_params['min_identity'] = 0
    cblaster_params['min_coverage'] = 0
    cblaster_params['max_evalue'] = 1
    cblaster_params['require'] = []
    cblaster_params['query_file'] = None
    cblaster_params['rid'] = None
    cblaster_params['entrez_query'] = ""

    session_dict['params'] = cblaster_params

    ### Organisms field
    cblaster_organisms = {}

    ## Group the session's cluster records in the same way as they will be structured in the session file,
    ## which will make processing much easier

    # Group by taxon ID and scaffold ID
    grouped_cl_df = clusters_df.groupby(['taxon_id', 'taxon_name', 'scaff'])

    ## Make the cblaster session fields inside out, i.e. populate organisms first with scaffolds (and other attributes),
    ## then populate the scaffolds with clusters and subjects, then populate the clusters with links to the subjects.
    for (txid, txname, scaff), clusters in grouped_cl_df:
        # Create a new organism instance if there's no one for this taxon ID
        if (txid, txname) in cblaster_organisms.keys():
            this_organism = cblaster_organisms[(txid, txname)]
        else:
            this_organism = {'name': txname,
                              'strain': "",
                              'scaffolds': {}}
            cblaster_organisms[(txid, txname)] = this_organism
        
        # Create a new scaffold instance of there is no one for this scaffold ID in this taxon
        if scaff in this_organism['scaffolds'].keys():
            this_scaffold = this_organism['scaffolds'][scaff]
        else:
            this_scaffold = {'accession': scaff,
                              'subjects': [],
                              'clusters': []}
            this_organism['scaffolds'][scaff] = this_scaffold
            
        # Populate this scaffold with clusters for each cluster identified on this scaffold
        # Keep track of the number of hits covered on this scaffold for proper references in the subject links
        nb_hits_covered = 0
        for cl in clusters.itertuples(index = False):
            cblaster_this_cluster = {}
            cblaster_this_cluster['subjects'] = []
            cblaster_this_cluster['intermediate_genes'] = []
            cblaster_this_cluster['score'] = cl.score
            cblaster_this_cluster['start'] = cl.start
            cblaster_this_cluster['end'] = cl.end
            cblaster_this_cluster['number'] = cl.number
            this_scaffold['clusters'].append(cblaster_this_cluster)
            
            # Populate this scaffold with subjects (hits) for each hit identified on this scaffold
            # and add the link with the cluster
            these_hits = cl.hits.split(',')
            hits_covered_this_cluster = len(these_hits)
            cblaster_this_cluster['indices'] = list(range(nb_hits_covered, nb_hits_covered + hits_covered_this_cluster))
            nb_hits_covered += hits_covered_this_cluster
            
            # First retrieve information for the hits from the hit table
            these_hits_df = hits_df[hits_df['db_id'].isin(these_hits)]
            these_hits_df = these_hits_df[these_hits_df['scaff'] == cl.scaff]
            these_hits_df = these_hits_df[[all([int(cdr) in pd.Interval(left = cl.start, 
                                                                        right = cl.end, 
                                                                        closed = 'both')
                                                for cdr in re.findall(r'\d+', crds)]
                                               )
                                           for crds in these_hits_df['coords']
                                           ]]
            
            for hit in these_hits_df.itertuples():
                cblaster_this_subject = {}
                cblaster_this_subject['id'] = None
                cblaster_this_subject['hits'] = [{'query': hit.query,
                                                  'subject': hit.db_id,
                                                  'identity': hit.seqid,
                                                  'coverage': hit.tcov,
                                                  'evalue': hit.evalue,
                                                  'bitscore': hit.score}]
                cblaster_this_subject['name'] = hit.db_id
                cblaster_this_subject['ipg'] = hit.db_id
                cblaster_this_subject['start'] = min([int(cdr) for cdr in re.findall(r'\d+', hit.coords)])
                cblaster_this_subject['end'] = max([int(cdr) for cdr in re.findall(r'\d+', hit.coords)])
                cblaster_this_subject['strand'] = int(f"{hit.strand}1")
                cblaster_this_subject['sequence'] = None
                this_scaffold['subjects'].append(cblaster_this_subject)

    ## Discard the taxon ID and scaffold ID nested index to get the cblaster session format
    cblaster_organisms = list(cblaster_organisms.values())
    cblaster_organisms_new = []
    for organism in cblaster_organisms:
        organism['scaffolds'] = list(organism['scaffolds'].values())
        cblaster_organisms_new.append(organism)
        
    session_dict['organisms'] = cblaster_organisms_new

    ### Construct the cblaster Session
    session = Session.from_dict(session_dict)
    
    return session

