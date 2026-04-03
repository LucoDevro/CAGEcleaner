#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import pandas as pd
import networkx as nx
from copy import deepcopy


LOG = logging.getLogger(__name__)


def correctLayouts(binary_df: pd.DataFrame) -> pd.DataFrame:
    # Calculate the complementary layout for each cluster
    original = list(zip(binary_df['Strand'],
                        binary_df['Layout_group']))
    complementary = list(zip(binary_df['Strand'].apply(lambda x: tuple([-i for i in x])),
                             binary_df['Layout_group'].apply(lambda x: tuple(reversed(x)))))
    orig_compl_pairs = list(zip(original, complementary))
    
    # Break backloops using a direct graph
    G = nx.DiGraph()
    G.add_edges_from(list(set(orig_compl_pairs)))
    G_pruned = G.copy()
    for u,v in G.edges():
        if G_pruned.has_edge(u,v) and G_pruned.has_edge(v,u):
            G_pruned.remove_edge(u,v)
    backloop_edges = list(G_pruned.edges())
    
    # Now correct complementary layouts
    correction_mapping = dict(backloop_edges)
    corrected = deepcopy(original)
    for orig_idx, orig in enumerate(corrected):
        if orig in correction_mapping.keys():
            corrected[orig_idx] = correction_mapping[orig]
    corrected_layouts = [ly for _,ly in corrected]
    binary_df['Layout_group'] = corrected_layouts
    
    return binary_df

