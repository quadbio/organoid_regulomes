#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc

from matching import get_cost_knn_graph, mcmf


# FUNC
def interface():
    parser = argparse.ArgumentParser(description='Extracts UMIs supported by enough reads using a GMM.')

    parser.add_argument('H5AD',
                    type=str,
                    metavar='<f>',
                    help='H5AD file with cell embeddings.')

    parser.add_argument('-o', '--out',
                    dest='out',
                    type=str,
                    default='./matches.tsv',
                    metavar='<f>',
                    help='Output file with cell to cell matches.')

    parser.add_argument('-l', '--layer',
                    dest='layer',
                    type=str,
                    default='matrix',
                    metavar='<str>',
                    help='Layer in the H5AD file.')

    parser.add_argument('-s', '--split-col',
                    dest='split_col',
                    type=str,
                    default='modality',
                    metavar='<str>',
                    help='Metadata column to split data.')

    parser.add_argument('-k',
                    dest='k',
                    type=int,
                    default='10',
                    metavar='<int>',
                    help='Number of neighbors for kNN.')

    parser.add_argument('-p', '--penalty',
                    dest='penalty',
                    type=int,
                    default='99',
                    metavar='<int>',
                    help='Penalty percentile for null matches.')

    parser.add_argument('--capacity',
                    dest='capacity',
                    choices=['uniform', 'inf', 'top', '1to1'],
                    type=str,
                    default='uniform',
                    metavar='<str>',
                    help='Capacity method.')

    parser.add_argument('-c', '--cpus', '--jobs',
                    dest='cpus',
                    type=int,
                    default='1',
                    metavar='<int>',
                    help='Number of cores to use.')

    args = parser.parse_args()
    return args


# MAIN
if __name__ == '__main__':
    __spec__ = None
    args = interface()


    # Load loom and extract matrix
    sys.stdout.write('Reading data.\n')
    adata = sc.read(args.H5AD, sparse=False, X_name=args.layer)
    try:
        embed_mat = pd.DataFrame(adata.X.todense())
    except:
        embed_mat = pd.DataFrame(adata.X)
    embed_mat.columns = list(adata.var.index.values)
    embed_mat.index = list(adata.obs.index.values)

    sys.stdout.write('Splitting groups.\n')
    # Get column to split by
    split_vals = adata.obs.loc[:, args.split_col]
    # Get larger group
    group1_name = split_vals.value_counts().idxmax()
    group2_name = split_vals[split_vals!=group1_name].unique()[0]

    # Split groups
    group1 = embed_mat[split_vals==group1_name]
    group2 = embed_mat[split_vals!=group1_name]

    sys.stdout.write('Computing digraph.\n')
    # Calculate cost graph
    cost_graph = get_cost_knn_graph(
        source = group1,
        target = group2,
        knn_k = args.k,
        knn_n_jobs = args.cpus,
        null_cost_percentile = args.penalty,
        capacity_method = args.capacity
    )

    sys.stdout.write('Finding optimal matches.\n')
    # Find bipartite matches
    g1_idx, g2_idx = mcmf(cost_graph)
    # Remove NaNs
    any_nans = (np.isnan(g1_idx) | np.isnan(g2_idx))
    g1_idx = g1_idx[~any_nans].astype(int)
    g2_idx = g2_idx[~any_nans].astype(int)
    matches = pd.DataFrame({
        group1_name: group1.index[g1_idx],
        group2_name: group2.index[g2_idx]
    })
    matches.to_csv(args.out, sep='\t', index=False)
