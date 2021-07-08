#!/usr/bin/env python3

import os
import sys
import h5py
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.mixture import GaussianMixture


# FUNC
def interface():
    parser = argparse.ArgumentParser(description='Extracts UMIs supported by enough reads using a GMM.')

    parser.add_argument('CELLRANGER_DIR',
                    type=str,
                    metavar='<dir>',
                    help='Directory with cellranger output.')

    parser.add_argument('-o', '--out-dir',
                    dest='out',
                    type=str,
                    default='./',
                    metavar='<f>',
                    help='Output directory.')

    parser.add_argument('--guide_id',
                    dest='guide_id',
                    type=str,
                    default='_gene',
                    metavar='<str>',
                    help='Identifier to distinguish guide mRNAs from other genes.')

    parser.add_argument('-p', '--cutoff',
                    dest='cutoff',
                    type=float,
                    default='0.5',
                    metavar='<float>',
                    help='Probability cutoff for GMM.')

    parser.add_argument('-n', '--name',
                    dest='name',
                    type=str,
                    metavar='<str>',
                    help='Sample name to add to barcode.')

    args = parser.parse_args()
    return args


def get_guide_umis(mol_info_file, guide_id='_gene', min_count=1):
    # Get arrays from H5
    with h5py.File(mol_info_file, 'r') as f:
        read_counts = np.array(f['count'])
        bc_idx = np.array(f['barcode_idx'])
        barcodes = np.array(f['barcodes'])
        read_bc = np.array(barcodes[bc_idx], dtype='str')
        feat_idx = np.array(f['feature_idx'])
        features = np.array(f['features']['name'])
        read_feature = np.array(features[feat_idx], dtype='str')

    # Make dataframe
    reads_df = pd.DataFrame({
        'count': read_counts,
        'cell': read_bc,
        'feature': read_feature
    })

    # Initial filter for read count > min_count
    reads_df = reads_df[reads_df['count'] > min_count]
    # Filter only guide mRNAs
    reads_df = reads_df[reads_df['feature'].str.contains(guide_id)]
    # Remove guide ID
    reads_df['feature'] = reads_df['feature'].str.replace(guide_id, '')
    # Log-transform counts
    reads_df['log_count'] = np.log10(reads_df['count'])
    return reads_df


def filter_umis_gmm(guide_umis, cutoff=1e-10):
    # Init GMM
    gmm = GaussianMixture(n_components=2)
    log_counts = np.array(guide_umis['log_count'])
    # Fit and predict class probs
    gmm_fit = gmm.fit(log_counts.reshape(-1, 1))
    gmm_prob = gmm.predict_proba(log_counts.reshape(-1, 1))
    gmm_class = gmm.predict(log_counts.reshape(-1, 1))
    # Find high class
    pos_mean = np.mean(log_counts[gmm_class==0])
    neg_mean = np.mean(log_counts[gmm_class==1])
    high_group = np.argmax([pos_mean, neg_mean])
    # Assign to df
    guide_umis['is_high'] = gmm_prob[:, high_group] >= cutoff
    return guide_umis


def plot_hist(data, file=None, vline=None):
    f, ax = plt.subplots()
    hist = sns.distplot(data, kde=False, bins=100, color='black')
    if vline:
        plt.axvline(vline, c='black')
    hist.set(xlabel='Log10 read counts', ylabel='# UMIs')
    sns.despine(offset=5, trim=True)
    if file:
        hist_fig = hist.get_figure()
        hist_fig.savefig(file)
    else:
        plt.show()


# MAIN
if __name__ == '__main__':
    __spec__ = None
    args = interface()

    CELLRANGER_DIR = os.path.abspath(args.CELLRANGER_DIR)
    out_dir = os.path.abspath(args.out)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    min_count = 0 if args.cutoff == 0 else 1

    # Get Molecule info and fit GMM
    mol_info_path = os.path.join(CELLRANGER_DIR, 'outs/molecule_info.h5')
    guide_umis = get_guide_umis(
        mol_info_path,
        guide_id = args.guide_id,
        min_count = min_count
    )
    if args.name:
        guide_umis['cell'] = args.name + '_' + guide_umis['cell'].astype(str)

    if args.cutoff > 0:
        guide_umis = filter_umis_gmm(guide_umis, cutoff=args.cutoff)
        boundary = guide_umis['log_count'][guide_umis['is_high']].min()
    else:
        guide_umis['is_high'] = True
        boundary = None

    # Plot histogram
    plot_out = os.path.join(out_dir, 'read_counts.pdf')
    plot_hist(guide_umis['log_count'], file=plot_out, vline=boundary)

    # Write guide assignments
    out_path = os.path.join(out_dir, 'umi_counts.tsv')
    guide_umis.to_csv(out_path, sep='\t', index = False)

    guide_counts = (guide_umis[guide_umis['is_high']]
        .groupby(['cell', 'feature'])
        .size()
        .reset_index(name='count')
    )
    out_path = os.path.join(out_dir, 'guide_counts.tsv')
    guide_counts.to_csv(out_path, sep='\t', index = False)

    guide_counts_mat = (guide_counts
        .pivot(index='cell', columns='feature', values='count')
        .reset_index()
        .fillna(0)
    )
    out_path = os.path.join(out_dir, 'guide_count_mat.tsv')
    guide_counts_mat.to_csv(out_path, sep='\t', index=False)
