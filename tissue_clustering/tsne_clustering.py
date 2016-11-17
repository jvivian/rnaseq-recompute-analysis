#!/usr/bin/env python2.7
"""
John Vivian
November, 2016
"""
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from itertools import combinations
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


def create_classification_vector(df_paths):
    vector, label = [], []
    for i, path in enumerate(df_paths):
        with open(path, 'r') as f:
            line = f.readline().strip().split('\t')
            vector.extend([x for x in line])
            label.extend([i for _ in line])
    return vector, np.array(label)


def create_combined_df(df_paths, pc_genes):
    dfs = [pd.read_csv(x, sep='\t', index_col=0) for x in df_paths]
    df = pd.concat(dfs, axis=1)
    df = remove_nonprotein_coding_genes(df, pc_genes)
    return df.T


def find_protein_coding_genes(gencode_path):
    pc_genes = set()
    with open(gencode_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[line.index('gene_type') + 1] == '"protein_coding";':
                    pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
    return list(pc_genes)


def remove_nonprotein_coding_genes(df, pc_genes):
    return df.ix[pc_genes]


def run_pca(df):
    pca = PCA(n_components=2)
    return pca.fit(df).transform(df)


def run_tsne(df):
    model = TSNE(n_components=2)
    return model.fit_transform(np.array(df))


def plot_dimensionality_reduction(ax, x, files, label, title, alpha=0.5):
    cm = plt.get_cmap('Accent')
    color_set = (cm(1. * i / len(files)) for i in xrange(len(files)))
    names = [os.path.basename(z).split('.tsv')[0] for z in files]
    length = [z for z in xrange(len(files))]
    for color, i, target_name in zip(color_set, length, names):
        ax.scatter(x[label == i, 0], x[label == i, 1], alpha=alpha, color=color, label=target_name)
    ax.legend(loc='best', fontsize=6)
    ax.set_title(title)


def determine_distance(x, label, files):
    smallest_distance = sys.maxint
    smallest_set = None
    for a, b in combinations(set(label), 2):
        dist = np.sqrt((x[label == a, 0].mean() - x[label == b, 0].mean()) ** 2 +
                       (x[label == a, 1].mean() - x[label == b, 1].mean()) ** 2)
        s1 = os.path.basename(files[a]).split('.tsv')[0]
        s2 = os.path.basename(files[b]).split('.tsv')[0]
        if 'tumor' in s1 or 'tumor' in s2:
            print 'Distance for: {} {}: {} '.format(s1, s2, dist)
            if dist < smallest_distance:
                smallest_set = (s1, s2)
                smallest_distance = dist
    print '\nSmallest Distance: ' + ', '.join(smallest_set) + ': ' + str(smallest_distance)


def remove_outliers(x, label, vector, files, title, percentile=98, alpha=0.5):
    # Create outlier mask
    outlier_mask = (x[:, 0] < np.percentile(x[:, 0], percentile)) & (x[:, 1] < np.percentile(x[:, 1], percentile))
    # Refit / transform with the mask applied
    x_mask = run_tsne(x[outlier_mask])
    # Apply mask to labels
    label_mask = label[outlier_mask]
    f, ax = plt.subplots()
    plot_dimensionality_reduction(ax, x_mask, files, label_mask, title=title, alpha=alpha)
    return x_mask, label_mask, set(vector) - set(np.array(vector)[outlier_mask])


def split_tcga_tumor_normal(directory):
    files = [os.path.join(directory, x) for x in os.listdir(directory)
             if not x.endswith('normal.tsv') and not x.endswith('tumor.tsv')]
    tcga_dfs = []
    for f in files:
        with open(f, 'r') as f_out:
            line = f_out.readline().strip().split('\t')
            if all('TCGA-' in x for x in line):
                tcga_dfs.append(f)
    # Create df splits
    for f in tcga_dfs:
        df = pd.read_csv(f, sep='\t', index_col=0)
        tumor = [x for x in df.columns if x.endswith('-01')]
        normal = [x for x in df.columns if x.endswith('-11')]
        if normal:
            df[normal].to_csv(f[:-4] + '_normal.tsv', sep='\t')
        df[tumor].to_csv(f[:-4] + '_tumor.tsv', sep='\t')
    return [os.path.join(directory, x) for x in os.listdir(directory) if os.path.join(directory, x) not in tcga_dfs]
