import errno
import os
import pickle

import pandas as pd


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def add_gene_names(df_path, gene_map_pickle):
    """Assumes ensemble gene names are the index"""
    gene_map = pickle.load(open(gene_map_pickle, 'rb'))
    df = pd.read_csv(df_path, sep='\t', index_col=0)

    gene_names = [gene_map[x] if x in gene_map.keys() else x for x in df.index]
    df['geneId'] = df.index
    df.index = gene_names

    return df
