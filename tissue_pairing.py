#!/usr/bin/env python2.7
"""
John Vivian
September, 2016
"""
import os

import pandas as pd


def concat_frames(gtex_df, tcga_df, name, path='/mnt/tissues/'):
    """
    Concantenate dataframes produces by tissue_preprocessing
    Creates two files:
        1. num_samples which contains the number of gtex and tcga samples
        2. A concantenated dataframe (axis=1) of gtex_df and tcga_df

    :param pd.DataFrame gtex_df: GTex dataframe
    :param pd.DataFrame tcga_df: TCGA dataframe
    :param str name: Name of tissue
    :param str path: Path to where directory and tissue will be created
    """
    output_dir = os.path.join(path, name)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(path + name + '/num_samples', 'w') as f:
        f.write('{}\t{}\n'.format(len(gtex_df.columns), len(tcga_df.columns)))
    df = pd.concat([gtex_df, tcga_df], axis=1)
    combined_name = os.path.join(path, name, name + '_combined.tsv')
    df.to_csv(combined_name, sep='\t')


# Liver
gl = pd.read_csv('/mnt/tissues/liver/Liver.tsv', sep='\t', index_col=0)
tl = pd.read_csv('/mnt/tissues/liver/Liver_hepatocellular_carcinoma.tsv', sep='\t', index_col=0)
concat_frames(gl, tl, name='liver')

# Thyroid
gt = pd.read_csv('/mnt/tissues/thyroid/Thyroid.tsv', sep='\t', index_col=0)
tt = pd.read_csv('/mnt/tissues/thyroid/Thyroid_carcinoma.tsv', sep='\t', index_col=0)
concat_frames(gt, tt, name='thyroid')

# Prostate
gp = pd.read_csv('/mnt/tissues/prostate/Prostate.tsv', sep='\t', index_col=0)
tp = pd.read_csv('/mnt/tissues/prostate/Prostate_adenocarcinoma.tsv', sep='\t', index_col=0)
concat_frames(gp, tp, name='prostate')

# Stomach
gs = pd.read_csv('/mnt/tissues/stomach/Stomach.tsv', sep='\t', index_col=0)
ts = pd.read_csv('/mnt/tissues/stomach/Stomach_adenocarcinoma.tsv', sep='\t', index_col=0)
concat_frames(gs, ts, name='stomach')
