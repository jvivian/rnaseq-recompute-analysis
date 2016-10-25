#!/usr/bin/env python2.7
"""
John Vivian
September, 2016
"""
import argparse
import os

import numpy as np
import pandas as pd
from tqdm import tqdm


def process_raw_xena_df(df):
    """
    Returns a dataframe in the format: samples x genes/isoforms

    :param pd.DataFrame df: Raw expression dataframe pulled from Xena
    :return: pd.DataFrame
    """
    df = df.T
    df.columns = df.iloc[0]
    print 'Samples: {}\tGenes: {}'.format(len(df.index), len(df.columns))
    return df.drop('sample', axis=0)


def prepare_for_de(df):
    """
    Returns gene/isoforms by samples

    :param pd.DataFrame df: Dataframe to be transposed
    :return: pd.DataFrame
    """
    df = df.T
    df.index.name = None
    return df


def create_subframe(df, samples, name, output_dir):
    """
    Creates a subframe from a dataframe given a set of samples
    Applies reverse normalization to get expected counts.

    :param pd.DataFrame df: Expression dataframe
    :param pd.Series samples: Selection of samples to subselect from
    :param str name: Name of output dataframe
    :param str output_dir: Path to output dataframe
    """
    sub = df[df.index.isin(samples)]
    sub = prepare_for_de(sub)
    sub = sub.apply(lambda x: (2**x) - 1)  # Reverse of Xena normalization
    sub.to_csv(os.path.join(output_dir, name + '.tsv'), sep='\t')


def create_subframes(gtex_metadata, tcga_metadata, gtex_expression, tcga_expression, output_dir):
    """
    Create subframes for every tissue

    :param str gtex_metadata: Path to GTEx metadata
    :param str tcga_metadata: Path to TCGA metadata
    :param str gtex_expression: Path to GTEx expression table
    :param str tcga_expression: Path to TCGA expression table
    """
    # GTEx metadata
    gtex_meta = pd.read_csv(gtex_metadata, delimiter='\t')
    gtex_meta.replace('<not provided>', np.nan, inplace=True)
    gtex_meta.columns = [x[:-2] for x in gtex_meta.columns]

    # TCGA metadata
    # Fix naming convention - Xena's barcode consists of only 15 characters
    tcga_meta = pd.read_csv(tcga_metadata, delimiter='\t')
    tcga_meta.barcode = [x[:15] for x in tcga_meta.barcode]

    # Read in expression dataframes
    gt = pd.read_csv(gtex_expression, delimiter='\t')
    gt = process_raw_xena_df(gt)

    tc = pd.read_csv(tcga_expression, delimiter='\t')
    tc = process_raw_xena_df(tc)

    for tissue in tqdm(gtex_meta.body_site.unique()):
        subtype = gtex_meta[gtex_meta.body_site == tissue]
        name = '_'.join(' '.join(tissue.split('-')).split())
        create_subframe(gt, samples=subtype.Sample_Name, name=name, output_dir=output_dir)

    for tissue in tqdm(tcga_meta.disease_name.unique()):
        subtype = tcga_meta[tcga_meta.disease_name == tissue]
        name = '_'.join(' '.join(tissue.split('-')).split())
        create_subframe(tc, samples=subtype.barcode, name=name, output_dir=output_dir)


def concat_frames(gtex_df, tcga_df, output_dir):
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
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    df = pd.concat([gtex_df, tcga_df], axis=1)
    combined_name = os.path.join(output_dir, 'combined-gtex-tcga.tsv')
    df.to_csv(combined_name, sep='\t')
