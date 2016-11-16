#!/usr/bin/env python2.7
"""
John Vivian
September, 2016
"""
import logging
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


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
    :param str output_dir: Path to output directory for tables
    """
    # GTEx metadata
    gtex_meta = pd.read_csv(gtex_metadata, delimiter='\t')
    gtex_meta.replace('<not provided>', np.nan, inplace=True)
    gtex_meta.columns = [x[:-2] for x in gtex_meta.columns]

    # TCGA metadata
    # Fix naming convention - Xena's barcode consists of only 15 characters
    tcga_meta = pd.read_csv(tcga_metadata, delimiter='\t')
    tcga_meta.barcode = [x[:15] for x in tcga_meta.barcode]

    log.info('Reading in in GTEx expression table')
    gt = pd.read_csv(gtex_expression, delimiter='\t')
    gt = process_raw_xena_df(gt)

    log.info('Reading in TCGA expression table')
    tc = pd.read_csv(tcga_expression, delimiter='\t')
    tc = process_raw_xena_df(tc)

    log.info('Creating GTEx tissue dataframes')
    for tissue in tqdm(gtex_meta.body_site.unique()):
        subtype = gtex_meta[gtex_meta.body_site == tissue]
        name = '_'.join(' '.join(tissue.split('-')).split())
        create_subframe(gt, samples=subtype.Sample_Name, name=name, output_dir=output_dir)

    log.info('Creating TCGA tissue dataframes')
    for tissue in tqdm(tcga_meta.disease_name.unique()):
        subtype = tcga_meta[tcga_meta.disease_name == tissue]
        name = '_'.join(' '.join(tissue.split('-')).split())
        create_subframe(tc, samples=subtype.barcode, name=name, output_dir=output_dir)


def concat_frames(gtex_df_paths, tcga_df_path, output_path):
    """
    Concantenate tissue dataframes to create a single combined data frame

    :param list gtex_df_paths: Path(s) to GTEx dataframe(s)
    :param str tcga_df_path: Path to TCGA dataframe
    :param str output_path: Path to where directory and tissue will be created
    """
    log.debug('Combining: {}\t{}'.format(gtex_df_paths, tcga_df_path))
    gtex = [pd.read_csv(gtex_df, sep='\t', index_col=0) for gtex_df in gtex_df_paths]
    tcga = pd.read_csv(tcga_df_path, sep='\t', index_col=0)
    df = pd.concat(gtex + [tcga], axis=1)
    df.to_csv(output_path, sep='\t')


# TODO: Replace with GTFParser
def remove_nonprotein_coding_genes(df_path, gencode_path):
    """
    Removes non-protein coding genes which can skew normalization

    :param str df_path: Path to combined-gtex-tcga-counts.tsv dataframe
    :param str gencode_path: Path to gencode GTF
    """
    df = pd.read_csv(df_path, sep='\t', index_col=0)
    pc_genes = set()
    with open(gencode_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[line.index('gene_type') + 1] == '"protein_coding";':
                    pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
    # Subset list of protein-coding genes
    pc_genes = list(pc_genes)
    df = df.ix[pc_genes]
    output_path = os.path.join(os.path.dirname(df_path), 'combined-gtex-tcga-counts-protein-coding.tsv')
    df.to_csv(output_path, sep='\t')
