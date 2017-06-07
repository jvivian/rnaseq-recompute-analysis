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
import pickle

from utils import mkdir_p, flatten

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def create_candidate_pairs(df, root_dir):
    pair_file = os.path.join(os.path.dirname(__file__), 'candidate_tissues.csv')
    combined_metadata = os.path.join(root_dir, 'metadata/tcga_gtex_metadata_intersect.tsv')
    metadata = pd.read_csv(combined_metadata, index_col=0, sep='\t')

    # For each candidate tissue (and possible pairing) create a dataframe for downstream clustering
    with open(pair_file, 'r') as f:
        for line in f:
            line = line.strip().split(',') if ',' in line else [line.strip()]

            # Pull samples that match our tissue
            samples = flatten([list(metadata[metadata.tissue == x].index) for x in line])

            # Filter out samples that aren't in our expression dataset
            samples = [x for x in samples if x in df.columns]
            if samples:
                name = os.path.join(root_dir, 'data/candidate-tissues', '-'.join(line))
                mkdir_p(name)
                log.info('Subsetting and saving dataframe: {}'.format(os.path.basename(name)))
                df[samples].T.to_csv(os.path.join(name, 'exp.tsv'), sep='\t')  # Samples by genes for clustering

                log.info('Clustering: {}'.format(os.path.basename(name)))
                output_path = os.path.join(name, 'cluster.html')
                cluster_df(df[samples].T, root_dir, output_path=output_path)


def cluster_df(df, root_dir, output_path, title='Bokeh Plot', norm=True):
    """df needs to be samples by features"""
    from bokeh.charts import Scatter, output_file
    from bokeh.palettes import Category20c
    from sklearn.manifold import TSNE
    from sklearn.decomposition import TruncatedSVD

    log.info('Dataframe shape: {}'.format(df.shape))
    df = df.apply(lambda x: np.log2(x + 1)) if norm else df
    log.info('\tReducing feature space to 50 with TruncatedSVD')
    y = TruncatedSVD(n_components=50, random_state=0).fit_transform(np.array(df))
    log.info('\tRunning t-SNE')
    model = TSNE(n_components=2, random_state=1, perplexity=50, learning_rate=1000, verbose=2)
    z = model.fit_transform(np.array(y))

    # Load pickled map obojects
    tissue_map = pickle.load(open(os.path.join(root_dir, 'data/objects/tissue_map.pickle'), 'rb'))
    type_map = pickle.load(open(os.path.join(root_dir, 'data/objects/type_map.pickle'), 'rb'))

    # Preparing Bokeh Plot
    samples = df.index
    pdf = pd.DataFrame()
    pdf['sample'] = samples
    pdf['tissue'] = [tissue_map[x].capitalize() for x in samples]
    pdf['x'] = z[:, 0]
    pdf['y'] = z[:, 1]
    pdf['type'] = [type_map[x] for x in samples]

    tooltips = [
        ('Tissue', '@tissue'),
        ('Type', '@type'),
        ('Sample', '@sample')]

    p = Scatter(pdf, x='x', y='y', title=title,
                xlabel="1", ylabel="2",
                color='tissue',
                tooltips=tooltips,
                legend=True,
                plot_width=1024, plot_height=1024,
                palette=Category20c[20],
                responsive=True)

    p.title.align = 'center'
    output_file(output_path)


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
    output_path = os.path.join(output_dir, name + '.tsv')
    if not os.path.exists(output_path):
        sub = df[df.index.isin(samples)]
        sub = prepare_for_de(sub)
        sub = sub.apply(lambda x: (2**x) - 1)  # Reverse of Xena normalization
        sub.to_csv(output_path, sep='\t')


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


def filter_nonprotein_coding_genes(df, gencode_path):
    pc_genes = set()
    with open(gencode_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[line.index('gene_type') + 1] == '"protein_coding";':
                    pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
    pc_genes = list(pc_genes)
    return df.loc[pc_genes]
