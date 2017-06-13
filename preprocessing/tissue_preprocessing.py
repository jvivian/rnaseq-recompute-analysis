#!/usr/bin/env python2.7
"""
John Vivian
September, 2016
"""
import logging
import os
import pickle

import numpy as np
import pandas as pd
from utils import mkdir_p, flatten

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def create_tissue_pairs(df, root_dir, tsv_name='tcga-gtex-exp.tsv'):
    pair_file = os.path.join(os.path.dirname(__file__), 'candidate_tissues.csv')
    # For each candidate tissue create a dataframe
    tissues = []
    with open(pair_file, 'r') as f:
        for line in f:
            line = line.strip().split(',') if ',' in line else [line.strip()]
            tissues.append(line)
            samples = get_samples_for_tissue(df, root_dir, line)

            if samples:
                tissue = '-'.join(line)
                tissue_path = os.path.join(root_dir, 'data/tissue-pairs', tissue)
                mkdir_p(tissue_path)
                tsv_path = os.path.join(tissue_path, tsv_name)
                if not os.path.exists(tsv_path):
                    log.info('Subsetting and saving dataframe: {}'.format(os.path.basename(tissue_path)))
                    df[samples].to_csv(tsv_path, sep='\t')
                    tissues.append(tissue)
    return tissues


def cluster_df(df, root_dir, output_path, title='Bokeh Plot', norm=True, colorby='type'):
    """df needs to be samples by features"""
    from bokeh.charts import Scatter, save
    from bokeh.palettes import Category10
    from bokeh.embed import autoload_static
    from bokeh.resources import CDN
    from sklearn.manifold import TSNE
    from sklearn.decomposition import TruncatedSVD

    log.info('Dataframe shape: {}'.format(df.shape))
    log.info('Applying log2(x + 1) normalization to raw counts')
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
    pdf['tissue'] = [tissue_map[x].capitalize() for x in samples if x in tissue_map]
    pdf['x'] = z[:, 0]
    pdf['y'] = z[:, 1]
    pdf['type'] = [type_map[x] for x in samples if x in type_map]
    dataset = []
    for s in samples:
        if s.startswith('TCGA'):
            if s.endswith('11'):
                dataset.append('TCGA-Normal')
            else:
                dataset.append('TCGA-Tumor')
        else:
            dataset.append('GTEX')
    pdf['dataset'] = dataset

    # Determine number of colors for palette
    if 3 < len(pdf.type.unique()) < 10:
        num_colors = len(pdf.type.unique())
    elif len(pdf.type.unique()) >= 10:
        num_colors = 10
    else:
        num_colors = 3

    tooltips = [
        ('Tissue', '@tissue'),
        ('Type', '@type'),
        ('Sample', '@sample')]

    log.info('Creating Bokeh Scatterplot')
    p = Scatter(pdf, x='x', y='y', title=title,
                xlabel="1", ylabel="2",
                color=colorby,
                tooltips=tooltips,
                legend=True,
                plot_width=1024, plot_height=1024,
                palette=Category10[num_colors],
                responsive=True)

    log.info('Outputting HTML to: {}'.format(output_path))
    p.title.align = 'center'
    save(p, output_path, title=title)

    # Save Javascript and HTMl tag versions of plots
    output_dir = os.path.join(os.path.dirname(output_path), 'javascript')
    mkdir_p(output_dir)
    name = os.path.basename(output_path).split('.')[0] + '.js'
    output_name = os.path.join(output_dir, name)
    js, tag = autoload_static(p, CDN, "js/bokeh/{}".format(name))
    with open(output_name, 'w') as f:
        f.write(js)
    with open(os.path.join(output_dir, "tags.txt"), 'a') as f:
        f.write(tag)


def filter_nonprotein_coding_genes(df, root_dir):
    gencode_path = os.path.join(root_dir, 'metadata/gencode.v23.annotation.gtf')
    pc_genes = set()
    with open(gencode_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[line.index('gene_type') + 1] == '"protein_coding";':
                    pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
    pc_genes = list(pc_genes)
    return df.loc[pc_genes]


def filter_samples_by_metadata(df, root_dir):
    """Removes samples not in metadata"""
    metadata = pd.read_csv(os.path.join(root_dir, 'metadata/tcga_gtex_metadata_intersect.tsv'), sep='\t', index_col=0)
    met_samples = list(metadata.index)
    samples = [x for x in df.columns if x in met_samples]
    log.info('Reducing dataframe from {} to {} samples'.format(len(df.columns), len(samples)))
    return df[samples], samples


def get_samples_for_tissue(df, root_dir, samples):
    combined_metadata = os.path.join(root_dir, 'metadata/tcga_gtex_metadata_intersect.tsv')
    metadata = pd.read_csv(combined_metadata, index_col=0, sep='\t')
    # Pull samples that match our tissue
    samples = flatten([list(metadata[metadata.tissue == x].index) for x in samples])
    # Filter out samples that aren't in our expression dataset
    return [x for x in samples if x in df.columns]


def cluster_tissues(df, root_dir, tissues, title):
    for tissue in tissues:
        tissue_name = '-'.join(tissue)
        output_dir = os.path.join(root_dir, 'data/clustering', tissue_name)
        mkdir_p(output_dir)
        output_path = os.path.join(output_dir, '{}.html'.format(title))
        if not os.path.exists(output_path):
            log.info('Clustering: {}'.format(tissue))
            # Get samples that correspond to our tissue
            tissue_samples = get_samples_for_tissue(df, root_dir, samples=tissue)
            # Note the transpose of the matrix is passed to get: samples by genes
            cluster_df(df[tissue_samples].T, root_dir, output_path=output_path, title=tissue_name)


def cluster_entire_dataset(df, root_dir, base_title):
    output_dir = os.path.join(root_dir, 'data/clustering', 'all')
    mkdir_p(output_dir)
    for cluster_type in ['tissue', 'type', 'dataset']:
        output_path = os.path.join(output_dir, '{}-{}.html'.format(base_title, cluster_type))
        if not os.path.exists(output_path):
            log.info('Clustering entire dataset by: {}'.format(cluster_type))
            cluster_df(df.T, root_dir, output_path=output_path,
                       title='t-SNE Clustering of TCGA and GTEx by {}'.format(cluster_type), colorby='tissue')
