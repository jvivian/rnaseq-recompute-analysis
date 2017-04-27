import logging
import os
import pickle
from itertools import combinations

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import gmean
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from tqdm import tqdm

from experiments.AbstractExperiment import AbstractExperiment
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class TissueClustering(AbstractExperiment):

    def __init__(self, root_dir):
        super(TissueClustering, self).__init__(root_dir)
        self.experiment_dir = os.path.join(root_dir, 'experiments/tissue-clustering')
        self.tsne_dir = os.path.join(self.experiment_dir, 'tsne')
        self.pca_dir = os.path.join(self.experiment_dir, 'pca')
        self.pickles = os.path.join(self.experiment_dir, 'pickles')
        self.tcga = os.path.join(self.experiment_dir, 'tcga-only')
        self.tcga_matched = os.path.join(self.experiment_dir, 'tcga-matched')
        self.expression_dataframes = self.protein_coding_paths

    def setup(self):
        dirtree = [self.tsne_dir, self.pca_dir, self.tcga, self.tcga_matched] + \
                  [os.path.join(self.pickles, x) for x in ['pca', 'tsne', 'tcga-only', 'tcga-matched']]
        self.create_directories(dirtree)

    def run_experiment(self):
        log.info('Running t-SNE Clustering')
        self.run_clustering(mode='tsne')

        log.info('Running PCA Clustering')
        self.run_clustering(mode='pca')

        log.info('Running TCGA Clustering')
        self.run_clustering(mode='tcga-only')

        log.info('Running TCGA-Matched Clustering')
        self.run_clustering(mode='tcga-matched')

    def run_clustering(self, mode='tsne'):
        for tissue_path in tqdm(sorted(self.expression_dataframes)):
            tissue = os.path.basename(os.path.dirname(tissue_path))
            pickle_path = os.path.join(self.pickles, mode, tissue + '.pickle')
            if os.path.exists(pickle_path):
                log.info('Pickle file found, loading: ' + pickle_path)
                x, label = pickle.load(open(pickle_path, 'rb'))
            else:
                df = pd.read_csv(tissue_path, sep='\t', index_col=0)
                if mode == 'tcga-only':
                    tcga_cols = [x for x in df.columns if 'TCGA-' in x]
                    df = df[tcga_cols]
                elif mode == 'tcga-matched':
                    samples = list(df.columns)
                    barcodes = [x[:-3] for x in samples]
                    matched = list(set(x for x in barcodes if x + '-11' in samples and x + '-01' in samples))
                    if matched:
                        matched_vector = [f(x) for x in matched for f in (lambda f1: f1 + '-01', lambda f2: f2 + '-11')]
                        df = df[matched_vector]
                    else:
                        continue

                label = get_label(df)
                df = df.T  # Transpose so dataframe is samples x genes

                # log normalization
                # Also experimented with Size Factor Rescaling (from DESeq2) and Quantile Normalization
                # Log normalization seemed sufficient, is fast, and straight forward
                df = df.apply(lambda y: np.log(y + 1))

                # Cluster by method
                if mode == 'pca':
                    x = run_pca(df)
                else:
                    x = run_tsne(df)

                # Save relavant info in pickle
                info = (x, label)
                with open(pickle_path, 'w') as f:
                    pickle.dump(info, f)

            # Plot
            f, ax = plt.subplots()
            plot_dimensionality_reduction(ax, x, label, title=tissue)
            f.savefig(os.path.join(self.experiment_dir, mode, tissue + '.png'), format='png', dpi=300)
            plt.close()

    def teardown(self):
        pass


def get_label(df):
    labels = []
    for sample in df.columns:
        if 'GTEX-' in sample:
            labels.append(0)
        elif sample.endswith('01'):
            labels.append(1)
        elif sample.endswith('11'):
            labels.append(2)
        else:
            labels.append(3)
    return np.array(labels)


def size_factor_scale(df):
    # Calculate denominator
    # The geometric mean for all genes across all samples
    # We'll pad the counts with a +1 so the geometric mean is never 0
    df = df.apply(lambda x: x + 1)
    denom = gmean(df, axis=1)  # Length n
    temp = df.divide(denom, axis=0)
    size_factors = temp.median(axis=0)
    return df.divide(size_factors, axis=1).apply(lambda x: np.log2(x))


def plot_dimensionality_reduction(ax, x, label, title, alpha=0.5):
    names = ['GTEX', 'TCGA-T', 'TCGA-N', 'TCGA-O']
    length = [0, 1, 2, 3]
    cm = plt.get_cmap('Accent')
    color_set = (cm(1. * i / len(names)) for i in xrange(len(names)))
    for color, i, target_name in zip(color_set, length, names):
        ax.scatter(x[label == i, 0], x[label == i, 1], alpha=alpha, color=color, label=target_name)
    ax.legend(loc='best', fontsize=6)
    ax.set_title(title)
    log.info('Plotting: ' + title)

    # Cluster distance measurements
    tcga_tn_dist = None
    tcga_gtex_dist = None
    for a, b in combinations(set(label), 2):
        dist = np.sqrt((x[label == a, 0].mean() - x[label == b, 0].mean()) ** 2 +
                       (x[label == a, 1].mean() - x[label == b, 1].mean()) ** 2)
        log.info('Distance for: {:>8} {:>8} {:>8}'.format(names[a], names[b], dist))
        if sum([a, b]) == 2:  # GTEX / TCGA-N
            tcga_gtex_dist = dist
        elif sum([a, b]) == 3:  # TCGA-T / TCGA-N
            tcga_tn_dist = dist
    if tcga_gtex_dist < tcga_tn_dist:
        log.info('GTEx and TCGA Normal are more similar!')
    else:
        log.info('TCGA tumor / normal are more similar')


def run_pca(df):
    pca = PCA(n_components=2, random_state=1)
    return pca.fit(df).transform(df)


def run_tsne(df):
    model = TSNE(n_components=2)
    return model.fit_transform(np.array(df))
