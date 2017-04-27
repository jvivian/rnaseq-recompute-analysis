import logging
import os

import matplotlib

from experiments.tissue_clustering import TissueClustering

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class CBNTClustering(TissueClustering):

    def __init__(self, root_dir):
        super(CBNTClustering, self).__init__(root_dir)
        self.experiment_dir = os.path.join(root_dir, 'experiments/CBNT-Clustering')
        self.tsne_dir = os.path.join(self.experiment_dir, 'tsne')
        self.pca_dir = os.path.join(self.experiment_dir, 'pca')
        self.pickles = os.path.join(self.experiment_dir, 'pickles')
        self.expression_dataframes = self.cbnt_pc

    def setup(self):
        dirtree = [self.tsne_dir, self.pca_dir, self.tcga, self.tcga_matched] + \
                  [os.path.join(self.pickles, x) for x in ['pca', 'tsne']]
        self.create_directories(dirtree)

    def run_experiment(self):
        log.info('Running t-SNE Clustering')
        self.run_clustering(mode='tsne')

        log.info('Running PCA Clustering')
        self.run_clustering(mode='pca')
