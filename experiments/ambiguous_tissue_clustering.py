import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pickle
import shutil
import numpy as np

from tqdm import tqdm

from tissue_clustering.tsne_clustering import create_classification_vector
from tissue_clustering.tsne_clustering import create_combined_df
from tissue_clustering.tsne_clustering import find_protein_coding_genes
from tissue_clustering.tsne_clustering import plot_dimensionality_reduction
from tissue_clustering.tsne_clustering import run_pca
from tissue_clustering.tsne_clustering import run_tsne
from tissue_clustering.tsne_clustering import split_tcga_tumor_normal
from tissue_clustering.tsne_clustering import tissues
from experiments.abstractExperiment import AbstractExperiment
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
class AmbiguousTissueClustering(AbstractExperiment):

    def __init__(self, root_dir):
        super(AmbiguousTissueClustering, self).__init__(root_dir)
        self.root_dir = root_dir
        self.experiment_dir = os.path.join(root_dir, 'experiments/ambiguous-tissue-clustering')
        self.tsne_dir = os.path.join(self.experiment_dir, 'tsne')
        self.pca_dir = os.path.join(self.experiment_dir, 'pca')

    def setup(self):
        dirtree = [os.path.join(self.experiment_dir, x) for x in tissues] + [self.tsne_dir + self.pca_dir]
        self.create_directories(dirtree)

        log.info('Copying over requisite dataframes')
        tissue_dataframes = os.path.join(self.root_dir, 'data/tissue-dataframes')
        for tissue in tqdm(tissues):
            for tsv in tissues[tissue]:
                if not os.path.exists(os.path.join(self.experiment_dir, tissue, tsv)):
                    shutil.copy(os.path.join(tissue_dataframes, tsv), os.path.join(self.experiment_dir, tissue, tsv))

    def run_experiment(self):
        self.run_clustering(self.tsne_dir)
        self.run_clustering(self.pca_dir, name='pca')

    def teardown(self):
        pass

    def run_clustering(self, cluster_dir, name='tsne'):
        tsne_output = os.path.join(cluster_dir, name + '-cluster-data.pickle')
        if os.path.exists(tsne_output):
            log.info('Pickle file found, loading: ' + tsne_output)
            tsne = pickle.load(open(tsne_output, 'rb'))
        else:
            pc_genes = find_protein_coding_genes(self.gencode_path)
            tsne = {}
            for tissue in tqdm(tissues):
                files = split_tcga_tumor_normal(os.path.join(self.experiment_dir, tissue))
                vector, label = create_classification_vector(files)
                df = create_combined_df(files, pc_genes)

                # Log-normalizing expression data for PCA / t-SNE
                ln_df = df.apply(lambda y: np.log(y + 1))
                if name == 'tsne':
                    x = run_tsne(ln_df)
                else:
                    x = run_pca(ln_df)
                tsne[tissue] = (x, files, label)
            with open(tsne_output, 'w') as f:
                pickle.dump(tsne, f)

        log.info('Creating one large subplot')
        f, axes = plt.subplots(len(tissues), figsize=(8, 72))
        cm = plt.get_cmap('Accent')
        for i, tissue in enumerate(sorted(tsne)):
            matrix, files, label = tsne[tissue]
            color_set = (cm(1. * z / len(files)) for z in xrange(len(files)))
            names = [os.path.basename(x).split('.tsv')[0] for x in files]
            length = [x for x in xrange(len(files))]
            for color, l, target_name in zip(color_set, length, names):
                axes[i].scatter(matrix[label == l, 0], matrix[label == l, 1], alpha=0.5, color=color, label=target_name)
            axes[i].legend(loc='best', fontsize=6)
            axes[i].set_title(tissue)
        plt.savefig(os.path.join(cluster_dir, name + '-plots.pdf'), format='pdf')

        log.info('Creating per-tissue plots')
        for tissue in tsne:
            f, ax = plt.subplots()
            x, files, label = tsne[tissue]
            plot_dimensionality_reduction(ax, x, files, label, title=tissue)
            f.savefig(os.path.join(cluster_dir, tissue + '-' + name + '.png'), format='png', dpi=300)
            plt.close()
