import os
import logging

from experiments.tcga_tumor_vs_normal import TcgaTumorVsNormal
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class GTExVsTCGA(TcgaTumorVsNormal):

    def __init__(self, root_dir, cores):
        super(TcgaTumorVsNormal, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/gtex-vs-tcga')
        self.vector_dir = os.path.join(self.experiment_dir, 'vectors')
        self.results_dir = os.path.join(self.experiment_dir, 'results')
        self.plots_dir = os.path.join(self.experiment_dir, 'plots')
        self.script_path = None
        self.vectors = []

    def setup(self):
        dirtree = [self.vector_dir, self.results_dir] + [os.path.join(self.plots_dir, x) for x in self.tissues]
        self.create_directories(dirtree)

        self.script_path = write_script(self.deseq2_script, directory=self.experiment_dir)

        log.info('Writing out tissue vectors')
        for df in self.protein_coding_paths:
            tissue = os.path.basename(os.path.dirname(df))
            with open(df, 'r') as f_in:
                samples = [x.replace('-', '.') for x in f_in.readline().strip().split('\t')]
                tcga = [x for x in samples if x.startswith('TCGA') and x.endswith('01')]
                gtex = [x for x in samples if x.startswith('GTEX')]

                vector_path = os.path.join(self.vector_dir, tissue + '-vector')
                with open(vector_path, 'w') as f:
                    f.write('\n'.join(tcga + gtex))

                disease_vector = ['T'] * len(tcga) + ['N'] * len(gtex)
                disease_vector_path = os.path.join(self.vector_dir, tissue + '-disease')
                with open(disease_vector_path, 'w') as f:
                    f.write('\n'.join(disease_vector))
