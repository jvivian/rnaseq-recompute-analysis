import logging
import os
from abc import abstractmethod, ABCMeta

from utils import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class AbstractExperiment(object):

    __metaclass__ = ABCMeta

    def __init__(self, root_dir):
        super(AbstractExperiment, self).__init__()

        self.root_dir = root_dir
        self.gencode_path = os.path.join(root_dir, 'metadata/gencode.v23.annotation.gtf')
        self.tissue_dataframe_dir = os.path.join(root_dir, 'data/tissue-dataframes')
        self.tissue_pair_dir = os.path.join(root_dir, 'data/tissue-pairs')
        self.tissues = os.listdir(os.path.join(root_dir, 'data/tissue-pairs'))
        self.gene_map = os.path.join(root_dir, 'metadata/gene_map.pickle')

        self.protein_coding_paths = [
            os.path.join(self.tissue_pair_dir, x, 'combined-gtex-tcga-counts-protein-coding.tsv')
            for x in self.tissues]

    def create_directories(self, dirtree):
        """
        Iterates over dirtree to make directories if they do not already exist

        :param list[str] dirtree:
        """
        log.info('Creating directory trees')
        for path in dirtree:
            log.debug(path)
            mkdir_p(path)

    @abstractmethod
    def setup(self):
        raise NotImplementedError

    @abstractmethod
    def run_experiment(self):
        raise NotImplementedError

    @abstractmethod
    def teardown(self):
        raise NotImplementedError
