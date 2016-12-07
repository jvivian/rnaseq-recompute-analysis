import argparse
import logging
import os

from concurrent.futures import ThreadPoolExecutor

from deprecated_code.GTEX_matched import create_deseq2_inputs, run_deseq2, write_de_script
from utils import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main():
    """
    First run create_project.py

    H: Running DE for each tissue compared to all other tissues should produce tissue-specific genesets

    WARNING: Given each tissue is being run against all ~8k samples, this uses
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna-seq-directory', type=str, required=True,
                        help='Full path to rna-seq_analysis directory, created by running create_project.py')
    parser.add_argument('--cores', type=int, help='Number of cores to use when running R.', default=1)
    params = parser.parse_args()
    root_dir = params.rna_seq_directory
    experiment_dir = os.path.join(root_dir, 'experiments/gtex-matched-tissues')
    tissue_pairs = os.listdir(os.path.join(root_dir, 'data/tissue-pairs'))

    log.info('Creating subdirectories for experiment: ' + experiment_dir)
    [mkdir_p(os.path.join(experiment_dir, x, 'plots')) for x in tissue_pairs]
    [write_de_script(os.path.join(experiment_dir, x)) for x in tissue_pairs]

    # TODO: Apply GTEx RIN filter / outlier fitering once available
    log.info('Creating ordered dataframe of all GTEx tissues')
    pc_paths = [os.path.join(root_dir, 'data/tissue-pairs', x, 'combined-gtex-tcga-counts-protein-coding.tsv')
                for x in tissue_pairs]
    create_deseq2_inputs(pc_paths, experiment_dir)

    log.info('Starting DESeq2 run using {} cores'.format(params.cores))
    tissue_dirs = [os.path.join(experiment_dir, tissue) for tissue in tissue_pairs]
    with ThreadPoolExecutor(max_workers=params.cores) as executor:
        executor.map(run_deseq2, tissue_dirs)

    log.info('Adding mapped genes')


if __name__ == '__main__':
    main()
