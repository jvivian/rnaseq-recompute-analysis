import argparse
import logging
import os

from concurrent.futures import ThreadPoolExecutor

from experiments.matched_de import create_matched_dataframe, write_de_script, run_deseq2
from utility_functions import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main():
    """
    First run create_project.py

    Runs DESeq2 for all TCGA samples with a tumor-normal pair
    :return:
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna-seq-directory', type=str, required=True,
                        help='Full path to rna-seq_analysis directory, created by running create_project.py')
    parser.add_argument('--cores', type=int, help='Number of cores to use when running R.', default=1)
    params = parser.parse_args()
    # Create subdirectories
    root_dir = params.rna_seq_directory
    experiment_dir = os.path.join(root_dir, 'experiments/matched-TCGA-tumor-normal')
    tissue_pairs = os.listdir(os.path.join(root_dir, 'data/tissue-pairs'))
    [mkdir_p(os.path.join(experiment_dir, x, 'plots')) for x in tissue_pairs]
    # Create matched dataframes
    for tissue in tissue_pairs:
        input_df = os.path.join(root_dir, 'data/tissue-pairs', tissue, 'combined-gtex-tcga-counts-protein-coding.tsv')
        output_df = os.path.join(experiment_dir, tissue, 'matched-tcga-counts.tsv')
        paired = create_matched_dataframe(input_df, output_df)
        # Write out DESeq2 R script
        if paired:
            write_de_script(os.path.join(experiment_dir, tissue))
        else:
            tissue_pairs.remove(tissue)
    # Run DESeq2 in batches
    log.info('Starting DESeq2 run using {} number of cores'.format(params.cores))
    with ThreadPoolExecutor(max_workers=params.cores) as executor:
        executor.map(run_deseq2, [os.path.join(experiment_dir, tissue) for tissue in tissue_pairs])


if __name__ == '__main__':
    main()
