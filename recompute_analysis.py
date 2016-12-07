# coding: utf-8
import argparse
import logging

from experiments.gtex_one_vs_all import GtexOneVsAll
from experiments.clustering import Clustering

from utils import title_gtex_one_vs_all
from utils import title_clustering

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def runner(instance):
    instance.setup()
    instance.run_experiment()
    instance.teardown()


def main():
    """
    Launchpoint for all experiments associated with the CGL RNA-seq recompute analysis
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Parent parser arguments
    parser.add_argument('--project-dir', help='Full path to the project directory (rna-seq-analysis) created'
                                              'by launching the create_project.py or using the create-project'
                                              'entrypoint.')

    # GTEx One Versus All
    parser_gtex_ova = subparsers.add_parser('gtex-one-vs-all', parents=[parser], help='Run GTEx One Vs All comparison.')
    parser_gtex_ova.add_argument('--cores', help='The number of cores to utilize during run. '
                                                 'Warning: This experiment is memory intensive (~35G per core).')

    # GTEx Pairwise
    parser_gtex_pairwise = subparsers.add_parser('gtex-pairwise', parents=[parser], help='Run GTEx Pairwise Comparison')
    parser_gtex_pairwise.add_argument('--cores', required=True, help='Number of cores to utilize during run.')

    # Clustering
    subparsers.add_parser('clustering', parents=[parser], help='Run PCA / t-SNE clustering.')

    # Execution
    params = parser.parse_args()

    if params.command == 'gtex-one-vs-all':
        log.info(title_gtex_one_vs_all())
        runner(GtexOneVsAll(params.project_dir, params.cores))

    elif params.command == 'clustering':
        log.info(title_clustering())
        runner(Clustering(params.project_dir))

if __name__ == '__main__':
    main()
