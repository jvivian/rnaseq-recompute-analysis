# coding: utf-8
import argparse
import logging

from experiments.gtex_one_vs_all import GtexOneVsAll
from experiments.clustering import Clustering
from experiments.gtex_pairwise import GtexPairwise

from utils import title_gtex_one_vs_all, title_gtex_pairwise
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
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # GTEx One Versus All
    parser_gtex_ova = subparsers.add_parser('gtex-one-vs-all', parents=[parser], help='Run GTEx One Vs All comparison.')
    parser_gtex_ova.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')
    parser_gtex_ova.add_argument('--cores', required=True, help='The number of cores to utilize during run. Warning: '
                                                                'This experiment is memory intensive (~35G per core).')

    # GTEx Pairwise
    parser_gtex_pairwise = subparsers.add_parser('gtex-pairwise', help='Run GTEx Pairwise Comparison')
    parser_gtex_pairwise.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')
    parser_gtex_pairwise.add_argument('--cores', required=True, help='Number of cores to utilize during run.')

    # Clustering
    parser_clustering = subparsers.add_parser('clustering', parents=[parser], help='Run PCA / t-SNE clustering.')
    parser_clustering.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')

    # Execution
    params = parser.parse_args()

    if params.command == 'gtex-one-vs-all':
        log.info(title_gtex_one_vs_all())
        runner(GtexOneVsAll(params.project_dir, params.cores))

    elif params.command == 'clustering':
        log.info(title_clustering())
        runner(Clustering(params.project_dir))

    elif params.command == 'gtex-pairwise':
        log.info(title_gtex_pairwise())
        runner(GtexPairwise(params.project_dir, params.cores))

if __name__ == '__main__':
    main()
