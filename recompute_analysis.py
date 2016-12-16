# coding: utf-8
import argparse
import logging
import sys

from experiments.clustering import Clustering
from experiments.gtex_one_vs_all import GtexOneVsAll
from experiments.gtex_pairwise import GtexPairwise
from experiments.pairwise_gtex_vs_tcga import PairwiseTcgaVsGtex
from experiments.tcga_matched import TcgaMatched
from utils import title_clustering, cls
from utils import title_gtex_one_vs_all, title_gtex_pairwise, title_tcga_matched, title_pairwise_gtex_tcga

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def runner(instance):
    instance.setup()
    instance.run_experiment()
    instance.teardown()


def main():
    """
    Launchpoint for all experiments associated with the CGL RNA-seq recompute analysis

    First run create_project.py (or create-project entrypoint) to create the necessary
    prerequisite project directory
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # GTEx One Versus All
    parser_gtex_ova = subparsers.add_parser('gtex-one-vs-all', parents=[parser], help='Run GTEx One Vs All comparison.')
    parser_gtex_ova.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')
    parser_gtex_ova.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # GTEx Pairwise
    parser_gtex_pairwise = subparsers.add_parser('gtex-pairwise', help='Run GTEx Pairwise Comparison')
    parser_gtex_pairwise.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')
    parser_gtex_pairwise.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # Clustering
    parser_clustering = subparsers.add_parser('clustering', parents=[parser], help='Run PCA / t-SNE clustering.')
    parser_clustering.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis')

    # TCGA Matched
    parser_tcga_matched = subparsers.add_parser('tcga-matched', help='Run TCGA matched T/N analysis.')
    parser_tcga_matched.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    parser_tcga_matched.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # Paired TCGA v GTEx
    parser_pairwise = subparsers.add_parser('pairwise-gtex-tcga',
                                            help='Performs pairwise comparison between GTEx and TCGA')
    parser_pairwise.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    parser_pairwise.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        cls()
        parser.print_help()
        sys.exit(1)

    # Execution
    params = parser.parse_args()

    print params.command

    if params.command == 'gtex-one-vs-all':
        log.info(title_gtex_one_vs_all())
        runner(GtexOneVsAll(params.project_dir, params.cores))

    elif params.command == 'clustering':
        log.info(title_clustering())
        runner(Clustering(params.project_dir))

    elif params.command == 'gtex-pairwise':
        log.info(title_gtex_pairwise())
        runner(GtexPairwise(params.project_dir, params.cores))

    elif params.command == 'tcga-matched':
        log.info(title_tcga_matched())
        runner(TcgaMatched(params.project_dir, params.cores))

    elif params.command == 'pairwise-gtex-tcga':
        log.info(title_pairwise_gtex_tcga())
        runner(PairwiseTcgaVsGtex(params.project_dir, params.cores))

if __name__ == '__main__':
    main()
