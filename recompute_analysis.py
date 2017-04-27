# coding: utf-8
import argparse
import logging
import sys

from experiments.cbnt_clustering import CBNTClustering
from experiments.deseq2_time_test import DESeq2TimeTest
from experiments.gtex_vs_tcga import GTExVsTCGA
from experiments.partition_pearson_relationship import PartitionPearsonRelationship
from experiments.tcga_matched_negative_control import TcgaMatchedNegativeControl
from experiments.tcga_tumor_vs_normal import TcgaTumorVsNormal
from experiments.tcga_tumor_vs_normal_negative_control import TcgaNegativeControl
from experiments.tcga_tumor_vs_normal_paired import TcgaMatched
from experiments.tissue_clustering import TissueClustering
from utils import cls, title_tcga_matched

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def runner(instance):
    instance.setup()
    instance.run_experiment()
    instance.teardown()


def main():
    """
    Launchpoint for all experiments associated with the CGL RNA-seq recompute analysis

    WARNING: Some of these experiments take a rather long time to run (days / weeks)

    First run create_project.py (or create-project entrypoint) to create the necessary
    prerequisite project directory
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # TCGA Tumor Vs Normal
    parser_tcga = subparsers.add_parser('tcga-tumor-vs-normal', help='Run TCGA T/N Analysis')
    parser_tcga.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    parser_tcga.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # GTEX (Normal) Vs TCGA Tumor
    parser_gtex_tcga = subparsers.add_parser('gtex-vs-tcga', help='Run GTEX normals against TCGA tumor samples')
    parser_gtex_tcga.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    parser_gtex_tcga.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # TCGA Matched
    parser_tcga_matched = subparsers.add_parser('tcga-matched', help='Run TCGA T/N analysis, matching T/N patients.')
    parser_tcga_matched.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    parser_tcga_matched.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # Negative controls
    tcga_neg = subparsers.add_parser('tcga-neg-control', help='Performs negative control experiment by randomizing the '
                                                              'order of the disease vector relative to the input.')
    tcga_neg.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    tcga_neg.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    tcga_match_neg = subparsers.add_parser('tcga-matched-neg-control',
                                           help='Performs negative control experiment by randomizing disease'
                                                'and patient vector relative to the input.')
    tcga_match_neg.add_argument('--project-dir', help='Full path to project dir (rna-seq-analysis')
    tcga_match_neg.add_argument('--cores', required=True, type=int, help='Number of cores to utilize during run.')

    # Tissue Pair Clustering
    parser_tissue_clustering = subparsers.add_parser('tissue-clustering',
                                                     help='Run PCA / t-SNE clustering to examine similarity between'
                                                          ' TCGA tumor, normal, and GTEx samples.')
    parser_tissue_clustering.add_argument('--project-dir', required=True,
                                          help='Full path to project dir (rna-seq-analysis')

    # CBNT Clustering
    parser_cbnt_clustering = subparsers.add_parser('cbnt-clustering',
                                                   help='Run PCA / t-SNE clustering to examine similarity between'
                                                        ' TCGA tumor, normal, and GTEx samples.'
                                                        ' ComBat Nonparametric with Tissue Covariates')
    parser_cbnt_clustering.add_argument('--project-dir', required=True,
                                        help='Full path to project dir (rna-seq-analysis')

    # DeSeq2 Time Test
    parser_deseq2 = subparsers.add_parser('deseq2-time-test', help='Runs DeSeq2 with increasing number of samples and '
                                                                   'records how long it takes to run')
    parser_deseq2.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis)')
    parser_deseq2.add_argument('--cores', required=True, help='Number of cores to utilize during run.')

    # Partition PearsonR Relationship
    parser_ppr = subparsers.add_parser('partition-pearson', help='Collects individual partition results '
                                                                 'to determine the pearsonR gain'
                                                                 'when combining random groups')
    parser_ppr.add_argument('--project-dir', required=True, help='Full path to project dir (rna-seq-analysis)')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        cls()
        parser.print_help()
        sys.exit(1)

    # Execution
    params = parser.parse_args()
    cls()

    if params.command == 'tissue-clustering':
        log.info('Tissue Clustering')
        runner(TissueClustering(params.project_dir))

    elif params.command == 'partition-pearson':
        log.info('Saves individual partition results for downstream analysis')
        runner(PartitionPearsonRelationship(params.project_dir))

    elif params.command == 'cbnt-clustering':
        log.info('Clustering: ComBat nonparametric with tissue covariates')
        runner(CBNTClustering(params.project_dir))

    elif params.command == 'tcga-matched':
        log.info(title_tcga_matched())
        runner(TcgaMatched(params.project_dir, params.cores))

    elif params.command == 'tcga-tumor-vs-normal':
        log.info('TCGA Tumor Vs Normal')
        runner(TcgaTumorVsNormal(params.project_dir, params.cores))

    elif params.command == 'tcga-neg-control':
        log.info('TCGA Tumor vs Normal Negative Control')
        runner(TcgaNegativeControl(params.project_dir, params.cores))

    elif params.command == 'tcga-matched-neg-control':
        log.info('TCGA Matched Negative Control')
        runner(TcgaMatchedNegativeControl(params.project_dir, params.cores))

    elif params.command == 'deseq2-time-test':
        log.info('DESeq2 Time Test')
        runner(DESeq2TimeTest(params.project_dir, params.cores))

    elif params.command == 'gtex-vs-tcga':
        log.info('GTEx vs. TCGA')
        runner(GTExVsTCGA(params.project_dir, params.cores))


if __name__ == '__main__':
    main()
