#!/usr/bin/env python2.7
"""
John Vivian
October, 2016
"""
import argparse
import logging
import os
import sys

import pandas as pd
import synapseclient
from concurrent.futures import ThreadPoolExecutor
from synapseclient.exceptions import SynapseHTTPError

from preprocessing.tissue_preprocessing import filter_nonprotein_coding_genes, create_candidate_pairs, cluster_df
from utils import mkdir_p, cls

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Synapse inputs
inputs = {
    # Expression tables
    'syn7434140': 'data/xena',
    'syn7434253': 'data/xena',
    'syn9998692': 'data/xena',
    # Pickled Objects
    'syn9962451': 'data/objects',
    'syn9962453': 'data/objects',
    'syn9998690': 'data/objects',
    # Additional Metadata
    'syn7248852': 'metadata',
    'syn7248855': 'metadata',
    'syn7248851': 'metadata',
    'syn9962462': 'metadata',
    'syn9998691': 'metadata'}


def download_input_data(root_dir, user_name, cores):
    """
    Downloads input data for the project

    :param str root_dir: Project root directory
    :param str user_name: Synapse user name
    :param int cores: Number of cores to use
    """
    syn = synapseclient.Synapse()
    log.info('Attempting to login to Synapse')
    try:
        syn.login(user_name, os.environ['SYNAPSE_PASSWORD'])
    except KeyError:
        raise RuntimeError('User failed to supply an environment variable: "SYNAPSE_PASSWORD".')
    except SynapseHTTPError as e:
        raise RuntimeError('Failed to connect Synapse client, check password: ' + e.message)
    # Download input tables
    log.info('Downloading data')
    download_information = zip([syn for _ in xrange(len(inputs))],
                               [root_dir for _ in xrange(len(inputs))],
                               inputs.iteritems())

    with ThreadPoolExecutor(max_workers=cores) as executor:
        executor.map(synpase_download, download_information)


def synpase_download(blob):
    """Map function for downloading from Synapse"""
    syn, root_dir, info = blob
    syn_id, location = info
    download_dir = os.path.join(root_dir, location)
    if not os.path.exists(os.path.join(download_dir, syn.get(syn_id, downloadFile=False).name)):
        syn.get(syn_id, downloadLocation=download_dir)


def build(root_dir):
    log.info('Reading in expression dataframes')
    tcga_df = pd.read_csv(os.path.join(root_dir, 'data/xena/tcga_rsem_gene_counts.tsv'), sep='\t', index_col=0)
    gtex_df = pd.read_csv(os.path.join(root_dir, 'data/xena/gtex_rsem_gene_counts.tsv'), sep='\t', index_col=0)
    df = pd.concat([tcga_df, gtex_df], axis=1)
    df.index.name = None

    log.info('Reversing Xena normalization to get raw counts')
    df = df.apply(lambda x: (2 ** x) - 1)

    log.info('Retaining only protein-coding genes')
    gencode_path = os.path.join(root_dir, 'metadata/gencode.v23.annotation.gtf')
    df = filter_nonprotein_coding_genes(df, gencode_path=gencode_path)

    log.info('Filtering out samples with no corresponding metadata')
    metadata = pd.read_csv(os.path.join(root_dir, 'metadata/tcga_gtex_metadata_intersect.tsv'), sep='\t', index_col=0)
    met_samples = list(metadata.index)
    samples = [x for x in df.columns if x in met_samples]
    log.info('Reducing dataframe from {} to {} samples'.format(len(df.columns), len(samples)))
    df = df[samples]

    log.info('Saving Dataframe')
    output_name = os.path.join(root_dir, 'data/xena/tcga_gtex_counts_protein_coding.tsv')
    if not os.path.exists(output_name):
        df.to_csv(output_name, sep='\t')

    log.info('Creating and clustering candidate pairs')
    create_candidate_pairs(df, root_dir)

    output_dir = os.path.join(root_dir, 'data/clustering', 'all')
    mkdir_p(output_dir)
    for cluster_type in ['tissue', 'type', 'dataset']:
        output_path = os.path.join(output_dir, 'tSNE-clustering-{}.html'.format(cluster_type))
        log.info('Clustering entire dataset by: {}'.format(cluster_type))
        cluster_df(df.T, root_dir, output_path=output_path,
                   title='t-SNE Clustering of TCGA and GTEx by {}'.format(cluster_type), colorby='tissue')


def main():
    """
    Recreates the RNA-seq Recompute Analysis project structure.

    Two methods:
    1. Download all files needed for downstream experiments directly from Synapse

    2. Download raw input files and build project from scratch
            - Download input data / metadata from Synapse
        - Creates dataframes for GTEx and TCGA separated by body site or disease name
        - Pairs matching tissues together
        - Creates a subset of the combined dataframes containing only protein-coding genes

    REQUIRED: Your Synapse password must be stored in the environment variable: SYNAPSE_PASSWORD
    e.g.
    $ export SYNAPSE_PASSWORD=foobar
    $ python create_project.py --location=/home/ubuntu/ --username=foo@bar.com

    If you do not have a synapse account, create one for free in under a minute at www.synapse.org
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--location', type=str, help='Directory to create project: "rna-seq-analysis"')
    parser.add_argument('--username', type=str, help='Synapse username (email). Create account at Synpase.org and set '
                                                     'the password in the environment variable "SYNAPSE_PASSWORD".')
    parser.add_argument('--cores', default=2, type=int, help='Number of cores to use.')

    # Subparsers
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('download', help='Recommended. Downloads directly from Synapse')
    subparsers.add_parser('build', help='Builds project from scratch.')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        cls()
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    cls()

    # Create directories
    log.info('Creating project directory tree at: ' + args.location)
    root_dir = os.path.join(args.location, 'rna-seq-analysis')
    leaves = ['data/xena', 'data/objects', 'data/tissue-pairs', 'data/clustering', 'metadata', 'experiments']
    [mkdir_p(os.path.join(root_dir, x)) for x in leaves]

    # Commmand
    if args.command == 'download':
        download_input_data(root_dir, args.username, args.cores)
    elif args.command == 'build':
        download_input_data(root_dir, args.username, args.cores)
        build(root_dir)
    else:
        raise RuntimeError('Argument error, check commands: {}'.format(args))


if __name__ == '__main__':
    main()
