#!/usr/bin/env python2.7
"""
John Vivian
October, 2016
"""
import argparse
import logging
import os
import pickle
import sys

import pandas as pd
import synapseclient
from concurrent.futures import ThreadPoolExecutor
from synapseclient.exceptions import SynapseHTTPError

from preprocessing.tissue_preprocessing import cluster_entire_dataset
from preprocessing.tissue_preprocessing import cluster_tissues
from preprocessing.tissue_preprocessing import create_tissue_pairs, plot_num_samples_per_dataset
from preprocessing.tissue_preprocessing import filter_nonprotein_coding_genes
from preprocessing.tissue_preprocessing import filter_samples_by_metadata
from utils import mkdir_p, cls

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Synapse inputs
inputs = {
    # Expression tables
    'syn7434140': 'data/xena',
    'syn7434253': 'data/xena',
    'syn10002767': 'data/xena',
    # Pickled Objects
    'syn9962451': 'data/objects',
    'syn9962453': 'data/objects',
    'syn9998690': 'data/objects',
    # Additional Metadata
    'syn7248852': 'metadata',
    'syn7248855': 'metadata',
    'syn7248851': 'metadata',
    'syn9962462': 'metadata',
    'syn9998691': 'metadata',
    'syn10131007': 'metadata',
    'syn10142937': 'metadata',
    'syn10142930': 'metadata'}


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
    output_name = os.path.join(root_dir, 'data/xena/tcga_gtex_counts_protein_coding.tsv')
    if not os.path.exists(output_name):
        log.info('Reading in expression dataframes')
        tcga_df = pd.read_csv(os.path.join(root_dir, 'data/xena/tcga_rsem_gene_counts.tsv'), sep='\t', index_col=0)
        gtex_df = pd.read_csv(os.path.join(root_dir, 'data/xena/gtex_rsem_gene_counts.tsv'), sep='\t', index_col=0)
        df = pd.concat([tcga_df, gtex_df], axis=1)
        df.index.name = None

        log.info('Reversing Xena normalization to get raw counts')
        df = df.apply(lambda x: (2 ** x) - 1)

        log.info('Retaining only protein-coding genes')
        df = filter_nonprotein_coding_genes(df, root_dir)

        log.info('Filtering out samples with no corresponding metadata')
        df, samples = filter_samples_by_metadata(df, root_dir)

        log.info('Creating candidate pairs')
        tissues = create_tissue_pairs(df, root_dir)

        log.info('Clustering raw tissue pairs (for comparison)')
        cluster_tissues(df, root_dir, tissues)

        # Clustering of entire Dataset
        cluster_entire_dataset(df, root_dir)

        log.info('Saving Dataframe')
        df.to_csv(output_name, sep='\t')

    # Read in DESeq2 Normalized Dataset
    log.info('\n\nCreating tissue pairs and clusters with DESeq2 normalized expression values')
    df = pd.read_csv(os.path.join(root_dir, 'data/xena/deseq2_normalized_tcga_gtex_counts.tsv'),
                     sep='\t', index_col=0, low_memory=True)

    log.info('Creating tissue pairs from normalized data')
    tissues = create_tissue_pairs(df, root_dir, sub_dir='normalized')

    log.info('Plotting dataset by number of samples')
    plot_num_samples_per_dataset(df, tissues, root_dir)

    log.info('Clustering normalized tissue pairs')
    cluster_tissues(df, root_dir, tissues, sub_dir='normalized')

    # Clustering of entire Dataset
    cluster_entire_dataset(df, root_dir, sub_dir='normalized')

    # Cluster in a reduced gene space using 409 genes from UCSF
    gene_map = pickle.load(open(os.path.join(root_dir, 'data/objects/gene_map.pickle'), 'rb'))
    genes = [gene_map[x] if x in gene_map else x for x in df.index]
    df.index = genes

    ucsf_path = os.path.join(root_dir, 'metadata/UCSF-RNAPanel-Final-412-genes.csv')
    ucsf_genes = [x.strip() for x in open(ucsf_path, 'r').readlines()]
    ucsf_genes = [x for x in ucsf_genes if x in genes]

    log.info('Clustering UCSF gene subset (409 genes)')
    cluster_tissues(df.T[ucsf_genes].T, root_dir, tissues, sub_dir='UCSF-subset')

    cluster_entire_dataset(df.T[ucsf_genes].T, root_dir, sub_dir='UCSF-subset')


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
    leaves = ['data/xena', 'data/objects', 'data/tissue-pairs', 'data/tsne-clustering', 'metadata', 'experiments']
    [mkdir_p(os.path.join(root_dir, x)) for x in leaves]

    # Download, build tissue pairs, and cluster
    download_input_data(root_dir, args.username, args.cores)
    build(root_dir)


if __name__ == '__main__':
    main()
