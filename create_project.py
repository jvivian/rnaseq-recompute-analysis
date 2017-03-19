#!/usr/bin/env python2.7
"""
John Vivian
October, 2016
"""
import argparse
import logging
import os
import shutil
import sys

import synapseclient
from concurrent.futures import ThreadPoolExecutor
from synapseclient.exceptions import SynapseHTTPError
from tqdm import tqdm

from preprocessing.tissue_preprocessing import create_subframes, concat_frames, remove_nonprotein_coding_genes
from utils import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Synapse inputs
gtex_counts = 'syn7434140'
tcga_counts = 'syn7434253'
gtex_metadata = 'syn7248852'
tcga_metadata = 'syn7248855'
gencode_metadata = 'syn7248851'
paired_table = 'syn7541065'
gene_map = 'syn7801438'
# Subdirectories to create
leaves = ['data/xena-tables/gtex', 'data/xena-tables/tcga', 'data/tissue-pairs',
          'data/tissue-dataframes', 'metadata', 'experiments']


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
    log.info('Downloading input data')
    metadata_dir = os.path.join(root_dir, 'metadata')
    download_information = [(syn, gtex_counts, os.path.join(root_dir, 'data/xena-tables/gtex')),
                            (syn, tcga_counts, os.path.join(root_dir, 'data/xena-tables/tcga')),
                            (syn, gtex_metadata, metadata_dir),
                            (syn, tcga_metadata, metadata_dir),
                            (syn, gencode_metadata, metadata_dir),
                            (syn, paired_table, metadata_dir),
                            (syn, gene_map, metadata_dir)]
    with ThreadPoolExecutor(max_workers=cores) as executor:
        executor.map(synpase_download, download_information)


def create_paired_tissues(root_dir):
    """
    Creates paired tissue dataframes

    :param str root_dir: Project root directory
    """
    log.info('Creating paired tissues')
    with open(os.path.join(root_dir, 'metadata/tissue-pairings.tsv'), 'r') as f:
        for line in tqdm(f):
            if line:
                dirname, gtex, tcga = line.strip().split('\t')
                # Some GTEx tissues need to be combined in the final dataframe
                gtex = gtex.split(',') if ',' in gtex else [gtex]
                tissue_dir = os.path.join(root_dir, 'data/tissue-pairs', dirname)
                combined_path = os.path.join(tissue_dir, 'combined-gtex-tcga-counts.tsv')
                if not os.path.exists(combined_path):
                    mkdir_p(tissue_dir)
                    gtex_dfs = [os.path.join(root_dir, 'data/tissue-dataframes/', g) for g in gtex]
                    tcga_df = os.path.join(root_dir, 'data/tissue-dataframes/', tcga)
                    # Create combined dataframe and group tissues together
                    concat_frames(gtex_df_paths=gtex_dfs, tcga_df_path=tcga_df, output_path=combined_path)
                    # Copy input dataframe NAMES over for clarity
                    for name in gtex_dfs + [tcga_df]:
                        with open(os.path.join(tissue_dir, os.path.basename(name)), 'w') as f:
                            f.write('\n')
                    # Create dataframe of just protein-coding genes
                    gencode_path = os.path.join(root_dir, 'metadata/gencode.v23.annotation.gtf')
                    remove_nonprotein_coding_genes(df_path=combined_path, gencode_path=gencode_path)


def synpase_download(blob):
    """Map function for downloading from Synapse"""
    syn, syn_id, location = blob
    syn.get(syn_id, downloadLocation=location)


def main():
    """
    Recreates the RNA-seq Recompute Analysis project structure.
        - Downloads input data / metadata from Synapse
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
    parser.add_argument('--location', type=str, help='Directory to create project.')
    parser.add_argument('--username', type=str, help='Synapse username (email). Create account at Synpase.org and set '
                                                     'the password in the environment variable "SYNAPSE_PASSWORD".')
    parser.add_argument('--cores', type=int, help='Number of cores to use when running R.', default=1)
    parser.add_argument('--no-download', action='store_true', help='Flag for disabling downloading from Synapse')
    params = parser.parse_args()

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    log.info('Creating project directory tree at: ' + params.location)
    root_dir = os.path.join(params.location, 'rna-seq-analysis')
    [mkdir_p(os.path.join(root_dir, x)) for x in leaves]

    if not params.no_download:
        download_input_data(root_dir=root_dir, user_name=params.username, cores=params.cores)
    else:
        log.info('--no-download enabled, skipping downloads from Synapse. Warning: missing files'
                 'will cause failures downstream.')

    gtex_metadata_path = os.path.join(root_dir, 'metadata/gtex-table.txt')
    tcga_metadata_path = os.path.join(root_dir, 'metadata/tcga-summary.tsv')
    gtex_xena_path = os.path.join(root_dir, 'data/xena-tables/gtex/gtex_gene_counts')
    tcga_xena_path = os.path.join(root_dir, 'data/xena-tables/tcga/tcga_gene_counts')
    tissue_dataframe_path = os.path.join(root_dir, 'data/tissue-dataframes')
    log.info('Creating tissue dataframes at: ' + tissue_dataframe_path)
    create_subframes(gtex_metadata=gtex_metadata_path, tcga_metadata=tcga_metadata_path,
                     tcga_expression=tcga_xena_path, gtex_expression=gtex_xena_path, output_dir=tissue_dataframe_path)
    # Create paired tissue directories
    create_paired_tissues(root_dir)


if __name__ == '__main__':
    main()
