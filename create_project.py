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

import pandas as pd
import synapseclient
from concurrent.futures import ThreadPoolExecutor
from synapseclient.exceptions import SynapseHTTPError
from tqdm import tqdm

from preprocessing.tissue_preprocessing import create_subframes, concat_frames, remove_nonprotein_coding_genes
from utility_functions import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# Synapse inputs
gtex_counts = 'syn7434140'
tcga_counts = 'syn7434253'
gtex_metadata = 'syn7248852'
tcga_metadata = 'syn7248855'
gencode_metadata = 'syn7248851'
paired_table = 'syn7437125'
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
    download_information = [(syn, gtex_counts, os.path.join(root_dir, 'data/xena-tables/gtex')),
                            (syn, tcga_counts, os.path.join(root_dir, 'data/xena-tables/tcga')),
                            (syn, gtex_metadata, os.path.join(root_dir, 'metadata')),
                            (syn, tcga_metadata, os.path.join(root_dir, 'metadata')),
                            (syn, gencode_metadata, os.path.join(root_dir, 'metadata')),
                            (syn, paired_table, os.path.join(root_dir, 'metadata'))]
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
                tissue_dir = os.path.join(root_dir, 'data/tissue-pairs', dirname)
                mkdir_p(tissue_dir)
                gtex_df = os.path.join(root_dir, 'data/tissue-dataframes/', gtex)
                tcga_df = os.path.join(root_dir, 'data/tissue-dataframes/', tcga)
                # Create combined dataframe and group tissues together
                combined_path = os.path.join(tissue_dir, 'combined-gtex-tcga-counts.tsv')
                concat_frames(gtex_df_path=gtex_df, tcga_df_path=tcga_df, output_path=combined_path)
                shutil.copy(gtex_df, os.path.join(tissue_dir, gtex))
                shutil.copy(tcga_df, os.path.join(tissue_dir, tcga))
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

    REQUIRED: Your Synapse password must be stored in the environment variable: SYNAPSE_PASSWORD
    e.g.
    $ export SYNAPSE_PASSWORD=foobar
    $ python create_project.py --location /home/ubuntu/
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--location', type=str, help='Directory to create project.')
    parser.add_argument('--username', type=str, help='Synapse username. Create account at Synpase.org and set the '
                                                     'password in the environment variable "SYNAPSE_PASSWORD".')
    parser.add_argument('--cores', type=int, help='Number of cores to use when running R.', default=1)
    params = parser.parse_args()

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    log.info('Creating project directory tree at: ' + params.location)
    root_dir = os.path.join(params.location, 'rna-seq-analysis')
    [mkdir_p(os.path.join(root_dir, x)) for x in leaves]

    download_input_data(root_dir=root_dir, user_name=params.username, cores=params.cores)

    gtex_metadata_path = os.path.join(root_dir, 'metadata/gtex-table.txt')
    tcga_metadata_path = os.path.join(root_dir, 'metadata/tcga-summary.tsv')
    gtex_xena_path = os.path.join(root_dir, 'data/xena-tables/gtex/gtex_gene_counts')
    tcga_xena_path = os.path.join(root_dir, 'data/xena-tables/tcga/tcga_gene_counts')
    tissue_dataframe_path = os.path.join(root_dir, 'data/tissue-dataframes')
    log.info('Creating tissue dataframes at: ' + tissue_dataframe_path)
    create_subframes(gtex_metadata=gtex_metadata_path, tcga_metadata=tcga_metadata_path,
                     tcga_expression=tcga_xena_path, gtex_expression=gtex_xena_path, output_dir=tissue_dataframe_path)
    # Deal with cervix - special case where we're combining the GTEx tissue to match TCGA
    c1 = pd.read_csv(os.path.join(tissue_dataframe_path, 'Cervix_Ectocervix.tsv'))
    c2 = pd.read_csv(os.path.join(tissue_dataframe_path, 'Cervix_Endocervix.tsv'))
    df = pd.concat([c1, c2], axis=1)
    df.to_csv(os.path.join(tissue_dataframe_path, 'Cervix_endo_and_ecto.tsv'), sep='\t')
    # Create paired tissue directories
    create_paired_tissues(root_dir)


if __name__ == '__main__':
    main()
