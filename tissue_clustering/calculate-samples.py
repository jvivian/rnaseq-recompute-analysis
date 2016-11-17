import os

import argparse
import pandas as pd


def calculate_samples(df_dir):
    """
    Calculates number of samples for tissues in tissue-pairings.tsv

    :param str df_dir: Path to the tissue-dataframes directory (rna-seq-analysis/data/tissue-dataframes)
    :return: Dataframe containing counts
    :rtype: pd.DataFrame
    """
    # Collect tissue pairing information
    tissues = {}
    with open('tissue-pairings.tsv', 'r') as f:
        for line in f:
            if line:
                dirname, gtex, tcga = line.strip().split('\t')
                gtex = gtex.split(',') if ',' in gtex else [gtex]
                tissues[dirname] = [gtex, tcga]

    # Calculate number of samples based on headers
    gtex_counts, tumor_counts, normal_counts = [], [], []
    for tissue in sorted(tissues.keys()):
        gtex, tcga = tissues[tissue]
        gtex_count, tcga_count, tcga_norm = 0, 0, 0
        for g in gtex:
            with open(os.path.join(df_dir, g), 'r') as f:
                line = f.readline().strip().split('\t')
                gtex_count += len(line)
        with open(os.path.join(df_dir, tcga), 'r') as f:
            line = f.readline().strip().split('\t')
            tumor = [l for l in line if l.endswith('-01')]
            normal = [l for l in line if l.endswith('-11')]
            tcga_count += len(tumor)
            tcga_norm += len(normal)

        # Store counts
        for counts, count in zip([gtex_counts, tumor_counts, normal_counts], [gtex_count, tcga_count, tcga_norm]):
            counts.append(count)

    # Create TSV
    df = pd.DataFrame()
    df['gtex'] = gtex_counts
    df['tcga_tumor'] = tumor_counts
    df['tcga_normal'] = normal_counts
    df.index = sorted(tissues.keys())
    df.to_csv('sample-counts.tsv', sep='\t')
    return df


def main():
    """Calculates number of samples for each tissue"""
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna-seq-directory', type=str, required=True,
                        help='Full path to rna-seq_analysis directory, created by running create_project.py')
    params = parser.parse_args()
    root_dir = params.rna_seq_directory

    calculate_samples(os.path.join(root_dir, 'data/tissue-dataframes'))

if __name__ == '__main__':
    main()
