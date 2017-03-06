import os
import pandas as pd


def pairwise_sig_top_genes(df, tissues, n=100, cutoff=0.90):
    """
    Returns significant genes and top genes for a given n and cutoff

    :param dict df: Dictionary of Dataframe results from a pairwise experiment
    :param list tissues: List of tissues
    :param int n: Number of values to consider "top"
    :param float cutoff: A value between 0 and 1 that determines what percentage of samples must have a
                            significant gene for it to be considered significant within that tissue's cohort.
    :return: Four dictionaries for sig/top genes up and down regulated
    :rtype: tuple(dict, dict, dict, dict)
    """
    sig_genes_up = {}
    sig_genes_down = {}
    top_genes_up = {}
    top_genes_down = {}

    for t in tissues:
        df_up = df[t][df[t].fc > 0]
        df_down = df[t][df[t].fc < 0]

        top_genes_up[t] = set(df_up.head(n).index)
        top_genes_down[t] = set(df_down.head(n).index)

        sig_genes_up[t] = set(df_up[df_up.pval_counts > cutoff].index)
        sig_genes_down[t] = set(df_down[df_down.pval_counts > cutoff].index)

    return sig_genes_up, sig_genes_down, top_genes_up, top_genes_down
