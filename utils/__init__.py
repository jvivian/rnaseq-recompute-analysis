# coding: utf-8
import errno
import os
import pickle
import textwrap

import pandas as pd
import logging

import subprocess

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def quantile_normalize(df):
    """
    Quantile normalization is a technique for making two distributions identical in statistical properties.
    Assuming genes by samples
    from: http://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe

    Gene x Sample matrix
    """
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def dedupe(items):
    seen = set()
    for item in items:
        if item not in seen:
            yield item
            seen.add(item)


def add_gene_names(df_path, gene_map_path):
    """
    Adds gene names to results.tsv from DESeq2

    :param str df_path: Path to dataframe from DESeq2
    :param str gene_map_path: Path to a pickled dictionary that maps geneId to geneName
    :return: Dataframe containing the geneNames as the index and the geneId as an appended column
    :rtype: pd.DataFrame
    """
    gene_map = pickle.load(open(gene_map_path, 'rb'))
    df = pd.read_csv(df_path, sep='\t', index_col=0)

    gene_names = [gene_map[x] if x in gene_map.keys() else x for x in df.index]
    df['geneId'] = df.index
    df.index = gene_names

    return df


def write_script(script_func, directory, name='deseq2.R'):
    """
    Writes out script (usually an R script for DESeq2) in a given directory

    :param function script_func: Function that returns text block containing script
    :param str directory: Directory where script will be written
    :param str name: Name of script
    :return: Path to script
    :rtype: str
    """
    deseq_script_path = os.path.join(directory, name)
    log.info('Writing out script: ' + deseq_script_path)
    with open(deseq_script_path, 'w') as f:
        f.write(script_func())
    return deseq_script_path


def run_deseq2(blob):
    """
    Function for running DESeq2 in batches
    Designed for use with ThreadPoolExecutor's map

    :param tuple(str, list[str]) blob:
    """
    script_path, args = blob
    log.debug(str(args))
    rt = subprocess.check_call(['Rscript', script_path] + list(args))
    if not rt == 0:
        raise RuntimeError('Run failed!')
    else:
        log.info('Run has finished successfully')
    return 'yay!'


def cls():
    os.system('cls' if os.name=='nt' else 'clear')


def title_gtex_one_vs_all():
    return(textwrap.dedent("""
     ██████╗████████╗███████╗██╗  ██╗     ██████╗ ███╗   ██╗███████╗    ██╗   ██╗███████╗     █████╗ ██╗     ██╗
    ██╔════╝╚══██╔══╝██╔════╝╚██╗██╔╝    ██╔═══██╗████╗  ██║██╔════╝    ██║   ██║██╔════╝    ██╔══██╗██║     ██║
    ██║  ███╗  ██║   █████╗   ╚███╔╝     ██║   ██║██╔██╗ ██║█████╗      ██║   ██║███████╗    ███████║██║     ██║
    ██║   ██║  ██║   ██╔══╝   ██╔██╗     ██║   ██║██║╚██╗██║██╔══╝      ╚██╗ ██╔╝╚════██║    ██╔══██║██║     ██║
    ╚██████╔╝  ██║   ███████╗██╔╝ ██╗    ╚██████╔╝██║ ╚████║███████╗     ╚████╔╝ ███████║    ██║  ██║███████╗███████╗
     ╚═════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝     ╚═════╝ ╚═╝  ╚═══╝╚══════╝      ╚═══╝  ╚══════╝    ╚═╝  ╚═╝╚══════╝╚══════╝
        """))


def title_clustering():
    return(textwrap.dedent("""
         ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ ██╗███╗   ██╗ ██████╗
        ██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗██║████╗  ██║██╔════╝
        ██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝██║██╔██╗ ██║██║  ███╗
        ██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗██║██║╚██╗██║██║   ██║
        ╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║██║██║ ╚████║╚██████╔╝
         ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝ ╚═════╝
        """))


def title_gtex_pairwise():
    return(textwrap.dedent("""
         ██████╗████████╗███████╗██╗  ██╗    ██████╗  █████╗ ██╗██████╗ ██╗    ██╗██╗███████╗███████╗
        ██╔════╝╚══██╔══╝██╔════╝╚██╗██╔╝    ██╔══██╗██╔══██╗██║██╔══██╗██║    ██║██║██╔════╝██╔════╝
        ██║  ███╗  ██║   █████╗   ╚███╔╝     ██████╔╝███████║██║██████╔╝██║ █╗ ██║██║███████╗█████╗
        ██║   ██║  ██║   ██╔══╝   ██╔██╗     ██╔═══╝ ██╔══██║██║██╔══██╗██║███╗██║██║╚════██║██╔══╝
        ╚██████╔╝  ██║   ███████╗██╔╝ ██╗    ██║     ██║  ██║██║██║  ██║╚███╔███╔╝██║███████║███████╗
         ╚═════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝    ╚═╝     ╚═╝  ╚═╝╚═╝╚═╝  ╚═╝ ╚══╝╚══╝ ╚═╝╚══════╝╚══════╝
        """))


def title_tcga_matched():
    return(textwrap.dedent("""
        ████████╗ ██████╗ ██████╗  █████╗     ███╗   ███╗ █████╗ ████████╗ ██████╗██╗  ██╗███████╗██████╗
        ╚══██╔══╝██╔════╝██╔════╝ ██╔══██╗    ████╗ ████║██╔══██╗╚══██╔══╝██╔════╝██║  ██║██╔════╝██╔══██╗
           ██║   ██║     ██║  ███╗███████║    ██╔████╔██║███████║   ██║   ██║     ███████║█████╗  ██║  ██║
           ██║   ██║     ██║   ██║██╔══██║    ██║╚██╔╝██║██╔══██║   ██║   ██║     ██╔══██║██╔══╝  ██║  ██║
           ██║   ╚██████╗╚██████╔╝██║  ██║    ██║ ╚═╝ ██║██║  ██║   ██║   ╚██████╗██║  ██║███████╗██████╔╝
           ╚═╝    ╚═════╝ ╚═════╝ ╚═╝  ╚═╝    ╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝╚══════╝╚═════╝
    """))


def title_pairwise_gtex_tcga():
    return(textwrap.dedent("""
         ██████╗████████╗███████╗██╗  ██╗    ██╗   ██╗███████╗    ████████╗ ██████╗ ██████╗  █████╗
        ██╔════╝╚══██╔══╝██╔════╝╚██╗██╔╝    ██║   ██║██╔════╝    ╚══██╔══╝██╔════╝██╔════╝ ██╔══██╗
        ██║  ███╗  ██║   █████╗   ╚███╔╝     ██║   ██║███████╗       ██║   ██║     ██║  ███╗███████║
        ██║   ██║  ██║   ██╔══╝   ██╔██╗     ╚██╗ ██╔╝╚════██║       ██║   ██║     ██║   ██║██╔══██║
        ╚██████╔╝  ██║   ███████╗██╔╝ ██╗     ╚████╔╝ ███████║       ██║   ╚██████╗╚██████╔╝██║  ██║
         ╚═════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝      ╚═══╝  ╚══════╝       ╚═╝    ╚═════╝ ╚═════╝ ╚═╝  ╚═╝
    """))
