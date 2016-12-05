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
    p = subprocess.Popen(['Rscript'] + list(args), stderr=subprocess.PIPE)
    out, err = p.communicate()
    if not p.returncode == 0:
        raise RuntimeError('Run failed: {}'.format(err))
    else:
        log.info('Run has finished successfully')
    return 'yay!'


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
