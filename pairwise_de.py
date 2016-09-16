#!/usr/bin/env python2.7
"""
John Vivian
September, 2016
"""
import argparse
import errno
import itertools
import logging
import os
import subprocess
import textwrap
from collections import defaultdict
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

log = logging.getLogger(__name__)
# See main for variable explanation
pc_df_path = None


def remove_nonprotein_coding_genes(df, gencode_path='/mnt/gencode.v23.annotation.gtf'):
    pc_genes = set()
    with open(gencode_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                line = line.split()
                if line[line.index('gene_type') + 1] == '"protein_coding";':
                    pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
    pc_genes = list(pc_genes)
    return df.ix[pc_genes]


def run_edger(sample):
    print 'Running sample: ' + sample
    try:
        subprocess.check_call(['Rscript', pc_df_path, sample])
    except Exception as e:
        print e
        return e.message
    return 'yay!'


# TODO: Get CPM cutoff to make sense... 25% IQR for CPM didn't really work
def generate_edger_script(tissue_df):
    tissue_dir = os.path.dirname(tissue_df)
    mds_dir = os.path.join(tissue_dir, 'plots_mds')
    ma_dir = os.path.join(tissue_dir, 'plots_ma')
    pairwise_dir = os.path.join(tissue_dir, 'pairwise_results')
    for d in [mds_dir, ma_dir, pairwise_dir]:
        mkdir_p(d)

    return textwrap.dedent("""
    plot_MDS <- function(name, f){
        output_name = paste('{mds_dir}/', name, '.pdf', sep='')
        pdf( output_name , width = 7 , height = 7 )
        title = paste('MDS Plot for', name, '\nNumber of Genes:', dim(f)[1])
        plotMDS( f , main = title, pch = c(rep(21, gtex_count), rep(25, tcga_count)), cex=0.75,
                col=c(rep('blue', gtex_count), rep('red', tcga_count)))
        dev.off()
        return(f)
    }

    # Belongs in a .R function for batching
    library(edgeR)
    args <- commandArgs(trailingOnly = TRUE)
    sample <- args[1]
    gtex_count <- length(colnames(n[grepl('GTEX', colnames(n))]))
    tcga_count <- length(colnames(n)) - gtex_count
    n <- read.table('{tissue_df}', sep='\t', header=1, row.names=1)
    cat('Sample: ', sample, '\t')
    # Define GTEx and TCGA Dataframes
    gtex <- n[0:gtex_count]
    tcga <- n[(gtex_count + 1) : (gtex_count + tcga_count)]

    # Filtering protocol: at least 1cpm in 90% of gtex OR >= 1 in the single tcga
    # Smaller CPM thresholds are usually appropriate for larger libraries. As a general
    # rule, a good threshold can be chosen by identifying the CPM that corresponds to a
    # count of 10
    # cutoff = cpm(10, mean(y$samples$lib.size))
    # Construct DGE object, normalize, and produce MDS plot
    group <- c(rep('gtex', gtex_count), rep('tcga', 1 ))
    y <- DGEList(counts=gtex, group=group)
    gtex_filter <- rowMeans(cpm(gtex)>1) >= gtex_count * .90
    tcga_filter <- rowMeans(cpm(tcga[sample])>1) >= 1
    filter <- gtex_filter + tcga_filter
    filter <- data.frame(filter)
    filter$filter <- as.logical(filter$filter)
    # Add sample from TCGA to GTEx DF
    gtex[sample] <- tcga[sample]
    # Apply filter
    gtex <- gtex[filter$filter,]
    cat('Genes after filtering: ', dim(gtex)[1], '\n')

    # Calculate normalization factors
    y <- calcNormFactors(y)
    f <- plot_MDS(sample, y)

    # Establish design and estimate dispersion
    # Currently the design model is a binary group: GTEx / TCGA
    design <- model.matrix(~group)
    f <- estimateDisp( f, design )

    # Generate normalized counts dataframe
    nc <- cpm(f, normalized.lib.sizes=FALSE)
    output_name <- paste('{tissue_dir}/norm_count_tables/', sample, '.tsv', sep='')
    write.table(nc, output_name, quote=FALSE, sep='\t', col.names=NA)

    # Fit the Quasi-Likelihood GLM
    fit <- glmQLFit(f, design)
    qlf <- glmQLFTest(fit, coef=2)

    # Create QLF table
    qlf_sort <- qlf$table[order(qlf$table$logFC),]

    # Find DE Genes
    de_qlf <- rownames(qlf_sort[abs(qlf_sort$logFC) > 2,])
    summary(de <- decideTestsDGE(qlf))
    detags <- rownames(f)[as.logical(de)]

    # Generate MA plot
    title = paste('{ma_dir}/', sample, '.pdf', sep='')
    pdf(title, width=7, height=7)
    plotSmear(qlf, de.tags = detags)
    abline(h=c(-2, 2), col="blue")
    dev.off()

    # Write out table
    output_name <- paste('{pairwise_dir}/', sample, '.tsv', sep='')
    write.table(qlf_sort, output_name, quote=FALSE, sep='\t', col.names=NA)
    """.format(tissue_dir=tissue_dir, mds_dir=mds_dir, ma_dir=ma_dir, pairwise_dir=pairwise_dir, tissue_df=tissue_df))


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def collect_pairwise_tables(fpath, cutoff=.90):
    pvals, fc, cpm = defaultdict(list), defaultdict(list), defaultdict(list)

    # Open each file and append values
    for f in tqdm(os.listdir(fpath)):
        df = pd.read_csv(os.path.join(fpath, f), sep='\t', index_col=0)
        for gene in df.index:
            pvals[gene].append(df.loc[gene]['PValue'])
            fc[gene].append(df.loc[gene]['logFC'])
            cpm[gene].append(df.loc[gene]['logCPM'])

    # Only retain genes that are present in 90% of samples, as we care about genes that matter at a cohort level.
    num_samples = len(os.listdir(fpath))
    for gene in pvals.keys():
        if len(pvals[gene]) < int(num_samples * cutoff):
            pvals.pop(gene)
            fc.pop(gene)
            cpm.pop(gene)
    len(pvals.keys())

    # Rank genes
    ranked = pd.DataFrame()
    genes = pvals.keys()
    ranked['pval'] = [np.mean(pvals[x]) for x in genes]
    ranked['pval_counts'] = [sum([1 for y in pvals[x] if y < 0.001]) for x in genes]
    ranked['pval_std'] = [np.std(pvals[x]) for x in genes]
    ranked['fc'] = [np.mean(fc[x]) for x in genes]
    ranked['fc_std'] = [np.std(fc[x]) for x in genes]
    ranked['cpm'] = [np.mean(cpm[x]) for x in genes]
    ranked['cpm_std'] = [np.std(cpm[x]) for x in genes]
    ranked['num_samples'] = [len(pvals[x]) for x in genes]
    ranked['gene'] = genes
    ranked.index = genes
    ranked.sort_values('pval_counts', inplace=True, ascending=False)

    return ranked, genes


def add_mapped_genes(gene_map, ranked, genes):
    id_map = pd.read_table(gene_map, sep='\t')
    gene_mappings = {x: y for x, y in itertools.izip(id_map['geneId'], id_map['geneName'])}
    mapped_genes = []
    for gene in genes:
        try:
            new_gene = gene_mappings[gene]
        except KeyError:
            new_gene = gene
        mapped_genes.append(new_gene)
    ranked['gene_name'] = mapped_genes
    return ranked


def main():
    """
    RNA-seq Pairwise Differential Expression Methodology:

    1. Read in combined dataframe for a particular tissue
    2. Remove all genes except protein-coding genes. We'll use the gencode annotation file,
        which was used in the recompute, to pull out the set of genes that are protein coding.
        HG38 includes many non-protein coding genes which may affect how normalization is performed.
    3. Run RScript to perform differential expression
        I. Create dataframe of GTEx and single TCGA sample
        II.

    :return:
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tissue-df', type=str, required=True,
                        help='Tissue dataframe: genes/isoforms by sample names. GTEx columns should come first followed'
                             'by TCGA columns.')
    parser.add_argument('--gene-map', type=str, default='/mnt/metadata/attrs.tsv',
                        help='File containing map information. Must have at least 2 column names: '
                             'geneID and geneName.')
    parser.add_argument('--cores', default=8, type=int, help='Number of cores to use.')
    params = parser.parse_args()

    tissue_dir = os.path.dirname(params.tissue_df)

    # Read in tissue dataframe and collect TCGA sample names
    df = pd.read_csv(params.tissue_df, sep='\t', index_col=0)
    tcga_cols = filter(lambda x: 'TCGA' in x, df.columns)
    tcga_cols = [x.replace('-', '.') for x in tcga_cols]

    # Remove non-protein coding genes
    log.info('Creating dataframe with non-protein coding genes removed.')
    pc_df = remove_nonprotein_coding_genes(df)
    global pc_df_path  # run_edger function needs a non-mapped variable passed to it
    pc_df_path = os.path.join(tissue_dir, os.path.splitext(os.path.basename(params.tissue_df)[0], '_pc.tsv'))
    pc_df.to_csv(pc_df_path, sep='\t')

    # Write out edgeR script
    with open(os.path.join(tissue_dir, 'edgeR-pairwise-DE.R'), 'w') as f:
        f.write(generate_edger_script(params.tissue_df))

    # Multiprocess edgeR
    log.info('Beginning pairwise differential expression: Using {} cores'.format(params.cores))
    pool = Pool(params.cores)
    pool.map(run_edger, tcga_cols)

    # Read in result tables
    log.info('Compiling pairwise differntial r')
    ranked, genes = collect_pairwise_tables(fpath=os.path.join(tissue_dir, 'pairwise_results'))
    ranked = add_mapped_genes(params.gene_map, ranked=ranked, genes=genes)
    results_dir = os.path.join(tissue_dir, 'results')
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    # pval / fc/ cpm histogram
    f, axarr = plt.subplots(3, figsize=(8, 6))
    sns.distplot(ranked.pval, ax=axarr[0], color='r')
    sns.distplot(ranked.fc, ax=axarr[1], color='b')
    sns.distplot(ranked.cpm, ax=axarr[2])
    axarr[0].set_ylabel('P-values')
    axarr[1].set_ylabel('Log FC')
    axarr[2].set_ylabel('Log CPM')
    plt.savefig(os.path.join(results_dir, 'pval_fc_cpm_histogram.pdf'))


# TODO: Plots, output results, ranking methods, etc
