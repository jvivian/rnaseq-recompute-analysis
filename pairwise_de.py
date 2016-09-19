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


class Analysis(object):

    def __init__(self, tissue_df, cores, gene_map, gencode_path='/mnt/gencode.v23.annotation.gtf'):
        self.df_path = tissue_df
        self.cores = cores
        self.gene_map = gene_map
        self.gencode_path = gencode_path
        self.tissue_dir = os.path.dirname(self.df_path)
        self.pairwise_dir = os.path.join(self.tissue_dir, 'pairwise_results')
        self.results_dir = os.path.join(self.tissue_dir, 'results')
        self.edger_script = os.path.join(self.tissue_dir, 'edgeR-pairwise-DE.R')
        self.pc_path = os.path.join(self.tissue_dir, os.path.splitext(os.path.basename(self.df_path))[0] + '_pc.tsv')
        # Read in dataframe and store tcga_names
        self.df = pd.read_csv(self.df_path, sep='\t', index_col=0)
        self.tcga_names = [name.replace('-', '.') for name in self.df.columns if 'TCGA' in name]

        self.ranked = pd.DataFrame()
        self.genes = {}

    def process_combined_df(self):
        # Read in tissue dataframe and collect TCGA sample names
        log.info('Reading in tissue dataframe: ' + self.df_path)
        self.df = self._remove_nonprotein_coding_genes()
        self.df.to_csv(self.pc_path, sep='\t')
        return self.df

    def _remove_nonprotein_coding_genes(self):
        log.info('Creating dataframe with non-protein coding genes removed.')
        pc_genes = set()
        with open(self.gencode_path, 'r') as f:
            for line in f.readlines():
                if not line.startswith('#'):
                    line = line.split()
                    if line[line.index('gene_type') + 1] == '"protein_coding";':
                        pc_genes.add(line[line.index('gene_id') + 1].split('"')[1])
        pc_genes = list(pc_genes)
        return self.df.ix[pc_genes]

    def run_pairwise_edger(self):
        # Write out edgeR script
        with open(self.edger_script, 'w') as f:
            f.write(self._generate_edger_script())

        # Multiprocess edgeR
        log.info('Beginning pairwise differential expression: Using {} cores'.format(self.cores))
        pool = Pool(self.cores)
        pool.map(self._run_edger, self.tcga_names)

    def _run_edger(self, sample):
        """
        Function used in pool.map to run Rscript

        :param str sample: TCGA sample
        """
        print 'Running sample: ' + sample
        try:
            subprocess.check_call(['Rscript', self.edger_script, sample])
        except Exception as e:
            print e
            return e.message
        return 'yay!'

    def read_results(self):
        # Read in result tables
        log.info('Compiling pairwise differntial r')
        self.ranked = self._rank_results()
        self.ranked = self._add_mapped_genes()
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # pval / fc/ cpm histogram
        f, axarr = plt.subplots(3, figsize=(8, 6))
        sns.distplot(self.ranked.pval, ax=axarr[0], color='r')
        sns.distplot(self.ranked.fc, ax=axarr[1], color='b')
        sns.distplot(self.ranked.cpm, ax=axarr[2])
        axarr[0].set_ylabel('P-values')
        axarr[1].set_ylabel('Log FC')
        axarr[2].set_ylabel('Log CPM')
        plt.savefig(os.path.join(self.results_dir, 'pval_fc_cpm_histogram.pdf'))

    def _rank_results(self, cutoff=.90):
        pvals, fc, cpm = defaultdict(list), defaultdict(list), defaultdict(list)

        # Open each file and append values
        for f in tqdm(os.listdir(self.pairwise_dir)):
            df = pd.read_csv(os.path.join(self.pairwise_dir, f), sep='\t', index_col=0)
            for gene in df.index:
                pvals[gene].append(df.loc[gene]['PValue'])
                fc[gene].append(df.loc[gene]['logFC'])
                cpm[gene].append(df.loc[gene]['logCPM'])

        # Only retain genes that are present in 90% of samples, as we care about genes that matter at a cohort level.
        num_samples = len(os.listdir(self.pairwise_dir))
        for gene in pvals.keys():
            if len(pvals[gene]) < int(num_samples * cutoff):
                pvals.pop(gene)
                fc.pop(gene)
                cpm.pop(gene)
        len(pvals.keys())

        # Rank genes
        self.genes = pvals.keys()
        self.ranked['pval'] = [np.mean(pvals[x]) for x in self.genes]
        self.ranked['pval_counts'] = [sum([1 for y in pvals[x] if y < 0.001]) for x in self.genes]
        self.ranked['pval_std'] = [np.std(pvals[x]) for x in self.genes]
        self.ranked['fc'] = [np.mean(fc[x]) for x in self.genes]
        self.ranked['fc_std'] = [np.std(fc[x]) for x in self.genes]
        self.ranked['cpm'] = [np.mean(cpm[x]) for x in self.genes]
        self.ranked['cpm_std'] = [np.std(cpm[x]) for x in self.genes]
        self.ranked['num_samples'] = [len(pvals[x]) for x in self.genes]
        self.ranked['gene'] = self.genes
        self.ranked.index = self.genes
        self.ranked.sort_values('pval_counts', inplace=True, ascending=False)
        self.ranked.to_csv(os.path.join(self.tissue_dir, 'ranked_df.tsv'), sep='\t')

        return self.ranked

    def _add_mapped_genes(self):
        id_map = pd.read_table(self.gene_map, sep='\t')
        gene_mappings = {x: y for x, y in itertools.izip(id_map['geneId'], id_map['geneName'])}
        mapped_genes = []
        for gene in self.genes:
            try:
                new_gene = gene_mappings[gene]
            except KeyError:
                new_gene = gene
            mapped_genes.append(new_gene)
        self.ranked['gene_name'] = mapped_genes
        return self.ranked

    def _generate_edger_script(self):
        mds_dir = os.path.join(self.tissue_dir, 'plots_mds')
        ma_dir = os.path.join(self.tissue_dir, 'plots_ma')

        for d in [mds_dir, ma_dir, self.pairwise_dir]:
            mkdir_p(d)

        return textwrap.dedent("""
        plot_MDS <- function(name, f){{
            output_name = paste('{mds_dir}/', name, '.pdf', sep='')
            pdf( output_name , width = 7 , height = 7 )
            title = paste('MDS Plot for', name, '\\nNumber of Genes:', dim(f)[1])
            plotMDS( f , main = title, pch = c(rep(21, gtex_count), rep(25, tcga_count)), cex=0.75,
                    col=c(rep('blue', gtex_count), rep('red', tcga_count)))
            dev.off()
            return(f)
        }}

        # Belongs in a .R function for batching
        library(edgeR)
        args <- commandArgs(trailingOnly = TRUE)
        sample <- args[1]
        n <- read.table('{pc_path}', sep='\\t', header=1, row.names=1)
        gtex_count <- length(colnames(n[grepl('GTEX', colnames(n))]))
        tcga_count <- length(colnames(n)) - gtex_count
        cat('Sample: ', sample, '\\t')
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
        gtex_filter <- rowMeans(cpm(gtex)>1) >= gtex_count * .90
        tcga_filter <- rowMeans(cpm(tcga[sample])>1) >= 1
        filter <- gtex_filter + tcga_filter
        filter <- data.frame(filter)
        filter$filter <- as.logical(filter$filter)
        # Add sample from TCGA to GTEx DF
        gtex[sample] <- tcga[sample]
        # Apply filter
        gtex <- gtex[filter$filter,]
        cat('Genes after filtering: ', dim(gtex)[1], '\\n')
        y <- DGEList(counts=gtex, group=group)

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
        write.table(nc, output_name, quote=FALSE, sep='\\t', col.names=NA)

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
        """.format(tissue_dir=self.tissue_dir, mds_dir=mds_dir, ma_dir=ma_dir, pairwise_dir=self.pairwise_dir,
                   pc_path=self.pc_path))


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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
    # params = parser.parse_args()
    pass
