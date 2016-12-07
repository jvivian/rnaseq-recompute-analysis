import logging
import os
import subprocess
import textwrap

import pandas as pd
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def dedupe(items):
    seen = set()
    for item in items:
        if item not in seen:
            yield item
            seen.add(item)


def create_deseq2_inputs(input_dfs, experiment_dir):
    log.info('Reading in dataframes')
    dfs = []
    type_vector = []
    for df_path in tqdm(input_dfs):
        df = pd.read_csv(df_path, sep='\t', index_col=0)
        type_vector.extend(os.path.basename(os.path.split(df_path)[0]) for _ in df.columns)
        dfs.append(df)

    # log.info('Writing out type vector')
    # vector_path = os.path.join(experiment_dir, 'type-vector.txt')
    # if not os.path.exists(vector_path):
    #     with open(vector_path, 'w') as f:
    #         f.write('\n'.join(type_vector))

    log.info('Writing out tissue vector')
    for tissue in dedupe(type_vector):
        vector = [x if x == tissue else 'Not-' + tissue for x in type_vector]
        try:
            with open(os.path.join(experiment_dir, tissue, 'tissue-vector.txt'), 'w') as f:
                f.write('\n'.join(vector))
        except:
            print experiment_dir, tissue

    log.info('Combining and outputting dataframe')
    gtex_tissue = os.path.join(experiment_dir, 'gtex-tissue.tsv')
    if not os.path.exists(gtex_tissue):
        df = pd.concat(dfs, axis=1)
        df.to_csv(gtex_tissue, sep='\t')


def generate_match_de():
    return textwrap.dedent("""
    library('DESeq2')

    # Argument parsing
    args <- commandArgs(trailingOnly = TRUE)
    df_path <- args[1]
    # type_path <- args[2]
    tissue_path <- args[2]

    script.dir <- dirname(tissue_path)

    # Read in tables / vectors
    # design = ~ type + condition
    # The specific tissue is the variable of interest
    # and we're controlling for all the other tissues
    n <- read.table(df_path, sep='\\t', header=1, row.names=1)
    # type <- read.table(type_path)
    tissue <- read.table(tissue_path)

    # Create matrix vectors
    # type_vector <- type$V1 # type
    tissue_vector <- tissue$V1 # condition

    # DESeq2 preprocessing
    # Rounding the countData since DESeQ2 only accepts integer counts
    countData <- round(n)
    colData <- data.frame(tissue=tissue_vector, row.names=colnames(countData))
    y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ tissue)

    # Run DESeq2
    y <- DESeq(y)
    res <- results(y)
    summary(res)

    # Write out table
    resOrdered <- res[order(res$padj),]
    res_name <- paste(script.dir, 'results.tsv', sep='/')
    write.table(as.data.frame(resOrdered), file=res_name, col.names=NA, sep='\\t',  quote=FALSE)

    # MA Plot
    ma_name <- paste(script.dir, 'plots', 'MA.pdf', sep='/')
    pdf(ma_name, width=7, height=7)
    plotMA(res, main='DESeq2')
    dev.off()

    # Dispersion Plot
    disp_name <- paste(script.dir, 'plots', 'dispersion.pdf', sep='/')
    pdf(disp_name, width=7, height=7)
    plotDispEsts( y, ylim = c(1e-6, 1e1) )
    dev.off()

    # PVal Hist
    hist_name <- paste(script.dir, 'plots', 'pval-hist.pdf', sep='/')
    pdf(hist_name, width=7, height=7)
    hist( res$pvalue, breaks=20, col="grey" )
    dev.off()

    # Ratios plots
    qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
    bins <- cut( res$baseMean, qs )
    levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
    ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
    ratio_name <- paste(script.dir, 'plots', 'ratios.pdf', sep='/')
    pdf(ratio_name, width=7, height=7)
    barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
    dev.off()
    """)


def write_de_script(tissue_dir):
    deseq_script_path = os.path.join(tissue_dir, 'deseq2.R')
    with open(deseq_script_path, 'w') as f:
        f.write(generate_match_de())
    return deseq_script_path


def run_deseq2(tissue_dir):
    experiment_dir = os.path.split(tissue_dir)[0]
    deseq_path = os.path.join(tissue_dir, 'deseq2.R')
    df_path = os.path.join(experiment_dir, 'gtex-tissue.tsv')
    # type_vector = os.path.join(experiment_dir, 'type-vector.txt')
    tissue_vector = os.path.join(tissue_dir, 'tissue-vector.txt')
    log.info('Running sample: ' + tissue_vector)
    p = subprocess.Popen(['Rscript', deseq_path, df_path, tissue_vector]) #, stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()
    if not p.returncode == 0:
        raise RuntimeError('DESeq run failed for {}!\n\n\nstdout:\n{}stderr:\n{}\n\n\n'.format(tissue_vector, out, err))
    else:
        log.info('Sample has finished successfully: ' + tissue_vector)
    return 'yay!'
