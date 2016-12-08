import logging
import os
import subprocess
from subprocess import PIPE
import textwrap

import pandas as pd

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def create_matched_dataframe(input_df, output_df):
    df = pd.read_csv(input_df, sep='\t', index_col=0)
    barcodes = [x[:-3] for x in df.columns]
    matched = list(set(x for x in barcodes if x + '-11' in df.columns and x + '-01' in df.columns))
    if matched:
        # Write out the patients here
        patient_path = os.path.join(os.path.dirname(output_df), 'patient-pairs.txt')
        with open(patient_path, 'w') as f:
            f.write('\n'.join([x for x in matched for _ in (0, 1)]))
        matched = [f(x) for x in matched for f in (lambda f1: f1 + '-01', lambda f2: f2 + '-11')]
        matched_df = df[matched]
        matched_df.to_csv(output_df, sep='\t')
        return True
    else:
        log.info(input_df + ' no matching samples found beetween barcodes: -01 and -11.')
        return False


def generate_match_de():
    return textwrap.dedent("""
    library('DESeq2')

    # Argument parsing
    args <- commandArgs(trailingOnly = TRUE)
    df_path <- args[1]
    patient_path <- args[2]
    script.dir <- dirname(df_path)

    # Read in tables / patients
    n <- read.table(df_path, sep='\\t', header=1, row.names=1)
    patients <- read.table(patient_path)

    # Create matrix vectors
    disease_vector <- rep(c('T', 'N'), length(patients$V1)/2)
    patient_vector <- patients$V1

    # DESeq2 preprocessing
    # Rounding the countData since DESeQ2 only accepts integer counts
    # The design matrix is conditioned on the two vectors: patient and condition
    countData <- round(n)
    colData <- data.frame(condition=disease_vector, patient=patient_vector, row.names=colnames(countData))
    y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ patient + condition)

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


def write_de_script(directory):
    deseq_script_path = os.path.join(directory, 'deseq2.R')
    with open(deseq_script_path, 'w') as f:
        f.write(generate_match_de())
    return deseq_script_path


def run_deseq2(tissue_dir):
    tissue_name = os.path.basename(tissue_dir)
    log.info('Running sample: ' + tissue_name)
    patient_path = os.path.join(tissue_dir, 'patient-pairs.txt')
    df_path = os.path.join(tissue_dir, 'matched-tcga-counts.tsv')
    deseq_path = os.path.join(tissue_dir, 'deseq2.R')
    p = subprocess.Popen(['Rscript', deseq_path, df_path, patient_path], stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()
    if not p.returncode == 0:
        raise RuntimeError('DESeq run failed for {}!\n\n\nstdout:\n{}stderr:\n{}\n\n\n'.format(tissue_name, out, err))
    else:
        log.info('Sample has finished successfully: ' + tissue_name)
    return 'yay!'
