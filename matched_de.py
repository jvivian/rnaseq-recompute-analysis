import os
import textwrap
import logging
import subprocess

import errno
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

from pairwise_de import PairwiseAnalysis

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


remove_protein_coding_genes = PairwiseAnalysis.remove_nonprotein_coding_genes


# REFACTOR OUT AT SOME POINT
gencode_path = '/mnt/metadata/gencode.v23.annotation.gtf'


def create_matched_dataframe(df_path, dest):
    """
    Accepts <tissue>_combined.tsv (or combined_pc.tsv) and returns <tissue>_tcga_matched_pairs.tsv

    Checks if the tsv has protein-coding genes in it.

    :param str df_path:
    :param str dest:
    :return:
    """
    df = pd.read_csv(df_path, sep='\t', index_col=0)
    if '_pc.tsv' not in df_path:
        log.info('Removing protein coding genes')
        output_path = os.path.splitext(df_path)[0] + '_pc.tsv'
        df = remove_protein_coding_genes(df, gencode_path, output_path=output_path)
    barcodes = [x[:-3] for x in df.columns]
    matched = list(set(x for x in barcodes if x + '-11' in df.columns and x + '-01' in df.columns))
    # Write out the patients here
    log.info('Writing out patient information to: ' + dest)
    patient_path = os.path.join(dest, 'patients')
    with open(patient_path, 'w') as f:
        f.write('\n'.join([x for x in matched for _ in (0, 1)]))
    matched = [f(x) for x in matched for f in (lambda f1: f1 + '-01', lambda f2: f2 + '-11')]
    matched_df = df[matched]
    log.info('Writing out matched dataframe to: ' + dest)
    match_df_path = os.path.join(dest, 'matched_df.tsv')
    matched_df.to_csv(match_df_path, sep='\t')
    return match_df_path, patient_path


def generate_match_de():
    return textwrap.dedent("""
    library('DESeq2')

    # Argument parsing
    args <- commandArgs(trailingOnly = TRUE)
    df_path <- args[1]
    patient_path <- args[2]
    script.dir <- dirname(sys.frame(1)$ofile)

    # Read in tables / patients
    n <- read.table(df_path, sep='\t', header=1, row.names=1)
    patients <- read.table(patient_path)

    # Create matrix vectors
    disease_vector <- rep(c('T', 'N'), length(patients$V1)/2)
    patient_vector <- patients$V1

    # DESeq2 preprocessing
    countData <- round(n)
    colData <- data.frame(condition=disease_vector, patient=patient_vector, row.names=colnames(countData))
    y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ patient + condition)

    # Run DESeq2
    y <- DESeq(y)
    res <- results(y)
    summary(res)

    # Write out table
    resOrdered <- res[order(res$padj),]
    res_name <- paste(script.dir, 'results.csv', sep='/')
    write.csv(as.data.frame(resOrdered), file=res_name)

    # MA Plot
    ma_name <- paste(script.dir, 'MA.pdf', sep='/')
    pdf(ma_name, width=7, height=7)
    plotMA(res, main='DESeq2')
    dev.off()

    # Dispersion Plot
    disp_name <- paste(script.dir, 'disp.pdf', sep='/')
    pdf(disp_name, width=7, height=7)
    plotDispEsts( y, ylim = c(1e-6, 1e1) )
    dev.off()

    # PVal Hist
    hist_name <- paste(script.dir, 'pval-hist.pdf', sep='/')
    pdf(hist_name, width=7, height=7)
    hist( res$pvalue, breaks=20, col="grey" )
    dev.off()

    # Ratios plots
    qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
    bins <- cut( res$baseMean, qs )
    levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
    ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
    ratio_name <- paste(script.dir, 'ratios.pdf', sep='/')
    pdf(ratio_name, width=7, height=7)
    barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
    dev.off()
    """)


def write_de_script(directory):
    deseq_script_path = os.path.join(directory, 'matched_deseq2.R')
    with open(deseq_script_path, 'w') as f:
        f.write(generate_match_de())
    return deseq_script_path


def run_deseq2(deseq_path, matched_df_path, patient_path):
    """

    :param deseq_path:
    :param matched_df_path:
    :param patient_path:
    :return:
    """
    log.info('Running sample: ' + matched_df_path)
    p = subprocess.Popen(['Rscript', deseq_path, matched_df_path, patient_path])
    out, err = p.communicate()
    if not p.returncode == 0:
        raise RuntimeError('DESeq run failed!\n\n\nstdout:\n{}stderr:\n{}\n\n\n'.format(out, err))
    return 'yay!'


def run_deseq2_wrapper(blob):
    deseq_path, matched_df_path, patient_path = blob
    run_deseq2(deseq_path, matched_df_path, patient_path)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def main():
    tissues = [os.path.join('/mnt/tissues', x) for x in os.listdir('/mnt/tissues')]
    deseq2_info = []
    for tissue in tissues:
        df_path = [os.path.join(tissue, x) for x in os.listdir(tissue) if '_combined' in x]
        if df_path:
            print tissue
            match_dir = os.path.join(tissue, 'matched_TCGA_DE')
            mkdir_p(match_dir)
            deseq_path = write_de_script(match_dir)
            df_path = df_path[0] if len(df_path) == 1 else [x for x in df_path if 'pc.tsv' in x][0]
            print df_path
            match_df_path, patient_path = create_matched_dataframe(df_path, dest=match_dir)
            deseq2_info.append((deseq_path, match_df_path, patient_path))

    with ThreadPoolExecutor(max_workers=8) as executor:
        executor.map(run_deseq2_wrapper, deseq2_info)


if __name__ == '__main__':
    main()
