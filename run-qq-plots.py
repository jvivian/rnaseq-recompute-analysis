import os
import textwrap
import logging
import argparse

import subprocess

from concurrent.futures import ThreadPoolExecutor

from utility_functions import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def generate_qq():
    return textwrap.dedent("""
    library(MASS)

    args <- commandArgs(trailingOnly = TRUE)
    df_path <- args[1]
    output_dir <- args[2]

    n <- read.table(df_path, sep='\t', row.names = 1, header=1)

    tcga <- n[,grepl('TCGA', names(n))]
    gtex <- n[,grepl('GTEX', names(n))]

    for(i in 1:5){
        # Select random gene
        random_gene <- row.names(n[sample(nrow(n), 1),])

        # QQ Plots
        pdf_name <- paste('qq-', random_gene, '.png', sep='')
        qq_name <- paste(output_dir, pdf_name, sep='/')
        png(qq_name, width=7, height=7, units='in', res=300)
        # GTEx
        par(mfrow=c(1,2))
        data <- unlist(round(gtex[random_gene,]))
        params = fitdistr(data, "Negative Binomial")
        plot(pnbinom(sort(data),  size=params$estimate['size'], mu=params$estimate['mu']), ppoints(data), xlab='Theoretical Distribution', ylab='Actual Distribution', main='GTEx')
        abline(0,1)

        # TCGA
        data <- unlist(round(tcga[random_gene,]))
        params = fitdistr(data, "Negative Binomial")
        plot(pnbinom(sort(data),  size=params$estimate['size'], mu=params$estimate['mu']), ppoints(data), xlab='Theoretical Distribution', ylab='Actual Distribution', main='TCGA')
        abline(0,1)

        dev.off()
    }""")


def write_de_script(directory):
    script_path = os.path.join(directory, 'qq-plots.R')
    with open(script_path, 'w') as f:
        f.write(generate_qq())
    return script_path


def run_qq(inputs):
    input_df, output_dir = inputs
    log.info('Plotting GTEx/TCGA Q-Q plots for: ' + os.path.basename(output_dir))
    write_de_script(output_dir)
    script_path = os.path.join(output_dir, 'qq-plots.R')
    p = subprocess.Popen(['Rscript', script_path, input_df, output_dir])
    out, err = p.communicate()
    if not p.returncode == 0:
        raise RuntimeError('plotting failed!\n\n\nstdout:\n{}stderr:\n{}\n\n\n'.format(out, err))
    return 'yay!'


def main():
    """
    First run create_project.py

    Runs DESeq2 for all TCGA samples with a tumor-normal pair
    :return:
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna-seq-directory', type=str, required=True,
                        help='Full path to rna-seq_analysis directory, created by running create_project.py')
    parser.add_argument('--cores', type=int, help='Number of cores to use when running R.', default=1)
    params = parser.parse_args()
    # Create subdirectories
    root_dir = params.rna_seq_directory
    experiment_dir = os.path.join(root_dir, 'experiments/qq-plots')
    tissue_pair_dir = os.path.join(root_dir, 'data/tissue-pairs')
    tissue_pairs = os.listdir(tissue_pair_dir)
    [mkdir_p(os.path.join(experiment_dir, x)) for x in tissue_pairs]
    # Create inputs for plotting
    input_dfs = [os.path.join(tissue_pair_dir, x, 'combined-gtex-tcga-counts-protein-coding.tsv') for x in tissue_pairs]
    output_dirs = [os.path.join(experiment_dir, x) for x in tissue_pairs]
    # Run plots
    with ThreadPoolExecutor(max_workers=params.cores) as executor:
        executor.map(run_qq, zip(input_dfs, output_dirs))


if __name__ == '__main__':
    main()
