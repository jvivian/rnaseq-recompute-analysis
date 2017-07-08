import logging
import os
import shutil
import textwrap
from subprocess import Popen, PIPE
from experiments.abstractExperiment import AbstractExperiment
from utils import write_script
from random import sample
import pandas as pd
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class DESeq2TimeTest(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(DESeq2TimeTest, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/deseq2-time-test')
        self.vector_dir = os.path.join(self.experiment_dir, 'vectors')
        self.results_dir = os.path.join(self.experiment_dir, 'results')
        self.script_path = None
        self.df_b = os.path.join(self.tissue_pair_dir, 'Breast', 'tcga-gtex-exp.tsv')
        self.df_h = os.path.join(self.tissue_pair_dir, 'Skin-Head', 'tcga-gtex-exp.tsv')
        self.results = os.path.join(self.experiment_dir, 'results.tsv')
        self.df_path = os.path.join(self.experiment_dir, 'arbitary_df.tsv')
        self.num_samples = [1, 50, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]

    def setup(self):
        self.create_directories([self.vector_dir, self.results_dir])
        self.script_path = write_script(self.deseq2_script, directory=self.experiment_dir)

        log.info('Reading in expression dataframes')
        self.df_b = pd.read_csv(self.df_b, index_col=0, sep='\t')
        self.df_h = pd.read_csv(self.df_h, index_col=0, sep='\t')

        log.info('Creating large dataframe to sample from')
        df = pd.concat([self.df_b, self.df_h], axis=1)
        df.to_csv(self.df_path, sep='\t')


        log.info('Creating vectors')
        samples = list(df.columns)

        for i in self.num_samples:
            s = [x.replace('-', '.') for x in sample(samples, i)]
            with open(os.path.join(self.vector_dir, '{}-vector'.format(i)), 'w') as f:
                f.write('\n'.join(s) + '\n')

    def run_experiment(self):
        vectors = [os.path.join(self.vector_dir, x) for x in
                   sorted(os.listdir(self.vector_dir), key=lambda x: float(x.split('-')[0]))]  # Sort by # of samples
        for i, vector in zip(self.num_samples, vectors):
            log.info('Running {} number of samples'.format((i)))
            err = self.run_deseq2(vector)
            with open(self.results, 'a') as f:
                f.write('{}\t{}\n'.format(vector, err))

    def run_deseq2(self, vector):
        p = Popen(['time', 'Rscript', self.script_path] + [vector], stderr=PIPE, stdout=PIPE)
        out, err = p.communicate()
        log.info('Out: ' + out + '\n\n')
        log.info('Err: ' + err + '\n\n')
        return err

    def deseq2_script(self):
        return textwrap.dedent("""
            suppressMessages(library('DESeq2'))
            suppressMessages(library('data.table'))
            suppressMessages(library('BiocParallel'))
            register(MulticoreParam({cores}))

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- '{df_path}'
            vector_path <- args[1]
            tissue <- basename(substr(vector_path, 1, nchar(vector_path)-7))
            results_dir <- paste(dirname(dirname(vector_path)), 'results', sep='/')

            # Read in tables / patients
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            vector <- read.table(vector_path)$V1
            sub <- n[, colnames(n)%in%vector]
            setcolorder(sub, as.character(vector))

            # DESeq2 prep
            disease_vector <- rep(c('T', 'N'), length(vector)/2)
            countData <- round(sub)
            colData <- data.frame(disease=disease_vector, row.names=colnames(countData))
            y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ disease)

            # Run DESeq2
            y <- DESeq(y, parallel=TRUE)
            res <- results(y, parallel=TRUE)
            summary(res)

            # Write out table
            resOrdered <- res[order(res$padj),]
            res_name <- paste(tissue, 'results.tsv', sep='-')
            res_path <- paste(results_dir, res_name, sep='/')
            write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\\t',  quote=FALSE)
            """.format(cores=self.cores, df_path=self.df_path))

    def teardown(self):
        shutil.rmtree(self.results_dir)
