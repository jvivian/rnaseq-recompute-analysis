import logging
import os
import pickle
import textwrap
from collections import defaultdict

import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

from experiments.abstractExperiment import AbstractExperiment
from utils import run_deseq2
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class PairwiseGTEx(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(PairwiseGTEx, self).__init__(root_dir)
        self.cores = int(cores)
        self.experiment_dir = os.path.join(root_dir, 'experiments/pairwise-gtex')
        self.tissue_dirs = [os.path.join(self.experiment_dir, x) for x in self.tissues]
        self.output_df = os.path.join(self.experiment_dir, 'gtex-combined.tsv')
        self.script_path = None

    def setup(self):
        dirtree = [os.path.join(x, 'samples') for x in self.tissue_dirs] + \
                  [os.path.join(x, 'results') for x in self.tissue_dirs]
        self.create_directories(dirtree)

        self.script_path = write_script(self.deseq2_script, self.experiment_dir)

        log.info('Writing out combined GTEx dataframe')
        if not os.path.exists(self.output_df):
            dfs = [pd.read_csv(x, sep='\t', index_col=0).filter(regex=('GTEX')) for x in self.protein_coding_paths]
            df = pd.concat(dfs, axis=1)
            df.to_csv(self.output_df, sep='\t')
            all_samples_vector = set(df.columns)
        else:
            all_samples_vector = set(open(self.output_df, 'r').readline().strip().split('\t'))

        # A vector consists of all current tissue samples + one for each sample NOT in this set
        log.info('Writing out vectors')
        for df in tqdm(self.protein_coding_paths):
            tissue = os.path.basename(os.path.dirname(df))
            tissue_set = {x for x in open(df, 'r').readline().strip().split('\t') if 'GTEX' in x}
            for sample in (all_samples_vector - tissue_set):
                sample_vector = [x.replace('-', '.') for x in list(tissue_set) + [sample]]
                with open(os.path.join(self.experiment_dir, tissue, 'samples', sample), 'w') as f:
                    f.write('\n'.join(sample_vector))

    def run_experiment(self):
        vectors = [[os.path.join(x, 'samples', y)] for x in self.tissue_dirs
                   for y in os.listdir(os.path.join(x, 'samples'))]
        blob = zip([self.script_path for _ in xrange(len(vectors))], vectors)

        log.info('Starting DESeq2 Runs using {} cores'.format(self.cores))
        with ThreadPoolExecutor(max_workers=self.cores) as executor:
            executor.map(run_deseq2, blob)

        self.reduce()

    def reduce(self):
        """Reduce results for each tissue into a single dataframe with p-value counts"""
        log.info('Reducing results for each tissue into a single dataframe of p-value counts')
        gene_map = pickle.load(open(self.gene_map, 'rb'))
        for tissue_dir in tqdm(self.tissue_dirs):
            results_path = os.path.join(tissue_dir, 'results.tsv')
            ranked = pd.DataFrame()
            pvals = defaultdict(list)
            fc = defaultdict(list)
            results = [os.path.join(tissue_dir, 'results', x) for x in os.listdir(os.path.join(tissue_dir, 'results'))]
            for result in tqdm(results):
                df = pd.read_csv(result, index_col=0, sep='\t')
                for gene in df.index:
                    pvals[gene].append(df.loc[gene]['padj'])
                    fc[gene].append(df.loc[gene]['log2FoldChange'])

            genes = pvals.keys()
            ranked['num_samples'] = [len(pvals[x]) for x in genes]
            ranked['pval_counts'] = [sum([1 for y in pvals[x] if y < 0.001]) for x in genes]
            ranked['pval'] = [np.median(pvals[x]) for x in genes]
            ranked['pval_std'] = [round(np.std(pvals[x]), 4) for x in genes]
            ranked['fc'] = [round(np.median(fc[x]), 4) for x in genes]
            ranked['fc_std'] = [round(np.std(fc[x]), 4) for x in genes]

            gene_names = [gene_map[x] if x in gene_map.keys() else x for x in genes]
            ranked['gene_id'] = genes
            ranked.index = gene_names
            ranked.sort_values('pval_counts', inplace=True, ascending=False)

            ranked.to_csv(results_path)

    def teardown(self):
        pass

    def deseq2_script(self):
        return textwrap.dedent("""
        suppressMessages(library('DESeq2'))
        suppressMessages(library('data.table'))

        args <- commandArgs(trailingOnly = TRUE)
        df_path <- '{df_path}'
        vector_path <- args[1]
        vector_name <- basename(vector_path)  # Name of vector should be the GTEx sample
        vector_dir <- dirname(vector_path)
        results_dir <- paste(dirname(vector_dir), 'results', sep='/')

        n <- read.table(df_path, sep='\\t', header=1, row.names=1)
        vector <- read.table(vector_path)$V1

        # Subset dataframe
        sub <- n[, colnames(n)%in%vector]
        setcolorder(sub, as.character(vector))

        # Create condition vector
        condition <- c(rep('A', length(sub)-1), 'B')

        # DESeq2 Preprocessing
        countData <- round(sub)
        colData <- data.frame(condition=condition, row.names=colnames(countData))
        y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

        # Run DESeq2
        y <- DESeq(y)
        res <- results(y)

        # Write out table
        resOrdered <- res[order(res$padj),]
        results_name <- paste(results_dir, vector_name, sep='/')
        write.table(as.data.frame(resOrdered), file=results_name, col.names=NA, sep='\\t',  quote=FALSE)
        """.format(df_path=self.output_df))
