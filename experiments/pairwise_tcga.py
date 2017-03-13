import logging
import os
import pickle
import textwrap
from collections import defaultdict

import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

from experiments.AbstractExperiment import AbstractExperiment
from utils import run_deseq2
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class PairwiseTCGA(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(PairwiseTCGA, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/pairwise-tcga')
        self.tissue_dirs = [os.path.join(self.experiment_dir, x) for x in self.tissues]
        self.vector_dirs = [os.path.join(x, 'vectors') for x in self.tissue_dirs]
        self.results_dirs = [os.path.join(x, 'results') for x in self.tissue_dirs]
        self.script_path = None

    def setup(self):
        for subdir in ['vectors', 'results', 'masked-results', 'masked-genes', 'masks']:
            self.create_directories([os.path.join(x, subdir) for x in self.tissue_dirs])

        self.script_path = write_script(self.deseq2_script, directory=self.experiment_dir)

        log.info('Writing out vectors')
        for df in tqdm(self.protein_coding_paths):
            tissue = os.path.basename(os.path.dirname(df))
            tissue_dir = os.path.join(self.experiment_dir, tissue)
            with open(df, 'r') as f:
                samples = [x for x in f.readline().strip().split('\t')]
                tcga = [x.replace('-', '.') for x in samples if 'TCGA-' in x]
                tcga_t = [x for x in tcga if x.endswith('01')]
                tcga_n = [x for x in tcga if x.endswith('11')]
                if tcga_t and tcga_n:
                    for sample in tcga_t:
                        vector = tcga_n + [sample]
                        with open(os.path.join(tissue_dir, 'vectors', sample), 'w') as f_out:
                            f_out.write('\n'.join(vector))
                else:
                    self.results_dirs.remove(os.path.join(self.experiment_dir, tissue, 'results'))

    def run_experiment(self):
        df_vector_pairs = []
        for df in self.protein_coding_paths:
            tissue = os.path.basename(os.path.dirname(df))
            vector_dir = os.path.join(self.experiment_dir, tissue, 'vectors')
            vectors = [os.path.join(vector_dir, x) for x in os.listdir(vector_dir)]
            for vector in vectors:
                df_vector_pairs.append([df, vector])

        blob = zip([self.script_path for _ in xrange(len(df_vector_pairs))], df_vector_pairs)

        if not os.listdir(self.results_dirs[0]):
            log.info('Starting DESeq2 Runs using {} cores'.format(self.cores))
            with ThreadPoolExecutor(max_workers=self.cores) as executor:
                executor.map(run_deseq2, blob)
        else:
            log.info('Results detected, skipping DESeq2 Run')

        log.info('Reducing results for each tissue into a single dataframe sorted by p-value counts')
        self.combine_results()

    def combine_results(self):
        gene_map = pickle.load(open(self.gene_map, 'rb'))
        for tissue_dir in sorted(self.tissue_dirs):
            log.info('Processing ' + os.path.basename(tissue_dir))
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

            ranked.to_csv(results_path, sep='\t')

    def teardown(self):
        pass

    def deseq2_script(self):
        return textwrap.dedent("""
            suppressMessages(library('DESeq2'))
            suppressMessages(library('data.table'))

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- args[1]
            vector_path <- args[2]
            results_dir <- paste(dirname(dirname(vector_path)), 'results', sep='/')
            sample_name <- basename(vector_path)

            # Read in tables / patients
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            vector <- read.table(vector_path)$V1
            sub <- n[, colnames(n)%in%vector]
            setcolorder(sub, as.character(vector))

            # Create matrix vectors
            disease_vector <- c(rep('N', length(vector)-1), 'T')

            # DESeq2 preprocessing
            # Rounding the countData since DESeQ2 only accepts integer counts
            countData <- round(sub)
            colData <- data.frame(disease=disease_vector, row.names=colnames(countData))
            y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ disease)

            # Run DESeq2
            y <- DESeq(y)
            res <- results(y)
            summary(res)

            # Write out table
            resOrdered <- res[order(res$padj),]
            res_path <- paste(results_dir, sample_name, sep='/')
            write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\\t',  quote=FALSE)
            """)
