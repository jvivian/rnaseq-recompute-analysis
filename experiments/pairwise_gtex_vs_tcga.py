import logging
import os
import pickle
import textwrap
from collections import defaultdict
from collections import Counter

import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

from experiments.AbstractExperiment import AbstractExperiment
from utils import run_deseq2
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class PairwiseTcgaVsGtex(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(PairwiseTcgaVsGtex, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/pairwise-tcga-vs-gtex')
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
                gtex = [x.replace('-', '.') for x in samples if 'GTEX-' in x]
                tcga = [x.replace('-', '.') for x in samples if 'TCGA-' in x]
                if gtex and tcga:
                    for sample in tcga:
                        vector = gtex + [sample]
                        with open(os.path.join(tissue_dir, 'vectors', sample), 'w') as f_out:
                            f_out.write('\n'.join(vector))

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

        self.create_masks()

        log.info('Reducing results for each tissue into a single dataframe sorted by p-value counts')
        self.combine_results()
        log.info('Reducing results for normal tissues sorted by p-value counts')
        self.combine_results(output_name='normal-results.tsv', normal=True)
        log.info('Reducing masked reuslts for each tissue into a single dataframe sorted by p-value counts')
        self.combine_results(output_name='results-masked.tsv', result_dir='masked-results')

    def create_masks(self):
        log.info('Creating masks for matched samples')
        for result_dir in tqdm(self.results_dirs):
            masked_result_dir = os.path.join(os.path.dirname(result_dir), 'masked-results')
            masked_gene_dir = os.path.join(os.path.dirname(result_dir), 'masked-genes')
            results = [os.path.join(result_dir, x) for x in os.listdir(result_dir)]
            result_names = [os.path.basename(x) for x in results]

            # Trim last 3 characters to get TCGA barcode without the tumor / normal tag
            barcodes = [x[:-3]for x in result_names]
            matches = {item for item, count in Counter(barcodes).items() if count == 2}
            matches = {x for x in matches if x + '.01' in result_names and x + '.11' in result_names}

            # For each match, produce a masked result file in the masked-results directory
            for match in matches:
                df_norm = pd.read_csv(os.path.join(result_dir, match + '.11'), sep='\t', index_col=0)
                df_tumor = pd.read_csv(os.path.join(result_dir, match + '.01'), sep='\t', index_col=0)
                masked_genes = df_norm[df_norm.padj < 0.001].index
                for gene in masked_genes:
                    try:
                        df_tumor.drop(gene, inplace=True)
                    except ValueError:  # DE gene in normal didn't show up in tumor
                        pass
                df_tumor.to_csv(os.path.join(masked_result_dir, match + '.01'), sep='\t')

                # Convert genes from ensembl to gene names and add to masked-genes directory
                gene_map = pickle.load(open(self.gene_map, 'rb'))
                gene_names = [gene_map[x] if x in gene_map.keys() else x for x in masked_genes]

                with open(os.path.join(masked_gene_dir, match + '.01'), 'w') as f:
                    f.write('\n'.join(gene_names))

    def combine_results(self, output_name='results.tsv', result_dir='results', normal=True):
        gene_map = pickle.load(open(self.gene_map, 'rb'))
        for tissue_dir in sorted(self.tissue_dirs):
            log.info('Processing ' + os.path.basename(tissue_dir))
            results_path = os.path.join(tissue_dir, output_name)
            ranked = pd.DataFrame()
            pvals = defaultdict(list)
            fc = defaultdict(list)
            sample_suffix = '.11' if normal else '.01'
            results = [os.path.join(tissue_dir, 'results', x) for x in
                       os.listdir(os.path.join(tissue_dir, result_dir)) if x.endswith(sample_suffix)]
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
            disease_vector <- c(rep('G', length(vector)-1), 'T')

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
