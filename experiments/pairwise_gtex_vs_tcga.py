import logging
import os
import shutil
import textwrap

from concurrent.futures import ThreadPoolExecutor

from experiments.AbstractExperiment import AbstractExperiment
from utils import run_deseq2, add_gene_names
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class PairwiseTcgaVsGtex(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(PairwiseTcgaVsGtex, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/pairwise_tcga_vs_gtex')
        self.tissue_dirs = [os.path.join(self.experiment_dir, x) for x in self.tissues]
        self.vector_dirs = [os.path.join(x, 'vectors') for x in self.tissue_dirs]
        self.script_path = None

    def setup(self):
        for subdir in ['vectors', 'results', 'plots', 'masks']:
            self.create_directories([os.path.join(x, subdir) for x in self.tissue_dirs])

        self.script_path = write_script(self.deseq2_script, directory=self.experiment_dir)

        log.info('Writing out vectors')
        for df in self.protein_coding_paths:
            tissue = os.path.basename(os.path.dirname(df))
            tissue_dir = os.path.join(self.experiment_dir, tissue)
            with open(df, 'r') as f:
                samples = [x for x in f.readline().strip().split('\t')]
                gtex = [x for x in samples if 'GTEX-' in x]
                tcga = [x for x in samples if 'TCGA-' in x and x.endswith('-01')]
                if gtex and tcga:
                    for sample in tcga:
                        vector = gtex + [sample]
                        with open(os.path.join(tissue_dir, 'vectors', sample + '-vector'), 'w') as f_out:
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

        log.info('Starting DESeq2 Runs using {} cores'.format(self.cores))
        with ThreadPoolExecutor(max_workers=self.cores) as executor:
            executor.map(run_deseq2, blob)

        self.create_masks()

        self.combine_results()

    def create_masks(self):
        pass

    def combine_results(self):
        pass

    def teardown(self):
        pass

    def deseq2_script(self):
        return textwrap.dedent("""
            suppressMessages(library('DESeq2'))

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- args[1]
            vector_path <- args[2]
            results_dir <- paste(dirname(dirname(vector_path)), 'results', sep='/')
            plot_dir <- paste(dirname(dirname(vector_path)), 'plots', tissue, sep='/')

            # Read in tables / patients
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            vector <- read.table(vector_path)$V1
            sub <- n[, vector]

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
            res_name <- paste(tissue, 'results.tsv', sep='-')
            res_path <- paste(results_dir, res_name, sep='/')
            write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\\t',  quote=FALSE)
            """)
