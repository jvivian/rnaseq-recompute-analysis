import logging
import os
import textwrap
from subprocess import check_call

from experiments.AbstractExperiment import AbstractExperiment
from utils import add_gene_names
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class TcgaTumorVsNormal(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(TcgaTumorVsNormal, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/tcga-tumor-vs-normal')
        self.vector_dir = os.path.join(self.experiment_dir, 'vectors')
        self.results_dir = os.path.join(self.experiment_dir, 'results')
        self.plots_dir = os.path.join(self.experiment_dir, 'plots')
        self.script_path = None
        self.vectors = []

    def setup(self):
        dirtree = [self.vector_dir, self.results_dir] + [os.path.join(self.plots_dir, x) for x in self.tissues]
        self.create_directories(dirtree)

        self.script_path = write_script(self.deseq2_script, directory=self.experiment_dir)

        log.info('Writing out tissue vectors')
        for df in self.protein_coding_paths:
            tissue = os.path.basename(os.path.dirname(df))
            with open(df, 'r') as f_in:
                samples = [x for x in f_in.readline().strip().split('\t') if x.startswith('TCGA')]
                samples = [x.replace('-', '.') for x in samples if x.endswith('11') or x.endswith('01')]
                if any(x for x in samples if x.endswith('11')):
                    vector_path = os.path.join(self.vector_dir, tissue + '-vector')
                    with open(vector_path, 'w') as f:
                        f.write('\n'.join(samples))

                    disease_vector = ['T' if x.endswith('01') else 'N' for x in samples]
                    disease_vector_path = os.path.join(self.vector_dir, tissue + '-disease')
                    with open(disease_vector_path, 'w') as f:
                        f.write('\n'.join(disease_vector))

    def run_experiment(self):
        log.info("what")
        for df in sorted(self.protein_coding_paths):
            tissue = os.path.basename(os.path.dirname(df))
            tissue_vector = tissue + '-vector'
            disease_vector = tissue_vector.replace('-vector', '-disease')
            if os.path.exists(os.path.join(self.vector_dir, tissue_vector)):
                cmd = ['docker', 'run',
                       '--rm',
                       '--log-driver=none',
                       '-v', self.results_dir + ':/results',
                       '-v', self.vector_dir + ':/vectors',
                       '-v', self.plots_dir + '/' + tissue + ':/plots',
                       '-v', os.path.dirname(df) + ':/df',
                       '-v', self.experiment_dir + ':/data',
                       'jvivian/deseq2:1.14.1', '/data/deseq2.R',
                       '/df/' + os.path.basename(df),
                       '/vectors/' + tissue_vector,
                       '/vectors/' + disease_vector]
                log.info('Starting DESeq2 Run for {} using {} cores'.format(tissue, self.cores))
                check_call(cmd)

    def teardown(self):
        log.info('Adding gene names to results.')
        for df_path in [os.path.join(self.results_dir, x) for x in os.listdir(self.results_dir)]:
            df = add_gene_names(df_path, self.gene_map)
            df.to_csv(df_path, sep='\t')

    def deseq2_script(self):
        return textwrap.dedent("""
            suppressMessages(library('DESeq2'))
            suppressMessages(library('data.table'))
            suppressMessages(library('BiocParallel'))
            register(MulticoreParam({cores}))

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- args[1]
            vector_path <- args[2]
            disease_path <- args[3]
            tissue <- basename(substr(vector_path, 1, nchar(vector_path)-7))
            results_dir <- '/results'
            plot_dir <- '/plots'

            # Read in tables / patients
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            vector <- read.table(vector_path)$V1
            disease_vector <- read.table(disease_path)$V1
            sub <- n[, colnames(n)%in%vector]
            setcolorder(sub, as.character(vector))

            # DESeq2 preprocessing
            # Rounding the countData since DESeQ2 only accepts integer counts
            # The design matrix is conditioned on the two vectors: patient and condition
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

            # MA Plot
            ma_name <- paste(plot_dir, 'MA.pdf', sep='/')
            pdf(ma_name, width=7, height=7)
            plotMA(res, main='DESeq2')
            dev.off()

            # Dispersion Plot
            disp_name <- paste(plot_dir, 'dispersion.pdf', sep='/')
            pdf(disp_name, width=7, height=7)
            plotDispEsts( y, ylim = c(1e-6, 1e1) )
            dev.off()

            # PVal Hist
            hist_name <- paste(plot_dir, 'pval-hist.pdf', sep='/')
            pdf(hist_name, width=7, height=7)
            hist( res$pvalue, breaks=20, col="grey" )
            dev.off()

            # Ratios plots
            qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
            bins <- cut( res$baseMean, qs )
            levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
            ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
            ratio_name <- paste(plot_dir, 'ratios.pdf', sep='/')
            pdf(ratio_name, width=7, height=7)
            barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
            dev.off()
            """.format(cores=self.cores))
