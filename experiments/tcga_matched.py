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


class TcgaMatched(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(TcgaMatched, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/tcga-matched')
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
                samples = [x for x in f_in.readline().strip().split('\t')]
                barcodes = [x[:-3] for x in samples]
                matched = list(set(x for x in barcodes if x + '-11' in samples and x + '-01' in samples))
                if matched:
                    matched_vector = [f(x) for x in matched for f in (lambda f1: f1 + '-01', lambda f2: f2 + '-11')]
                    vector_path = os.path.join(self.vector_dir, tissue + '-vector')
                    with open(vector_path, 'w') as f:
                        f.write('\n'.join(matched_vector))
                    self.vectors.append(vector_path)
                else:
                    log.info('No matching TCGA samples found for tissue: ' + tissue)
                    shutil.rmtree(os.path.join(self.plots_dir, tissue))

    def run_experiment(self):
        vectors = []
        for df in self.protein_coding_paths:
            tissue_vector = os.path.join(self.vector_dir, os.path.basename(os.path.dirname(df)) + '-vector')
            if tissue_vector in self.vectors:
                vectors.append([df, tissue_vector])
        blob = zip([self.script_path for _ in xrange(len(vectors))], vectors)

        log.info('Starting DESeq2 Runs using {} cores'.format(self.cores))

        with ThreadPoolExecutor(max_workers=int(self.cores)) as executor:
            executor.map(run_deseq2, blob)

    def teardown(self):
        for df_path in [os.path.join(self.results_dir, x) for x in os.listdir(self.results_dir)]:
            df = add_gene_names(df_path, self.gene_map)
            df.to_csv(df_path, sep='\t')

    def deseq2_script(self):
        return textwrap.dedent("""
            library('DESeq2')

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- args[1]
            vector_path <- args[2]
            tissue <- basename(substr(vector_path, 1, nchar(vector_path)-7))
            results_dir <- paste(dirname(dirname(vector_path)), 'results', sep='/')
            plot_dir <- paste(dirname(dirname(vector_path)), 'plots', tissue, sep='/')

            # Read in tables / patients
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            vector <- read.table(vector_path)$V1
            sub <- n[, colnames(n)%in%vector]

            # Create matrix vectors
            disease_vector <- rep(c('T', 'N'), length(vector)/2)
            patient_vector <- gsub('.{3}$', '', vector) # Remove barcode from vector to get patient vector

            # DESeq2 preprocessing
            # Rounding the countData since DESeQ2 only accepts integer counts
            # The design matrix is conditioned on the two vectors: patient and condition
            countData <- round(sub)
            colData <- data.frame(disease=disease_vector, patient=patient_vector, row.names=colnames(countData))
            y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ patient + disease)

            # Run DESeq2
            y <- DESeq(y)
            res <- results(y)
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
            """)
