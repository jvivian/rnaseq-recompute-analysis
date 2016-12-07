import os
import textwrap

from concurrent.futures import ThreadPoolExecutor

from experiments.AbstractExperiment import AbstractExperiment
import logging
from tqdm import tqdm
import pandas as pd

from utils import dedupe, run_deseq2
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class GtexOneVsAll(AbstractExperiment):

    def __init__(self, root_dir, cores):
        super(GtexOneVsAll, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(root_dir, 'experiments/gtex-matched-tissues')
        self.gtex_tissue = os.path.join(self.experiment_dir, 'gtex-tissue.tsv')
        self.tissue_dirs = [os.path.join(self.experiment_dir, x) for x in self.tissues]

    def setup(self):
        dirtree = [os.path.join(x, 'plots') for x in self.tissue_dirs]
        self.create_directories(dirtree)

        log.info('Reading in dataframes')
        dfs = []
        type_vector = []
        for df_path in tqdm(self.protein_coding_paths):
            df = pd.read_csv(df_path, sep='\t', index_col=0)
            type_vector.extend(os.path.basename(os.path.split(df_path)[0]) for _ in df.columns)
            dfs.append(df)

        log.info('Writing out tissue vector')
        for tissue in dedupe(type_vector):
            vector = [x if x == tissue else 'Not-' + tissue for x in type_vector]
            with open(os.path.join(self.experiment_dir, tissue, 'tissue-vector.txt'), 'w') as f:
                f.write('\n'.join(vector))

        log.info('Combining and outputting dataframe')
        if not os.path.exists(self.gtex_tissue):
            df = pd.concat(dfs, axis=1)
            df.to_csv(self.gtex_tissue, sep='\t')

        log.info('Writing out DESeq2.R scripts')
        [write_script(self.deseq2_script, x) for x in self.tissue_dirs]

    def run_experiment(self):
        df_path = [self.gtex_tissue for _ in self.tissues]
        script_path = [os.path.join(x, 'deseq2.R') for x in self.tissue_dirs]
        blob = zip(script_path, zip(df_path, self.tissue_dirs))

        with ThreadPoolExecutor(max_workers=self.cores) as executor:
            executor.map(run_deseq2, blob)

    def teardown(self):
        pass

    def deseq2_script(self):
        return textwrap.dedent("""
        library('DESeq2')

        # Argument parsing
        args <- commandArgs(trailingOnly = TRUE)
        df_path <- args[1]
        # type_path <- args[2]
        tissue_path <- args[2]

        script.dir <- dirname(tissue_path)

        # Read in tables / vectors
        # design = ~ type + condition
        # The specific tissue is the variable of interest
        # and we're controlling for all the other tissues
        n <- read.table(df_path, sep='\\t', header=1, row.names=1)
        # type <- read.table(type_path)
        tissue <- read.table(tissue_path)

        # Create matrix vectors
        # type_vector <- type$V1 # type
        tissue_vector <- tissue$V1 # condition

        # DESeq2 preprocessing
        # Rounding the countData since DESeQ2 only accepts integer counts
        countData <- round(n)
        colData <- data.frame(tissue=tissue_vector, row.names=colnames(countData))
        y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ tissue)

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
