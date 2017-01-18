
library('DESeq2')

# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
df_path <- args[1]
vector_path <- args[2]
tissue <- basename(substr(vector_path, 1, nchar(vector_path)-7))
results_dir <- paste(dirname(dirname(vector_path)), 'results', sep='/')
plot_dir <- paste(dirname(dirname(vector_path)), 'plots', tissue, sep='/')

# Read in tables / patients
n <- read.table(df_path, sep='\t', header=1, row.names=1)
vector <- read.table(vector_path)$V1
sub <- n[, vector]

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
write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\t',  quote=FALSE)

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
