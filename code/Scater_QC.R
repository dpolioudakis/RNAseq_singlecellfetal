# Damon Polioudakis
# 2017-06-22
# Make Scater QC plots

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())

require(scater)
require(Seurat)
require(gridExtra)
require(ggplot2)

# Load
# Dropseq
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")

# Pollen
pnExDF <- read.csv("../pollen_2015/data/htseq/Exprs_HTSCexon.csv")
# # Picard Sequencing Statistics - bulk RNAseq
# picStatsBuDF <- read.table("../../kriegstein_2015/metadata/PicardToolsQC.csv")

# Fluidigm LT
flExDF <- read.csv("../C196-001_002_DP/data/htseq/merged/Exprs_HTSCexon.csv"
  , row.names = 1)

# Fluidigm HT
fhExDF <- read.csv(
  "../HT-003_DP/data/htseq/GRCh37.75_NoERCC/Exprs_HTSCexon_FtMm3e4.csv")

# biomart gene symbols and ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Variables
outGraph <- "../analysis/graphs/Scater_QC_"
################################################################################

### Functions

# Convert ensembl IDs to gene symbols for expression matrix row names
# Sum counts for each cell for duplicated symbols
# Keep ensembl ID when no associated gene symbol
Rows_Ensembl_To_Gene <- function(exDF, bmDF) {
  # Use gene symbol if available, otherwise use ensembl id
  idx <- match(row.names(exDF), bmDF$ensembl_gene_id)
  exDF <- exDF[! is.na(idx), ]
  idx <- idx[! is.na(idx)]
  gs <- as.character(bmDF$hgnc_symbol[idx])
  ei <- as.character(bmDF$ensembl_gene_id[idx])
  gs[gs == ""] <- ei[gs == ""]
  
  # List of ensembl IDs associated with duplicated symbols
  du <- gs[duplicated(gs)]
  l <- lapply(du, function(du) {ei[gs %in% du]})
  
  # Sum counts for each cell for duplicated symbols
  l <- lapply(l, function(eis) {colSums(exDF[row.names(exDF) %in% eis, ])})
  names(l) <- du
  
  # Add gene symbols
  exDF$GENE_SYMBOL <- gs
  nrow(exDF)
  # Remove duplicated gene symbols
  exDF <- exDF[! duplicated(exDF$GENE_SYMBOL), ]
  nrow(exDF)
  
  # Add summed counts for duplicated gene symbols
  for (i in 1:length(names(l))){
    print(i)
    name <- names(l)[i]
    print(name)
    exDF[exDF$GENE_SYMBOL == name, -which(names(exDF) == "GENE_SYMBOL")] <- l[[name]]
  }
  
  # Check
  print(sum(l[["PRICKLE4"]]))
  print(rowSums(exDF[exDF$GENE_SYMBOL == "PRICKLE4", -ncol(exDF)]))
  
  # Move to rownames and remove column
  row.names(exDF) <- exDF$GENE_SYMBOL
  exDF <- exDF[, -which(names(exDF) == "GENE_SYMBOL")]
  
  return(exDF)
}
################################################################################

### Format

# Drop-seq
exDF <- centSO@raw.data
# centSO@meta.data <- centSO@meta.data
centSO@meta.data$Region <- "GZ"
centSO@meta.data$Region[grep("CP", row.names(centSO@meta.data))] <- "CP"
centSO@meta.data$Cell <- row.names(centSO@meta.data)
exDF <- exDF[ ,colnames(exDF) %in% centSO@meta.data$Cell]

# Remove ERCCs and STAR stats from tail
# Fluidigm HT
tail(fhExDF, 10)
fhExDF <- head(fhExDF, -5)
tail(fhExDF, 5)
# Fluidigm LT
tail(flExDF, 10)
flExDF <- head(flExDF, -97)
tail(flExDF, 5)

# Move gene names
# Pollen
row.names(pnExDF) <- pnExDF$X
pnExDF <- pnExDF[ ,-1]
# Fluidigm HT
row.names(fhExDF) <- fhExDF$X
fhExDF <- fhExDF[ ,-1]

# Clean chr number of ens IDs for bulk, Pollen
# Pollen
row.names(pnExDF) <- gsub("\\.[0-9]*", "", row.names(pnExDF))
# Fluidigm HT
row.names(fhExDF) <- gsub("\\.[0-9]*", "", row.names(fhExDF))
# Fluidigm LT
row.names(flExDF) <- gsub("\\.[0-9]*", "", row.names(flExDF))

# Convert ensembl IDs to gene symbols for expression matrix row names
# Sum counts for each cell for duplicated symbols
# Keep ensembl ID when no associated gene symbol
fhExDF <- Rows_Ensembl_To_Gene(exDF = fhExDF, bmDF = bmDF)
flExDF <- Rows_Ensembl_To_Gene(exDF = flExDF, bmDF = bmDF)
pnExDF <- Rows_Ensembl_To_Gene(exDF = pnExDF, bmDF = bmDF)
################################################################################

### QC

# We use these objects to form an SCESet object containing all of the necessary
# information for our analysis:
# Drop-seq
pd <- new("AnnotatedDataFrame", data = centSO@meta.data)
rownames(pd) <- pd$Cell
dsSCE <- newSCESet(countData = exDF, phenoData = pd)
# Fluidigm HT
fhSCE <- newSCESet(countData = fhExDF)
# Fluidigm LT
flSCE <- newSCESet(countData = flExDF)
# Pollen
pnSCE <- newSCESet(countData = pnExDF)

# # Subsetting is very convenient with this class. For example, we can filter out
# # features (genes) that are not expressed in any cells:
# keep_feature <- rowSums(exprs(dsSCE) > 0) > 0
# dsSCE <- dsSCE[keep_feature,]

# Now we have the expression data neatly stored in a structure that can be used
# for lots of exciting analyses.
# It is straight-forward to compute many quality control metrics:
dsSCE <- calculateQCMetrics(dsSCE)
fhSCE <- calculateQCMetrics(fhSCE)
flSCE <- calculateQCMetrics(flSCE)
pnSCE <- calculateQCMetrics(pnSCE)

varLabels(dsSCE)

# # The first step in the QC process is filtering out unwanted features. We will
# # typically filter out features with very low overall expression, and any others
# # that plots or other metrics indicate may be problematic. First we look at a
# # plot that shows the top 50 (by default) most-expressed features. By default,
# # “expression” is defined using the feature counts (if available), but tpm, cpm,
# # fpkmor the exprs values can be used instead, if desired.
# keep_feature <- rowSums(counts(dsSCE) > 0) > 4
# dsSCE <- dsSCE[keep_feature,]

## Plot QC

# First we look at a plot that shows the top 50 (by default) most-expressed features.
p1 <- plotQC(dsSCE, type = "highest-expression", exprs_values = "counts")
p1 <- p1 + ggtitle(paste0(
  "Drop-seq"
  , "\n", p1$labels$title))
p2 <- plotQC(fhSCE, type = "highest-expression", exprs_values = "counts")
p2 <- p2 + ggtitle(paste0(
  "Fluidigm HT"
  , "\n", p2$labels$title))
p3 <- plotQC(flSCE, type = "highest-expression", exprs_values = "counts")
p3 <- p3 + ggtitle(paste0(
  "Fluidigm LT"
  , "\n", p3$labels$title))
p4 <- plotQC(pnSCE, type = "highest-expression", exprs_values = "counts")
p4 <- p4 + ggtitle(paste0(
  "Pollen"
  , "\n", p4$labels$title))
png(paste0(outGraph, "HighestExpression.png")
  , width = 12, height = 7, units = "in", res = 300)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "Top 50 most expressed features")
dev.off()

# Another way to obtain an idea of the level of technical noise in the dataset
# is to plot the frequency of expression (that is, number of cells with
# expression for the gene above the defined threshold (default is zero)) against
# mean expression expression level . A set of specific features to plot can be
# defined, but need not be. By default, the function will look for defined
# feature controls (as supplied to calculateQCMetrics). If feature controls are
# found, then these will be plotted, if not then all features will be plotted.
pdf(paste0(outGraph, "ExprsFreqVsMean.pdf"))
plotQC(dsSCE, type = "exprs-freq-vs-mean")
## `geom_smooth()` using method = 'loess'
dev.off()

# Beyond these QC plots, we have a neat, general and flexible function for
# plotting two feature metadata variables:
png(paste0(outGraph, "nCellsExprsVsPctTotalCounts.png")
  , width = 12, height = 6, units = "in", res = 300)
plotFeatureData(dsSCE, aes(x = n_cells_exprs, y = pct_total_counts))
dev.off()

# We can see that there is a small number of features that are ubiquitously
# expressed expressed in all cells (n_cells_exprs) and account for a large
# proportion of all counts observed (pct_total_counts; more than 0.5% of all
# counts). The subsetting of rows of SCESet objects makes it easy to drop
# unwanted features.

# We also have neat functions to plot two cell metadata variables:
pdf(paste0(outGraph, "TotalCountsVsTotalFeatures.pdf"))
plotPhenoData(dsSCE, aes(x = total_counts, y = total_features,
  colour = log10_total_counts))
dev.off()

# See the plotQC options below. The various plotting functions enable
# visualisation of the relationship betwen experimental variables and the
# expression data.We can look at the relative importance of different
# explanatory variables with some of the plotQC function options. We can compute
# the median marginal R2 for each variable in pData(dsSCE) when fitting
# a linear model regressing exprs values against just that variable.The default
# approach looks at all variables in pData(object) and plots the top
# nvars_to_plot variables (default is 10).
pdf(paste0(outGraph, "VarExplained.pdf"))
plotQC(dsSCE, type = "expl")
dev.off()


normalizeExprs(dsSCE, method = )
