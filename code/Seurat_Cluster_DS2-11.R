# Damon Polioudakis
# 2017-05-28
# Clustering of Drop-seq cells by digital gene expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
require(irlba)
require(gridExtra)
require(cowplot)
require("ggdendro")
source("Function_Library.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)

# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

## Input data
# Digital gene expression

# DS-002-003-004-005-006-007-008-009-011: exDF and metDF
load("../analysis/analyzed_data/Expression_Matrix_Compile/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7-8-9-11.Rdata")
# load("../analysis/analyzed_data/Expression_Matrix_Compile_dge_FtMm250_Magic_DS-2-3-4-5-6-7-8-9-11_exon_FtMm250_200-3sdgd.Rdata")
# exDF <- as.data.frame(exDF)
# # Subsetting for testing
# idx <- sample(1:ncol(exDF), 2000)
# exDF <- exDF[ ,idx]
# metDF <- metDF[idx, ]

# BiomaRt ensembl ID, gene symbol table
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Seurat clustering of DS-002-003 to compare cluster identities
load("../analysis/analyzed_data/Cluster_Seurat/DS2-3/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

# Cell cycle markers from Macosko 2015 Table S2 to remove from variable gene
# list used for clustering
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)
cc.genes <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]

## Variables
graphCodeTitle <- "Seurat_Cluster_DS2-11.R"
# outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_"
# outTable <- "../analysis/tables/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_"
# outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40_"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC/Seurat_Cluster_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC/Seurat_Cluster_DS2-11_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_"
# outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC/Seurat_Cluster_DS2-11_"
# outTable <- "../analysis/tables/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC/Seurat_Cluster_DS2-11_"
# outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC_PC1to40/Seurat_Cluster_DS2-11_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outData)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Convert Ensembl ID to Gene symbol

print("### Convert Ensembl ID to Gene symbol")

# In case of multiple ensembl IDs mapping to the same gene symbol, counts are
# summed across each ensembl ID associated with that gene symbol. If there is no
# gene symbol then ensembl ID is kept

# Clean chr number of ens IDs
row.names(exDF) <- gsub("\\.[0-9]*", "", row.names(exDF))

# Use gene symbol if available, otherwise use ensembl id
idx <- match(row.names(exDF), bmDF$ensembl_gene_id)
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
sum(l[["PRICKLE4"]])
rowSums(exDF[exDF$GENE_SYMBOL == "PRICKLE4", -ncol(exDF)])

# Move to rownames and remove column
row.names(exDF) <- exDF$GENE_SYMBOL
exDF <- exDF[, -which(names(exDF) == "GENE_SYMBOL")]
################################################################################

### Filtering cells for further analysis

print("### Filtering cells for further analysis")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(exDF))
dense.size
sparse.size <- object.size(exDF)
sparse.size
dense.size/sparse.size

# Number genes detected
print("Mean number of genes detected with no cells filtered:")
print(mean(colSums(exDF >= 1)))
# Number cell input
print("Number of cells input:")
print(ncol(exDF))
filtersDF <- data.frame(FILTER = "None", REMAINING = ncol(exDF)
  , MEAN_GENES_DETECTED = mean(colSums(exDF >= 1)))

# Initialize the Seurat object with the raw (non-normalized data)
# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules
# (as in Macosko et al. Cell 2015)
centSO <- CreateSeuratObject(raw.data = exDF, min.cells = 3, min.genes = 200
  , normalization.method = "LogNormalize", scale.factor = 10000
  , project = "DS-2-3-4-5-6-7-8-9-11", do.scale = FALSE, do.center = FALSE)

# Add metadata
metDF$CELL <- colnames(exDF)
row.names(metDF) <- metDF$CELL
centSO <- AddMetaData(centSO, metadata = metDF)

print("Number of cells remaining after filtering cells <200 genes detected:")
# Number genes detected
print("Mean number of genes detected after filtering cells <200 genes detected:")
print(mean(colSums(as.matrix(centSO@data) > 0)))
print(ncol(centSO@data))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), "< 200 genes detected")
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED
    , mean(colSums(as.matrix(centSO@data) > 0))))

# While the setup function imposes a basic minimum gene-cutoff, you may want to
# filter out cells at this stage based on technical or biological parameters.
# Seurat allows you to easily explore QC metrics and filter cells based on any
# user-defined criteria. In the example below,  we visualize gene and molecule
# counts, plot their relationship, and exclude cells with a clear outlier number
# of genes detected as potential multiplets. Of course this is not a guarenteed
# method to exclude cell doublets, but we include this as an example of
# filtering user-defined outlier cells. We also filter cells based on the
# percentage of mitochondrial genes present.

# nGene and nUMI are automatically calculated for every object by Seurat. For
# non-UMI data, nUMI represents the sum of the non-normalized values within a
# cell

# We calculate the percentage of mitochondrial genes here and store it in
# percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a
# common scRNA-seq QC metric.
mito.genes <- grep("^MT-", rownames(centSO@data), value = TRUE)
percent.mito <- apply(
  expm1(centSO@data[mito.genes, ]), 2, sum)/apply(expm1(centSO@data), 2, sum)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC
# stats
centSO <- AddMetaData(centSO, percent.mito, "percent.mito")

# We filter out cells that have unique gene counts over 2,500
# Note that high.thresholds and low.thresholds can be used to define a 'gate', and can
# filter cells not only based on nGene but on anything in the object (as in
# GenePlot above)
# centSO <- SubsetData(centSO, subset.name = "nGene", high.thresholds = 2500)
# Filter cells with gene counts > 3 sd above the mean
high.thresholds <- round(
  mean(centSO@meta.data$nGene) + 3*sd(centSO@meta.data$nGene), 0)
# 3192
centSO <- FilterCells(centSO, subset.name = "nGene"
  , high.thresholds = high.thresholds)
print(paste0("Number of cells remaining after filtering cells with >"
  , high.thresholds, " (3 sd above mean) genes detected:"))
print(ncol(centSO@data))
print(paste0("Mean number of genes detected after after filtering cells with >"
  , high.thresholds, " (3 sd above mean) genes detected:"))
print(mean(colSums(as.matrix(centSO@data) > 0)))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), paste0("> ", high.thresholds, " genes detected"))
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED, mean(colSums(as.matrix(centSO@data) > 0))))

centSO <- FilterCells(centSO, subset.name = "percent.mito"
  , high.thresholds = 0.05)
print("Number of cells remaining after filtering cells with >5% of reads mapping to exons mapping to mitochondrial exons:")
print(ncol(centSO@data))
print("Number of genes detected after filtering cells with >5% of reads mapping to exons mapping to mitochondrial exons:")
print(mean(colSums(as.matrix(centSO@data) > 0)))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), ">5% Mt")
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED
    , mean(colSums(as.matrix(centSO@data) > 0))))

write.csv(filtersDF, file = paste0(outTable, "Filters.csv"), quote = FALSE)
################################################################################

### Assign Cell-Cycle scores

print("### Assign Cell-Cycle scores")

# First, we assign each cell a score, based on its expression of G2/M and S
# phase markers. These marker sets should be anticorrelated in their expression
# levels, and cells expressing neither are likely not cycling and in G1 phase.

# We assign scores in the CellCycleScoring function, which stores S and G2/M
# scores in object@meta.data, along with the predicted classification of each
# cell in either G2M, S or G1 phase. CellCycleScoring can also set the identity
# of the Seurat object to the cell-cycle phase by passing set.ident = TRUE (the
# original identities are stored as old.ident). Please note that Seurat does not
# use the discrete classifications (G2M/G1/S) in downstream cell cycle
# regression. Instead, it uses the quantitative scores for G2M and S phase.
# However, we provide our predicted classifications in case they are of
# interest.

# Uses AddModuleScore()
# Calculate the average expression levels of each program (cluster)
# on single cell level, subtracted by the aggregated expression of
# control gene sets. All analyzed genes are binned based on averaged
# expression, and the control genes are randomly selected from each
# bin.

centSO <- CellCycleScoring(object = centSO, s.genes = s.genes
  , g2m.genes = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = centSO@meta.data)

tb <- table(centSO@meta.data$Phase)
write.csv(tb, file = paste0(outTable, "CellCycle.csv"), quote = FALSE)

# # Visualize the distribution of cell cycle markers across
# JoyPlot(object = centSO, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"),
#   nCol = 2, do.return = TRUE)
# ggsave(paste0(outGraph, "QC_CellCycle_DensityPlot.png"), width = 8, height = 4)

# # Run PCA on cell cycle genes
# tmp <- RunPCA(object = centSO, pc.genes = c(s.genes, g2m.genes)
#   , do.print = FALSE, pcs.compute = 50, pcs.print = 1:5, genes.print = 5
#   , maxit = 500, weight.by.var = FALSE)
# PCAPlot(object = tmp)
# ggsave(paste0(outGraph, "QC_CellCycle_PCA.png"))
# rm(tmp)
################################################################################

### Regress out unwanted sources of variation

print("### Regress out unwanted sources of variation")

# Your single cell dataset likely contains ‘uninteresting’ sources of variation.
# This could include not only technical noise, but batch effects, or even
# biological sources of variation (cell cycle stage). As suggested in Buettner
# et al, NBT, 2015, regressing these signals out of the analysis can improve
# downstream dimensionality reduction and clustering. To mitigate the effect of
# these signals, Seurat constructs linear models to predict gene expression
# based on user-defined variables. The scaled z-scored residuals of these models
# are stored in the scale.data slot, and are used for dimensionality reduction
# and clustering.

# note that this overwrites centSO@scale.data. Therefore, if you intend to use
# ScaleData, you can set do.scale=F and do.center=F in the original object to
# save some time.

# Individual for regressing out
individual <- rep("Brain 2", length(centSO@data@Dimnames[[2]]))
# Brain 3
br3Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "3"])
individual[centSO@data@Dimnames[[2]] %in% br3Cells] <- "Brain 3"
# Brain 4
br4Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "4"])
individual[centSO@data@Dimnames[[2]] %in% br4Cells] <- "Brain 4"
# Brain 5
br5Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "5"])
individual[centSO@data@Dimnames[[2]] %in% br5Cells] <- "Brain 5"
# Add to Seurat object
names(individual) <- centSO@data@Dimnames[[2]]
centSO <- AddMetaData(centSO, individual, "individual")

# Lab library prep for regressing out
librarylab <- rep("Plath", length(centSO@data@Dimnames[[2]]))
dhgCells <- c(as.character(metDF$CELL)[metDF$LIBRARY == "Geschwind"])
librarylab[centSO@data@Dimnames[[2]] %in% dhgCells] <- "Geschwind"
names(librarylab) <- centSO@data@Dimnames[[2]]
centSO <- AddMetaData(centSO, librarylab, "librarylab")

# Center and scale
centSO <- ScaleData(centSO
  # , vars.to.regress = c("nUMI", "percent.mito"))
  # , vars.to.regress = c("nUMI", "percent.mito", "librarylab", "individual"))
  , vars.to.regress = c("nUMI", "librarylab", "individual"))
  # , vars.to.regress = c("nUMI", "librarylab", "individual", "S.Score", "G2M.Score"))
# , vars.to.regress = c("nUMI"))

# No center or scale

# Make dataframe of covariates to regress out
covDF <- data.frame(nUMI = centSO@meta.data$nUMI
  , librarylab = centSO@meta.data$librarylab
  , individual = centSO@meta.data$individual)
# covDF <- data.frame(nUMI = centSO@meta.data$nUMI
#   , librarylab = centSO@meta.data$librarylab
#   , individual = centSO@meta.data$individual
#   , S.Score = centSO@meta.data$S.Score
#   , G2M.Score = centSO@meta.data$G2M.Score)
# covDF <- data.frame(nUMI = centSO@meta.data$nUMI)
# Regress out confounding variables
RegressCovariates <- function (exM, covDF) {
  exRegCovM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
  rownames(exRegCovM) <- rownames(exM)
  colnames(exRegCovM) <- colnames(exM)
  # ncol(covDF)+1 when condition has 2 levels
  coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
  for (i in 1:nrow(exM)) {
    if (i%%1000 == 0) {print(paste("Done for", i, "genes..."))}
    mod <- lm(as.numeric(exM[i, ]) ~ ., data = covDF)
    # The full data - the undesired covariates
    exRegCovM[i,] <- coef(mod)[1] + mod$residuals
    # lmmod1 <- lm(as.numeric(exM[i, ]) ~ condition + age + sex + pmi, data = covDF)
  }
  return(exRegCovM)
}

noCentExM <- RegressCovariates(centSO@data, covDF)

# noCentSO <- centSO
# noCentSO <- ScaleData(centSO
#   # , vars.to.regress = c("nUMI", "percent.mito", "librarylab", "individual")
#   # , vars.to.regress = c("nUMI", "percent.mito")
#   , vars.to.regress = c("nUMI", "librarylab", "individual")
#   # , vars.to.regress = c("nUMI")
#   , do.scale = FALSE, do.center = FALSE)
################################################################################

### QC plots

print("### QC plots")

# Paper - violin plot of nGene, nUMI, percent MT
fmetDF <- centSO@meta.data
ggDF <- fmetDF[ ,c("nGene", "nUMI", "percent.mito")]
ggDF$percent.mito <- ggDF$percent.mito * 100
ggDF <- melt(ggDF)
means_DF <- aggregate(value~variable, ggDF, mean)
# ggplot
ggplot(ggDF, aes(y = value, x = variable)) +
  geom_violin(aes(fill = variable)) +
  geom_text(data = means_DF, aes(label = round(value,2), y = value)) +
  facet_wrap(~variable, ncol = 3, scales = "free") +
  ggplot_set_theme_publication
ggsave(paste0(outGraph, "QC_violinPlot_paper.pdf"), width = 8, height = 4)

# Mean nUMI and mean nGene for filtered data
df1 <- data.frame(Mean_nGene = mean(centSO@meta.data$nGene)
  , Mean_nUMI = mean(centSO@meta.data$nUMI)
)
write.csv(df1, file = paste0(outTable, "mean_nUMI_nGene.csv"))

# Calculate metadata for all cells (no filters)
metDF$nGene <- colSums(exDF > 0)
metDF$nUMI <- colSums(exDF)
mito.genes <- grep("^MT-", rownames(exDF), value = TRUE)
percent.mito <- apply(
  expm1(exDF[mito.genes, ]), 2, sum)/apply(expm1(exDF), 2, sum)
metDF$percent.mito <- percent.mito

# Violin plot of nGene, nUMI, percent MT
metDF$FILTERED <- FALSE
fmetDF <- centSO@meta.data
fmetDF$FILTERED <- TRUE
ggDF <- rbind(metDF[ ,c("nGene", "nUMI", "percent.mito", "FILTERED")]
  , fmetDF[ ,c("nGene", "nUMI", "percent.mito", "FILTERED")])
ggDF <- melt(ggDF)
# ggplot
ggplot(ggDF, aes(x = FILTERED, y = value)) +
  geom_violin(aes(fill = FILTERED)) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  facet_wrap(~variable, ncol = 3, scales = "free_y")
ggsave(paste0(outGraph, "QC_violinPlot.png"), width = 8, height = 4)
# png(paste0(outGraph, "QC_violinPlot.png"), width = 10, height = 5
#   , units = "in", res = 300)
# VlnPlot(centSO, c("nGene", "nUMI", "percent.mito"), ident.include = NULL
#   ,nCol = 3, point.size.use = 0.1, group.by = "orig.ident")
# dev.off()
ggDF <- fmetDF[ ,c("nGene", "nUMI", "percent.mito", "FILTERED")]
ggDF <- ggDF[ggDF$FILTERED == TRUE, ]
ggDF <- melt(ggDF)
# ggplot
ggplot(ggDF, aes(x = FILTERED, y = value)) +
geom_violin(aes(fill = variable)) +
# geom_jitter(size = 0.1, alpha = 0.1) +
facet_wrap(~variable, ncol = 3, scales = "free_y")
ggsave(paste0(outGraph, "QC_violinPlot2.png"), width = 6, height = 6)

# Variance and per gene histogram
# Raw counts
v1 <- apply(exDF, 1, var)
ggDF <- data.frame(v1)
p1 <- ggplot(ggDF, aes(x = v1)) +
  geom_histogram(binwdith = 0.01) +
  xlim(c(-0.1,2)) +
  coord_cartesian(ylim = c(0, 4e4)) +
  xlab("Variance") +
  ggtitle(paste0("Raw counts"
    , "\n"))
# Normalized
v2 <- apply(noCentExM, 1, var)
gg2DF <- data.frame(v2)
p2 <- ggplot(gg2DF, aes(x = v2)) +
  geom_histogram(binwdith = 0.01) +
  xlim(c(-0.1,2)) +
  coord_cartesian(ylim = c(0, 4e4)) +
  xlab("Variance") +
  ggtitle(paste0("Normalized expression"
    , "\nFiltered cells and genes"))
png(paste0(outGraph, "QC_VarianceHistogram.png"))
grid.arrange(p1, p2, ncol = 2, top = paste0(graphCodeTitle
  , "\n"
  , "\nHistogram of variance per gene"))
dev.off()

# CoV and per gene histogram
# Calculate coefficient of variation
Coefficent_Of_Variation <- function(mean, sd){
  (sd/mean)*100
}
# Raw counts
# Calculate mean expression for each gene
mnL <- rowMeans(exDF)
# Standard deviation for each gene
sdL <- apply(exDF, 1, sd)
# Calculate CoV for each gene
# CoV for each gene
Coefficent_Of_Variation(mnL[[1]], sdL[[1]])
cvL <- mapply(Coefficent_Of_Variation, mnL, sdL)
gg1DF <- data.frame(cvL)
# ggplot
p1 <- ggplot(gg1DF, aes(x = cvL)) +
  geom_histogram(binwidth = 100) +
  # xlim(c(-0.1,2)) +
  # coord_cartesian(ylim = c(0, 4e4)) +
  xlab("Coefficient of variation") +
  ggtitle(paste0("Raw counts"
    , "\n"))
# Normalized
# Calculate mean expression for each gene
mnL <- rowMeans(noCentExM)
# Standard deviation for each gene
sdL <- apply(noCentExM, 1, sd)
# Calculate CoV for each gene
# CoV for each gene
Coefficent_Of_Variation(mnL[[1]], sdL[[1]])
cvL <- mapply(Coefficent_Of_Variation, mnL, sdL)
gg2DF <- data.frame(cvL)
# ggplot
p2 <- ggplot(gg2DF, aes(x = v2)) +
  geom_histogram(binwidth = 0.1) +
  # xlim(c(-0.1,2)) +
  # coord_cartesian(ylim = c(0, 4e4)) +
  xlab("Coefficient of variation") +
  ggtitle(paste0("Normalized expression"
    , "\nFiltered cells and genes"))
png(paste0(outGraph, "QC_CovHistogram.png"))
grid.arrange(p1, p2, ncol = 2, top = paste0(graphCodeTitle
  , "\n"
  , "\nHistogram of coefficient of variation per gene"))
dev.off()

# GenePlot is typically used to visualize gene-gene relationships, but can be
# used for anything calculated by the object, i.e. columns in object@meta.data,
# PC scores etc.
png(paste0(outGraph, "QC_genePlot.png"), width = 7, height = 5
  , units = "in", res = 300)
par(mfrow = c(1, 2))
GenePlot(centSO, "nUMI", "percent.mito", cex.use = 0.2)
GenePlot(centSO, "nUMI", "nGene", cex.use = 0.2)
dev.off()

## Density plots
# nGene
ggplot(centSO@meta.data, aes(x = nGene)) +
  geom_density(aes(color = "Filtered")) +
  geom_density(data = metDF, aes(x = nGene, color = "All")) +
  labs(color = "Cells") +
  xlab("Number of genes detected") +
  ylab("Density") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nDensity plot - number of genes detected per cell (counts > 0)"
    , "\nFilters:"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > ", high.thresholds, " (3 SD) genes detected"
    , "\nRemove cells with MT transcripts > 5%"
    , "\nMean genes detected after filters: ", round(mean(centSO@meta.data$nGene),1)
    , "\n"
  ))
ggsave(paste0(outGraph, "QC_Density_nGene.pdf"), height = 5, width = 8)
# nUMI
ggplot(metDF, aes(x = nUMI)) +
  geom_density(aes(color = "All")) +
  geom_density(data = centSO@meta.data, aes(x = nUMI, color = "Filtered")) +
  coord_cartesian(xlim = c(0, 10000)) +
  labs(color = "Cells") +
  xlab("Number of UMIs mapping to exons") +
  ylab("Density") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nDensity plot - number of UMIs (reads after UMI collapse) mapping to exons"
    , "\nFilters:"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > ", high.thresholds, " (3 SD) genes detected"
    , "\nRemove cells with MT transcripts > 5%"
    , "\nMean nUMI after filters: ", round(mean(centSO@meta.data$nUMI),1)
    , "\n"
  ))
ggsave(paste0(outGraph, "QC_Density_nUMI.pdf"), height = 5, width = 8)

# Raw counts
# Randomly sample 500 cells
idx <- sample(1:ncol(exDF), 500, replace = TRUE)
ggDF <- exDF[ ,idx]
ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  geom_histogram(binwidth = 1) +
  # geom_density() +
  xlab("Raw counts") +
  ylab("Counts") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nHistogram - raw counts for each gene and each cell"
    , "\nRandomly sample 500 cells"
    , "\n"
  ))
ggsave(paste0(outGraph, "QC_Density_RawCounts.pdf"), height = 5, width = 8)
# Normalized expression values
idx <- sample(1:ncol(noCentExM), 500, replace = TRUE)
ggDF <- noCentExM[ ,idx]
ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  geom_histogram() +
  xlab("Normalized expression") +
  ylab("Counts") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nHistogram - normalized expression for each gene and each cell"
    , "\nRandomly sample 500 cells"
    , "\nFilters:"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > ", high.thresholds, " (3 SD) genes detected"
    , "\nRemove cells with MT transcripts > 5%"
    , "\n"
  ))
ggsave(paste0(outGraph, "QC_Density_NormExpression.pdf"), height = 5, width = 8)
# Normalized centered scaled expression values
idx <- sample(1:ncol(centSO@scale.data), 500, replace = TRUE)
ggDF <- centSO@scale.data[ ,idx]
ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  geom_histogram() +
  xlab("Normalized expression") +
  ylab("Counts") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nHistogram - normalized centered scaled expression for each gene and each cell"
    , "\nRandomly sample 500 cells"
    , "\nFilters:"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > ", high.thresholds, " (3 SD) genes detected"
    , "\nRemove cells with MT transcripts > 5%"
    , "\n"
  ))
rm(ggDF)
ggsave(paste0(outGraph, "QC_Density_NormExpressionCentScale.pdf"), height = 5, width = 8)
################################################################################

### Detection of variable genes across the single cells

print("### Detection of variable genes across the single cells")

# Seurat calculates highly variable genes and focuses on these for downstream
# analysis. FindVariableGenes calculates the average expression and dispersion
# for each gene, places these genes into bins, and then calculates a z-score for
# dispersion within each bin. This helps control for the relationship between
# variability and average expression. This function is unchanged from (Macosko
# et al.), but new methods for variable gene expression identification are
# coming soon. We suggest that users set these parameters to mark visual
# outliers on the dispersion plot, but the exact parameter settings may vary
# based on the data type, heterogeneity in the sample, and normalization
# strategy. The parameters here identify ~2,000 variable genes, and represent
# typical parameter settings for UMI data that is normalized to a total of 1e4
# molecules.

png(paste0(outGraph, "MeanVarPlot.png"), width = 7, height = 5
  , units = "in", res = 300)
centSO <- FindVariableGenes(centSO, mean.function = ExpMean
  , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3
  , y.cutoff = 0.5)
dev.off()

print("Number of variable genes used for clustering:")
length(centSO@var.genes)

## Remove CC genes from variable genes

# Cleanup CC marker data frame
ccDF <- data.frame(lapply(ccDF, as.character), stringsAsFactors=FALSE)
cc <- c(unlist(ccDF))
cc <- gsub(" *", "", cc)
idx <- centSO@var.genes %in% cc
centSO@var.genes <- centSO@var.genes[! idx]

print("Number of variable genes used for clustering after removing CC genes:")
length(centSO@var.genes)

save(centSO, noCentExM, metDF, high.thresholds, file = paste0(outData, "seuratO.Robj"))
################################################################################

### Perform linear dimensional reduction

print("### Perform linear dimensional reduction")

# Perform PCA on the scaled data. By default, the genes in object\@var.genes are
# used as input, but can be defined using pc.genes. We have typically found that
# running dimensionality reduction on genes with high-dispersion can improve
# performance. However, with UMI data - particularly after using ScaleData, we
# often see that PCA returns similar (albeit slower) results when run on much
# larger subsets of genes, including the whole transcriptome.
# centSO <- PCA(centSO, pc.genes = centSO@var.genes, do.print = TRUE
#   , pcs.print = 5, genes.print = 5)

# Run PCA with the IRLBA package (iteratively computes the top dimensions,
# dramatic increase in speed since we only use a fraction of the PCs anyways) if
# you see the warning "did not converge–results might be invalid!; try
# increasing maxit or fastpath=FALSE", try increasing maxit
centSO <- RunPCA(object = centSO, pc.genes = centSO@var.genes, pcs.compute = 50
  , pcs.print = 1:5, genes.print = 5, maxit = 500, weight.by.var = FALSE)

# ProjectPCA scores each gene in the dataset (including genes not included in
# the PCA) based on their correlation with the calculated components. Though we
# don't use this further here, it can be used to identify markers that are
# strongly correlated with cellular heterogeneity, but may not have passed
# through variable gene selection. The results of the projected PCA can be
# explored by setting use.full=T in the functions below.
# Saves in slot @pca.x.full
centSO <- ProjectPCA(centSO)

# Seurat provides several useful ways of visualizing both cells and genes that
# define the PCA, including **PrintPCA()**, **VizPCA()**, **PCAPlot()**, and
# **PCHeatmap()** Examine  and visualize PCA results a few different ways
pdf(paste0(outGraph, "PCAplot.pdf"))
PrintPCA(centSO, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
# Returns the top genes ranked by the score's absolute values
# From @pca.obj[[1]]$rotation
# VizPCA(centSO, 1:4)
PCAPlot(centSO, 1, 2, pt.size = 0.25)
# In particular **PCHeatmap()** allows for easy exploration of the primary
# sources of heterogeneity in a dataset, and can be useful when trying to decide
# which PCs to include for further downstream analyses. Both cells and genes are
# ordered according to their PCA scores. Setting cells.use to a number plots the
# 'extreme' cells on both ends of the spectrum, which dramatically speeds
# plotting for large datasets. Though clearly a supervised analysis, we find
# this to be a valuable tool for exploring correlated gene sets.
# PCHeatmap(centSO, pc.use = 1, do.balanced = TRUE)
# PCHeatmap(centSO, pc.use = 1:6, do.balanced = TRUE, label.columns = FALSE
#   , use.full = FALSE)
dev.off()

# Plot genes with highest PC loadings
# Genes with highest PC loadings
ggL <- lapply(1:8, function(pc) {
  df <- rbind(data.frame(PC = sort(centSO@dr$pca@gene.loadings[ ,pc])[1:10])
    , data.frame(PC = sort(centSO@dr$pca@gene.loadings[ ,pc], decreasing = TRUE)[1:10])
  )
  df$GENE <- factor(row.names(df), levels = row.names(df))
  gg <- ggplot(df, aes(x = PC, y = GENE)) +
    geom_point() +
    xlab("Loading") +
    ylab("Gene") +
    ggtitle(paste0("PC ", pc))
  return(gg)
})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nGenes with highest PC loadings"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PCAloadingPlots.pdf"), width = 13
  , height = 13, limitsize = FALSE)

save(centSO, noCentExM, metDF, high.thresholds, file = paste0(outData, "seuratO.Robj"))
################################################################################

### Determine statistically significant principal components

# To overcome the extensive technical noise in any single gene for scRNA-seq
# data, Seurat clusters cells based on their PCA scores, with each PC
# essentially representing a 'metagene' that combines information across a
# correlated gene set. Determining how many PCs to include downstream is
# therefore an important step.

# In Macosko et al, we implemented a resampling test inspired by the jackStraw
# procedure. We randomly permute a subset of the data (1% by default) and rerun
# PCA, constructing a 'null distribution' of gene scores, and repeat this
# procedure. We identify 'significant' PCs as those who have a strong enrichment
# of low p-value genes.

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
# centSO <- JackStraw(centSO, num.replicate = 100, do.print = FALSE)

# The **JackStrawPlot()** function provides a visualization tool for comparing
# the distribution of p-values for each PC with a uniform distribution (dashed
# line). 'Significant' PCs will show a strong enrichment of genes with low
# p-values (solid curve above the dashed line). In this case it appears that PCs
# 1-10 are significant.

# JackStrawPlot(centSO, PCs = 1:50)

# A more ad hoc method for determining which PCs to use is to look at a plot of
# the standard deviations of the principle components and draw your cutoff where
# there is a clear elbow in the graph. This can be done with **PCElbowPlot()**.
# In this example, it looks like the elbow would fall around PC 9.

pdf(paste0(outGraph, "PCElbowPlot.pdf"))
PCElbowPlot(centSO, num.pc = 50)
dev.off()

# # Calculate variance explained by each PC
# varExp <- (pcNmVar$sdev)^2 / sum(pcNmVar$sdev^2)
# topVar <- varExp[1:5]
# colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
#   , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

#PC selection -- identifying the true dimensionality of a dataset -- is an
#important step for Seurat, but can be challenging/uncertain for the user. We
#therefore suggest these three approaches to consider. The first is more
#supervised, exploring PCs to determine relevant sources of heterogeneity, and
#could be used in conjunction with GSEA for example. The second implements a
#statistical test based on a random null model, but is time-consuming for large
#datasets, and may not return a clear PC cutoff. The third is a heuristic that
#is commonly used, and can be calculated instantly. In this  example all three
#approaches yielded similar results, but we might have been justified in
#choosing anything between PC 7-10 as a cutoff. We followed the jackStraw  here,
#admittedly buoyed by seeing the PCHeatmap returning interpretable signals
#(including canonical dendritic cell markers) throughout these PCs. Though the
#results are only subtly affected by small shifts in this cutoff (you can test
#below), we strongly suggest always explore the PCs they choose to include
#downstream.
################################################################################

### Cluster the cells

# Seurat now includes an graph-based clustering approach compared to (Macosko et
# al.). Importantly, the distance metric which drives the clustering analysis
# (based on previously identified PCs) remains the same. However, our approach
# to partioning the cellular distance matrix into clusters has dramatically
# improved. Our approach was heavily inspired by recent manuscripts which
# applied graph-based clustering approaches to scRNA-seq data [SNN-Cliq, Xu and
# Su, Bioinformatics, 2015] and CyTOF data [PhenoGraph, Levine et al., Cell,
# 2015]. Briefly, these methods embed cells in a graph structure - for example a
# K-nearest neighbor (KNN) graph, with edges drawn between cells with similar
# gene expression patterns, and then attempt to partition this graph into highly
# interconnected ‘quasi-cliques’ or ‘communities’. As in PhenoGraph, we first
# construct a KNN graph based on the euclidean distance in PCA space, and refine
# the edge weights between any two cells based on the shared overlap in their
# local neighborhoods (Jaccard distance). To cluster the cells, we apply
# modularity optimization techniques [SLM, Blondel et al., Journal of
# Statistical Mechanics], to iteratively group cells together, with the goal of
# optimizing the standard modularity function.

# The FindClusters function implements the procedure, and contains a resolution
# parameter that sets the ‘granularity’ of the downstream clustering, with
# increased values leading to a greater number of clusters. We find that setting
# this parameter between 0.6-1.2 typically returns good results for single cell
# datasets of around 3K cells. Optimal resolution often increases for larger
# datasets. The clusters are saved in the object@ident slot.

# save.SNN=T saves the SNN so that the  SLM algorithm can be rerun using the
# same graph, but with a different resolution value (see docs for full details)
# centSO <- FindClusters(centSO, dims.use = 1:10, resolution = 0.6
# , print.output = 0, save.SNN = T)
################################################################################

### Run Non-linear dimensional reduction (tSNE)

print("### Run Non-linear dimensional reduction (tSNE)")

# Seurat continues to use tSNE as a powerful tool to visualize and explore these
# datasets. While we no longer advise clustering directly on tSNE components,
# cells within the graph-based clusters determined above should  co-localize on
# the tSNE plot. This is because the tSNE aims to place cells with similar local
# neighborhoods in high-dimensional space together in low-dimensional space. As
# input to the tSNE, we suggest using the same PCs as input to the clustering
# analysis, although computing the tSNE based on scaled gene expression is also
# supported using the genes.use argument.

## Cluster and tSNE with different parameters

## PCs

# PC 1-6
centSO <- RunTSNE(centSO, dims.use = 1:6, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:6, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg1 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg1 <- gg1 + ggtitle("PCs: 1-6")

# PC 1-10
centSO <- RunTSNE(centSO, dims.use = 1:10, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:10, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg2 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg2 <- gg2 + ggtitle("PCs: 1-10")

# PC 1-15
centSO <- RunTSNE(centSO, dims.use = 1:15, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:10, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg3 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg3 <- gg3 + ggtitle("PCs: 1-15")

# PC 1-20
centSO <- RunTSNE(centSO, dims.use = 1:20, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:20, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg4 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg4 <- gg4 + ggtitle("PCs: 1-20")

# PC 1-30
centSO <- RunTSNE(centSO, dims.use = 1:30, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:30, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg5 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg5 <- gg5 + ggtitle("PCs: 1-30")

# PC 1-40
centSO <- RunTSNE(centSO, dims.use = 1:40, do.fast = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:40, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
gg6 <- TSNEPlot(centSO, pt.size = 0.5, do.return = TRUE, do.label = TRUE)
gg6 <- gg6 + ggtitle("PCs: 1-40")

# plot_grid combine tSNE graphs
pg <- plot_grid(gg1, gg2, gg3, gg4, gg5, gg6, ncol = 2, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSeurat cluster and tSNE with different PCs"
  , "\n"
  , "\nRemove genes > 0 counts in < 3 cells"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nNormalize expression"
  , "\nRegress out covariates"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.25, 1))
ggsave(paste0(outGraph, "tSNE_PCsTest.png")
  , width = 12, height = 18)

## Resolution
Seurat_Resolution_Test <- function(seuratO, resolutions){
  # centSO <- RunTSNE(centSO, dims.use = 1:40, do.fast = TRUE)
  ggL <- lapply(resolutions, function(resolution) {
      so <- FindClusters(seuratO, dims.use = 1:40, resolution = resolution
        , print.output = 0, save.SNN = TRUE)
      gg <- TSNEPlot(so, pt.size = 0.01, do.return = TRUE, do.label = TRUE)
      gg <- gg + ggtitle(paste0("Resolution: ", resolution))
      return(gg)
    })
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 2, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different Seurat cluster resolutions"
    , "\n"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > 3192 (3 SD) genes detected"
    , "\nRemove genes detected in < 3 cells"
    , "\nNormalize expression"
    , "\nRegress out covariates"
    , "\n"))
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
  return(pg)
}
Seurat_Resolution_Test(
  seuratO = centSO
  , resolutions = c(
      0.2, 0.3, 0.4, 0.5, 0.54, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0)
)
ggsave(paste0(outGraph, "tSNE_ResolutionTest.png")
  , width = 12, height = 40)
Seurat_Resolution_Test(
  seuratO = centSO
  , resolutions = c(0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.6)
)
ggsave(paste0(outGraph, "tSNE_ResolutionTest0506.png")
  , width = 12, height = 22)

## Clustering parameters to save
centSO <- RunTSNE(centSO, dims.use = 1:40, do.fast = TRUE, force.recalc = TRUE)
centSO <- FindClusters(centSO, dims.use = 1:40, resolution = 0.54
  , print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

PrintFindClustersParams(object = centSO)

TSNEPlot(centSO, do.label = FALSE, pt.size = 0.6, do.return = TRUE
  , no.legend = TRUE)
ggsave(paste0(outGraph, "tSNE_NoLabels.png")
  , width = 7, height = 7)

## Plot tSNE graph colored by lanes, samples, or GZ CP
# Collect tSNE values
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

# metDF$LIBRARY <- "Plath"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N702"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N701"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-007" & metDF$NEXTERA == "N710"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-007" & metDF$NEXTERA == "N711"] <- "Geschwind"

# Add cluster identity
ggDF$CLUSTER <- centSO@ident
# Add metadata
ggDF <- data.frame(ggDF, metDF[match(row.names(ggDF), metDF$CELL), ])
ggDF$BRAIN <- as.factor(ggDF$BRAIN)

# Add DS-002-003 cluster identity
# Assign annotated cluster names to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
new.cluster.ids <- c(
  "Excitatory Upper Layer Neuron 1"
  , "Excitatory Neuron"
  , "Excitatory Upper Layer Neuron 2"
  , "Excitatory Deep Layer Neuron"
  , "Intermediate Progenitors"
  , "Interneuron"
  , "Mitotic Progenitors"
  , "oRG"
  , "Oligodendrocyte Precursor"
  , "Endothelial")
fetb@ident <- plyr::mapvalues(fetb@ident, from = current.cluster.ids
  , to = new.cluster.ids)
# Stash numerical cluster identities if want to use later
# fetb <- StashIdent(fetb, save.name = "Cluster_Numbers")
preCl <- fetb@ident
names(preCl) <- gsub("_.*", "", names(preCl))
idx <- match(ggDF$CELL, names(preCl))
ggDF$PREVIOUS_CLUSTER <- preCl[idx]

# ggplot clustering
ggC <- TSNEPlot(centSO, do.label = TRUE, pt.size = 0.1, do.return = TRUE
  , no.legend = FALSE, alpha = 0.5)
ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 7)))
ggC <- ggC + ggtitle("Seurat clustering")
# ggplot previous clustering
gg2DF <- ggDF[! is.na(ggDF$PREVIOUS_CLUSTER), ]
ggP <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PREVIOUS_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_point(data = gg2DF, size = 0.3, alpha = 0.5
    , aes(x = tSNE_1, y = tSNE_2, col = PREVIOUS_CLUSTER)) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  ggtitle("Seurat clustering from Donor 2")
# ggsave(paste0(outGraph, "tSNE_PreviousClustering_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Sample
ggB <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = BRAIN)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  ggtitle("Donor")
# ggsave(paste0(outGraph, "tSNE_Sample_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Seq Run
ggS <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = SEQ_RUN)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  ggtitle("Sequencing run")
# ggsave(paste0(outGraph, "tSNE_Lane_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Region
ggR <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = REGION)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  ggtitle("Region (GZ / CP)")
# ggplot Library
ggL <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = LIBRARY)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  ggtitle("Lab that did library prep")

# plot grid
pg <- plot_grid(ggC, ggR, ggP, ggB, ggL, align = "v", axis = "l", ncol = 1, nrow = 5)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSeurat cluster and tSNE colored by clusters or covariates"
  , "\n"
  , "\nRemove genes > 0 counts in < 3 cells"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nNormalize expression"
  , "\nRegress out covariates"
  , "\ntSNE + Seurat cluster PC 1-40"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
# save
ggsave(paste0(outGraph, "tSNE_Metadata_PC1-40.png"), width = 8, height = 28)

# Save as grid
png(paste0(outGraph, "tSNE_MetadataGrid_PC1-40.png")
  , width = 12, height = 16, units = "in", res = 300)
# plot_grid(ggC, ggP, ggB, ggS, ggL, ggR, ncol = 2)
layout <- rbind(c(1,1,2,2), c(3,3,3,NA), c(4,4,5,5))
grid.arrange(ggC, ggR, ggP, ggB, ggL, layout_matrix = layout)
dev.off()

# tSNE for paper
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# Add cluster identity
ggDF$CLUSTER <- centSO@ident
# Add metadata
ggDF <- data.frame(ggDF, metDF[match(row.names(ggDF), metDF$CELL), ])
# Plot
gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = CLUSTER)) +
  geom_point(size = 0.01, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(
    panel.grid.major = element_blank(),
    , panel.grid.minor = element_blank()
  ) +
  ggtitle("Seurat clusters")
gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = REGION)) +
  geom_point(size = 0.01, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(
    panel.grid.major = element_blank(),
    , panel.grid.minor = element_blank()
  ) +
  ggtitle("Region (GZ / CP)")
plot_grid(gg1, gg2, ncol = 2)
ggsave(paste0(outGraph, "tSNE_PC1-40_paper.png")
  , dpi = 150, width = 8, height = 3)

## Cell metadata as csv for paper
mdat_paper_DF <- centSO@meta.data
mdat_paper_DF <- mdat_paper_DF[ ,c("CELL", "res.0.54", "BRAIN", "REGION"
  , "NEXTERA", "LIBRARY", "nGene", "nUMI", "percent.mito"
  , "S.Score", "G2M.Score", "Phase")]
mdat_paper_DF$Cluster_number <- mdat_paper_DF$"res.0.54"
# Cluster annotations
cluster_annot <- c(
  "9" = "vRG"
  , "7" = "oRG"
  , "8" = "Cycling progenitor S phase"
  , "10" = "Cycling progenitor G2/M phase"
  , "2" = "IPC"
  , "0" = "Excitatory neuron new born migrating"
  , "1" = "Excitatory neuron"
  , "4" = "Excitatory neuron (collosal)"
  , "3" = "Deep layer excitatory neuron 1"
  , "13" = "Deep layer excitatory neuron 2"
  , "5" = "Interneuron (SST)"
  , "6" = "Interneuron (CALB2)"
  , "11" = "Oligodendrocyte precursor"
  , "12" = "Endothelial"
  , "14" = "Pericyte"
  , "15" = "Microglia"
  , "16" = "NA"
)
idx <- match(mdat_paper_DF$Cluster_number, names(cluster_annot))
mdat_paper_DF$Cluster <- cluster_annot[idx]
# Donor IDs
donor_annot <- c("2" = "368", "3" = "370", "4" = "371", "5" = "372")
idx <- match(mdat_paper_DF$BRAIN, names(donor_annot))
mdat_paper_DF$BRAIN <- donor_annot[idx]
mdat_paper_DF <- mdat_paper_DF[ ,c(1, 14:13, 3:12)]
# Gestation week
gw_annot <- c("368" = "17", "370" = "18", "371" = "17", "372" = "18")
idx <- match(mdat_paper_DF$BRAIN, names(gw_annot))
mdat_paper_DF$GW <- gw_annot[idx]
# Convert to percentage
mdat_paper_DF$percent.mito <- round(mdat_paper_DF$percent.mito * 100, 2)
# Order by cluster
mdat_paper_DF$Cluster_number <- factor(mdat_paper_DF$Cluster_number
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15,16))
mdat_paper_DF <- mdat_paper_DF[order(mdat_paper_DF$Cluster_number), ]
colnames(mdat_paper_DF) <- c("Cell", "Cluster", "Cluster_number", "Donor"
  , "Region", "Index", "Library", "Number_genes_detected", "Number_UMI"
  , "Percentage_mitochondrial", "S_phase_score", "G2M_phase_score", "Phase"
  , "Gestation_week")
mdat_paper_DF <- mdat_paper_DF[ ,c(1:5,14,6:13)]
write.csv(mdat_paper_DF, file = paste0(outTable, "Cell_Metadata_paper.csv")
  , quote = FALSE, row.names = FALSE)


# You can save the object at this point so that it can easily be loaded back in
# without having to rerun the computationally intensive steps performed above,
# or easily shared with collaborators.
save(centSO, noCentExM, metDF, high.thresholds, file = paste0(outData, "seuratO.Robj"))

# Heatmap of known markers genes
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
# kmDF <- kmDF[! is.na(kmDF$Gene.Symbol), ]
# # # Filter for gene symbols in data
# # kmDF <- kmDF[kmDF$Gene.Symbol %in% centSO@data@Dimnames[1][[1]], ]
# # pdf(paste0(outGraph, "DoHeatmap_KnownMarkers.pdf"), width = 10, height = 4 + 0.2*nrow(kmDF))
# # kmDFL <- split(kmDF, kmDF$Population.Marked)
# # sapply(names(kmDFL), function(kmLName) {
# #   # print(kmLName)
# #   df <- kmDFL[[kmLName]]
# #   # print(df)
# #   print(as.character(df$Gene.Symbol))
# #   DoHeatmap(centSO, genes.use = as.character(df$Gene.Symbol)
# #     , order.by.ident = TRUE
# #     , slim.col.label = TRUE
# #     , remove.key = FALSE
# #     , main = (kmLName)
# #   )
# # })
# # dev.off()
# pdf(paste0(outGraph, "DoHeatmap_KnownMarkers.pdf"), width = 16, height = 4 + 0.2*nrow(kmDF))
# DoHeatmap(centSO, genes.use = as.character(kmDF$Gene.Symbol)
#   , order.by.ident = TRUE
#   , slim.col.label = TRUE
#   , remove.key = FALSE
#   , main = "Known Markers")
# dev.off()
################################################################################

### Output down-sampled Seurat object for testing

# Sample cell IDs
cellIDs <- sample(colnames(centSO@scale.data), 2000)
# Subset
ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
ssCentSO@raw.data <- ssCentSO@raw.data[
  ,colnames(ssCentSO@raw.data) %in% cellIDs]
ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

save(ssCentSO, ssNoCentExM, metDF, high.thresholds
  , file = paste0(outData, "TEST_seuratO.Robj"))

# Sample cell IDs
cellIDs <- sample(colnames(centSO@scale.data), 5000)
# Subset
ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
ssCentSO@raw.data <- ssCentSO@raw.data[
  ,colnames(ssCentSO@raw.data) %in% cellIDs]
ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

save(ssCentSO, ssNoCentExM, metDF, high.thresholds
  , file = paste0(outData, "TEST5000_seuratO.Robj"))

# Sample cell IDs
cellIDs <- names(centSO@ident)[centSO@ident %in% c(7,8,9,10)]
# Subset
ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
ssCentSO@raw.data <- ssCentSO@raw.data[
  ,colnames(ssCentSO@raw.data) %in% cellIDs]
ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

save(ssCentSO, ssNoCentExM, metDF, high.thresholds
  , file = paste0(outData, "TESTcluster78910_seuratO.Robj"))

# Sample cell IDs
cellIDs <- names(centSO@ident)[centSO@ident %in% c(2,7,8,9,10)]
# Subset
ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
ssCentSO@raw.data <- ssCentSO@raw.data[
  ,colnames(ssCentSO@raw.data) %in% cellIDs]
ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

save(ssCentSO, ssNoCentExM, metDF, high.thresholds
  , file = paste0(outData, "TESTcluster278910_seuratO.Robj"))

# Sample cell IDs
cellIDs <- names(centSO@ident)[centSO@ident %in% c(0,2,7,8,9,10)]
# Subset
ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
ssCentSO@raw.data <- ssCentSO@raw.data[
  ,colnames(ssCentSO@raw.data) %in% cellIDs]
ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

save(ssCentSO, ssNoCentExM, metDF, high.thresholds
  , file = paste0(outData, "TESTcluster0278910_seuratO.Robj"))
################################################################################

### Feature plot of top PC scores

print("### Feature plot of top PC scores")

# Feature plots of top PC scores
ggL <- lapply(c(1:9), function(i) {
  # Collect tSNE values
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Add PC score
  ggDF$PC <- centSO@dr$pca@cell.embeddings[ ,i]
  gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PC)) +
    geom_point(size = 0.05, alpha = 0.5) +
    # guides(colour = guide_legend(override.aes = list(size = 7))) +
    scale_color_distiller(name = "PC score", type = "div"
      , palette = 5, direction = -1) +
    ggtitle(paste0("PC: ", i)) +
    theme(text = element_text(size = 14))
  return(gg)
})
# plot grid
pg <- plot_grid(plotlist = ggL, align = "v", axis = "l", ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSeurat cluster and tSNE colored by PC scores"
  , "\n"
  , "\nRemove genes > 0 counts in < 3 cells"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nNormalize expression"
  , "\nRegress out covariates"
  , "\ntSNE + Seurat cluster PC 1-40"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# save
ggsave(paste0(outGraph, "PCscore_FeaturePlot.png"), width = 15, height = 15)
################################################################################

### Heatmap of expression of genes with highest PC loadings sorted by cluster

# print("### Heatmap of expression of genes with highest PC loadings sorted by cluster")
#
# # Genes with highest PC loadings
# ldf <- lapply(1:20, function(pc) {
#   data.frame(GENES = names(sort(abs(centSO@dr$pca@gene.loadings[ ,pc])
#     , decreasing = TRUE)[1:10]), PC = paste0("PC", pc))
# })
# genesDF <- do.call("rbind", ldf)
#
# # # Genes with highest PC loadings
# # ldf <- lapply(1:40, function(pc) {
# #   data.frame(GENES = names(sort(abs(centSO@pca.obj[[1]]$rotation[ ,pc])
# #     , decreasing = TRUE)[1:10]), PC = paste0("PC", pc))
# # })
# # genesDF <- do.call("rbind", ldf)
#
# # # Genes with highest PC loadings - if PCAFast is used to calculate PCs
# # m <- centSO@pca.obj[[1]]$v
# # row.names(m) <- centSO@var.genes
# # ldf <- lapply(1:40, function(pc) {
# #   data.frame(GENES = row.names(m[order(-abs(m[,pc])), ])[1:10]
# #     , PC = paste0("PC", pc))
# # })
# # genesDF <- do.call("rbind", ldf)
#
# # Format for heatmap function
# colnames(genesDF) <- c("GENE", "GROUP")
# # Heatmap - mean centered scaled
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF = genesDF
#   , exprM = centSO@scale.data, seuratO = centSO
#   , lowerLimit = -1.5, upperLimit = 1.5
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:length(unique(centSO@ident))))
# # Change legend title
# ggL <- lapply(ggL, function(gg) {
#   gg + scale_fill_distiller(name = "Normalized\nExpression\nz-score")
#   return(gg)
# })
# # plot grid
# pg <- plot_grid(plotlist = ggL, align = "h", axis = "b", ncol = 3)
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n"
#   , "\nExpression of genes with highest PC loadings"
#   , "\nMean centered, variance scaled, normalized expression"
#   , "\nLimits set to -1.5 and 1.5"
#   , "\n"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
# # save
# ggsave(paste0(outGraph, "PC1-40_DoHeatmap_HighestLoading_CenterScale.png")
#   , width = 12, height = 4 + 0.16*nrow(genesDF))
#
# # Heatmap - not mean centered scaled
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF = genesDF
#   , exprM = noCentExM, seuratO = centSO
#   , lowerLimit = -1, upperLimit = 3
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:length(unique(centSO@ident))))
# # plot grid
# pg <- plot_grid(plotlist = ggL, align = "h", axis = "b", ncol = 3)
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n"
#   , "\nExpression of genes with highest PC loadings"
#   , "\nMean centered, variance scaled, normalized expression"
#   , "\nLimits set to -1 and 3"
#   , "\n"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
# # save
# ggsave(paste0(outGraph, "PC1-40_DoHeatmap_HighestLoading_NoCentScale.png")
#   , width = 12, height = 4 + 0.16*nrow(genesDF))
################################################################################

### Cluster QC and statistics

# nUMI and nGene per cluster violin plot
png(paste0(outGraph, "nGenenUMIperCluster_violinPlot_PC1-40.png")
  , width = 10, height = 5, units = "in", res = 300)
VlnPlot(centSO, c("nGene", "nUMI")
  , nCol = 2, point.size.use = 0.02)
dev.off()

# nGene boxplot
df <- data.frame(
  Number_Of_Genes_Detected = centSO@meta.data$nGene
  , Number_UMI = centSO@meta.data$nUMI
  , ClusterID = centSO@ident)
ggplot(df, aes(x = ClusterID, y = Number_Of_Genes_Detected)) +
  geom_boxplot() +
  xlab("Number of genes detected") +
  ylab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of genes detected per cluster"
    , "\n(> 0 counts)"))
ggsave(paste0(outGraph, "PC1-40_nGeneCluster_Boxplot.pdf"), width = 5, height = 5)

# Metrics per cluster
df <- data.frame(
  CLUSTER = names(table(centSO@ident))
  , NUMBER_CELLS = as.vector(table(centSO@ident))
  , NUMBER_CELLS_CP = as.vector(table(centSO@ident[centSO@meta.data$REGION == "CP"]))
  , NUMBER_CELLS_GZ = as.vector(table(centSO@ident[centSO@meta.data$REGION == "GZ"]))
  , MEAN_UMI = tapply(centSO@meta.data$nUMI, centSO@ident, mean)
  , MEAN_GENES_DETECTED = tapply(centSO@meta.data$nGene, centSO@ident, mean)
  , PERCENT_CELLS = (as.vector(table(centSO@ident)) /
    sum(as.vector(table(centSO@ident)))) * 100
)
write.csv(df, paste0(outTable, "PC1-40_Cluster_Metrics.csv")
  , quote = FALSE, row.names = FALSE)

# Format for paper
cluster_metrics_DF <- df
cluster_metrics_DF$Cluster_number <- df$CLUSTER
cluster_annot <- c(
  "9" = "vRG"
  , "7" = "oRG"
  , "8" = "Cycling progenitor S phase"
  , "10" = "Cycling progenitor G2/M phase"
  , "2" = "IPC"
  , "0" = "Excitatory neuron new born migrating"
  , "1" = "Excitatory neuron"
  , "4" = "Excitatory neuron (collosal)"
  , "3" = "Deep layer excitatory neuron 1"
  , "13" = "Deep layer excitatory neuron 2"
  , "5" = "Interneuron (SST)"
  , "6" = "Interneuron (CALB2)"
  , "11" = "Oligodendrocyte precursor"
  , "12" = "Endothelial"
  , "14" = "Pericyte"
  , "15" = "Microglia"
  , "16" = "NA"
)
idx <- match(cluster_metrics_DF$Cluster_number, names(cluster_annot))
cluster_metrics_DF$Cluster <- cluster_annot[idx]
# Order by cluster
cluster_metrics_DF$Cluster_number <- factor(cluster_metrics_DF$Cluster_number
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15,16))
cluster_metrics_DF <- cluster_metrics_DF[
  order(cluster_metrics_DF$Cluster_number), ]
# Order columns
cluster_metrics_DF <- cluster_metrics_DF[ ,c(9,8,2:7)]
# Round
cluster_metrics_DF[ ,6:8] <- round(cluster_metrics_DF[ ,6:8], 1)
colnames(cluster_metrics_DF) <- c("Cluster", "Cluster number"
  , "Number of cells", "Number of cells from CP", "Number of cells from GZ"
  , "Mean number UMI per cell", "Mean genes detected per cell"
  , "Percent of total cells")
write.csv(cluster_metrics_DF
  , file = paste0(outTable, "Cluster_metrics_paper.csv")
  , quote = FALSE, row.names = FALSE)

# Dot plot of percent of cell classes
v1 <- c(
  "Radial glia" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(7,9)])/sum(df$NUMBER_CELLS) * 100
  , "Cycling progenitor" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(8,10)])/sum(df$NUMBER_CELLS) * 100
  , "Intermediate progenitor" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(2)])/sum(df$NUMBER_CELLS) * 100
  , "Excitatory Neuron" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(0,1,3,4,13)])/sum(df$NUMBER_CELLS) * 100
  , "Interneuron" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(5,6)])/sum(df$NUMBER_CELLS) * 100
  , "Oligodendrocyte precursor" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(11)])/sum(df$NUMBER_CELLS) * 100
  , "Endothelial" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(12)])/sum(df$NUMBER_CELLS) * 100
  , "Pericyte" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(14)])/sum(df$NUMBER_CELLS) * 100
  , "Microglia" = sum(df$NUMBER_CELLS[df$CLUSTER %in% c(16)])/sum(df$NUMBER_CELLS) * 100
)
df1 <- data.frame(Class = names(v1), Percent = v1)
df1$Class <- factor(df1$Class, levels = c("Radial glia", "Cycling progenitor"
  , "Intermediate progenitor"
  , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
  , "Endothelial", "Pericyte", "Microglia"))
ggplot(df1, aes(x = Percent, y = Class, color = Class)) +
geom_point(size = 3) +
scale_color_brewer(type = "qual", palette = "Set3", direction = 1) +
scale_y_discrete(limits = rev(levels(df1$Class))) +
theme(legend.position = "none") +
ylab("") +
ggtitle(paste0(graphCodeTitle
  , "\nPercentages of cell classes"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Percent_CellType.png")
  , height = 2.5, width = 4)

# Bar plot of percent of cells in clusters
df1 <- df[df$CLUSTER != 16, ]
df1$CLUSTER <- factor(df1$CLUSTER
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
df1 <- df1[order(df1$CLUSTER), ]
df1$Cluster2 <- as.factor(c(0:15))
df1$Class = factor(c(
      rep("Radial glia", 2)
      , rep("Cycling progenitor", 2)
      , "Intermediate progenitor"
      , rep("Excitatory Neuron", 5)
      , rep("Interneuron", 2)
      , "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia"
    )
  , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
    , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia")
)
ggplot(df1, aes(x = Cluster2, y = PERCENT_CELLS, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  # scale_y_discrete(limits = rev(levels(df1$Cluster2))) +
  theme(legend.position = "none") +
  xlab("Cluster") +
  ylab("Percent") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  ggtitle(paste0(graphCodeTitle
    , "\nPercentages of cell classes"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Percent_CellType_Barplot.pdf")
  , height = 2.5, width = 6)

# GZ/CP ratio bar plot
df1 <- data.frame(Cluster = df$CLUSTER
  , GZCP_Log2_Ratio = log(df$NUMBER_CELLS_CP/df$NUMBER_CELLS_GZ, 2)
)
df1 <- df1[df1$Cluster != 16, ]
df1$Cluster <- factor(df1$Cluster
  , levels = c(c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)))
df1 <- df1[order(df1$Cluster), ]
df1$Cluster2 = as.factor(c(0:15))
df1$Class = factor(c(
      rep("Radial glia", 2)
      , rep("Cycling progenitor", 2)
      , "Intermediate progenitor"
      , rep("Excitatory Neuron", 5)
      , rep("Interneuron", 2)
      , "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia"
    )
  , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
    , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia")
)
ggplot(df1, aes(x = Cluster2, y = GZCP_Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 ratio number of cell GZ / CP for each cluster"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_GZCP_Ratio.png")
  , height = 3, width = 5.5)
# For paper
ggplot(df1, aes(x = Cluster2, y = GZCP_Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  theme(
    legend.position = "none"
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  )
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_GZCP_Ratio_paper.pdf")
  , height = 2.5, width = 6)

# Fraction of cells from each donor per cluster
ggDF <- centSO@meta.data
ggDF$Cluster <- ggDF$"res.0.54"
ggDF <- ggDF[ggDF$Cluster != 16, ]
# Reorder clusters
cluster_key_DF <- data.frame(
  Original = c(0:15)
  , Reorder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
idx <- match(as.numeric(ggDF$Cluster), cluster_key_DF$Original)
ggDF$Cluster2 <- cluster_key_DF$Reorder[idx]
# Calculate percent
ggDF <- table(ggDF$Cluster, ggDF$individual)
ggDF <- (ggDF / rowSums(ggDF)) * 100
# Format
ggDF <- melt(ggDF)
colnames(ggDF) <- c("Cluster", "Donor", "Percent_of_cells")
# Plot
ggplot(ggDF, aes(x = Cluster, y = Percent_of_cells, fill = Donor)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set2", direction = 1) +
  theme(
    legend.position = "none"
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  ylab("Percent of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\nFraction of cells from each donor per cluster"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Donor_Percent.png")
  , height = 2.5, width = 6)

# For paper
ggDF <- centSO@meta.data
ggDF$Cluster <- ggDF$"res.0.54"
ggDF <- ggDF[ggDF$Cluster != 16, ]
# Reorder clusters
cluster_key_DF <- data.frame(
  Original = c(0:15)
  , Reorder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
idx <- match(as.numeric(ggDF$Cluster), cluster_key_DF$Original)
ggDF$Cluster2 <- cluster_key_DF$Reorder[idx]
# Calculate percent
ggDF <- table(ggDF$Cluster2, ggDF$individual)
ggDF <- (ggDF / rowSums(ggDF)) * 100
# Format
ggDF <- melt(ggDF)
colnames(ggDF) <- c("Cluster", "Donor", "Percent_of_cells")
# Plot
ggplot(ggDF, aes(x = Cluster, y = Percent_of_cells, fill = Donor)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Set2", direction = 1) +
  theme(
    legend.position = "none"
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  ylab("Percent of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\nFraction of cells from each donor per cluster"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Donor_Percent_paper.pdf")
  , height = 2.5, width = 6)

# Number of cells per cluster versus mean genes detected per cluster
ggDF <- data.frame(table(centSO@ident)
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggplot(ggDF, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of cells in cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\nGenes detected per cluster versus number of cells per cluster"))
ggsave(paste0(outGraph, "PC1-40_nCellVsnGene_ScatterPlot.pdf")
  , width = 7, height = 7)

# Variance of gene expression per cluster
Plot_Expr_Variance_By_Cluster <- function(number_top_genes){
  expr_DF <- melt(noCentExM)
  idx <- match(expr_DF$Var2, names(centSO@ident))
  expr_DF$Cluster <- centSO@ident[idx]
  expr_DFL <- split(expr_DF, expr_DF$Cluster)
  var_DFL <- lapply(expr_DFL, function(ss_expr_DF){
    mn_expr_DF <- aggregate(value~Var1, ss_expr_DF, mean)
    mn_expr_DF[order(-mn_expr_DF$value), ]
    top_genes <- mn_expr_DF$Var1[1:number_top_genes]
    ss_expr_DF <- with(ss_expr_DF, ss_expr_DF[Var1 %in% top_genes, ])
    var_DF <- aggregate(value~Var1+Cluster, ss_expr_DF, var)
    return(var_DF)
  })
  var_DF <- do.call("rbind", var_DFL)
  ggplot(var_DF, aes(x = Cluster, y = value, fill = Cluster)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0, 1)) +
    xlab("Cluster") +
    ylab("Variance")
}
Plot_Expr_Variance_By_Cluster(number_top_genes = 1000) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nVariance of expression for"
    , "\ntop 1000 most highly expressed genes in each cluster"
    , "\nSubset dataset to 5000 cells"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Variance_Top1000.png")
  , width = 7, height = 7)
Plot_Expr_Variance_By_Cluster(number_top_genes = 35543) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nVariance of expression for all genes"
    , "\nSubset dataset to 5000 cells"))
ggsave(paste0(outGraph, "PC1-40_Cluster_Metrics_Variance.png")
  , width = 7, height = 7)

# Variance of gene expression by cell type
Plot_Expr_Variance_By_Cell_Type <- function(){
  expr_DF <- melt(noCentExM)
  idx <- match(expr_DF$Var2, names(centSO@ident))
  expr_DF$Cluster <- centSO@ident[idx]
  expr_DF$Cell_Type <- NA
  expr_DF$Cell_Type[expr_DF$Cluster %in% c(7,8,9,10)] <- "Progenitors"
  expr_DF$Cell_Type[expr_DF$Cluster %in% c(3,4)] <- "Maturing EN"
  expr_DFL <- split(expr_DF, expr_DF$Cell_Type)
  var_DFL <- lapply(expr_DFL, function(ss_expr_DF){
    mn_expr_DF <- aggregate(value~Var1, ss_expr_DF, mean)
    mn_expr_DF[order(-mn_expr_DF$value), ]
    top_genes <- mn_expr_DF$Var1[1:1000]
    ss_expr_DF <- with(ss_expr_DF, ss_expr_DF[Var1 %in% top_genes, ])
    var_DF <- aggregate(value~Var1+Cell_Type, ss_expr_DF, var)
    return(var_DF)
  })
  var_DF <- do.call("rbind", var_DFL)
  ggplot(var_DF, aes(x = Cell_Type, y = value, fill = Cell_Type)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0, 1)) +
    xlab("Cell types") +
    ylab("Variance") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nVariance of expression for top 1000 most highly expressed genes"
      , "\nProgenitors (vRG, oRG, cycling progenitors 7,8,9,10)"
      , "\nvs maturing excitatory neurons (Deep layer and callosal 3,4)"
      , "\nSubset dataset to 5000 cells"))
}
Plot_Expr_Variance_By_Cell_Type()
ggsave(paste0(outGraph, "PC1-40_Variance_ProgenitorvsNeuron.png")
  , width = 7, height = 7)

# TBR1 expression per cluster
ggDF <- noCentExM
ggDF <- data.frame(EXPRESSION = ggDF[row.names(ggDF) == "TBR1", ])
ggDF$CLUSTER <- centSO@ident
png(paste0(outGraph, "violinPlot_PC1-40_TBR1.png")
  , width = 7, height = 5, units = "in", res = 300)
ggplot(ggDF, aes(x = CLUSTER, y = EXPRESSION)) +
  geom_violin(aes(fill = CLUSTER)) +
  geom_jitter(size = 0.05, height = 0, alpha = 0.1) +
  theme(legend.position = "none") +
  ylab("Normalized expression") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nTBR1 expression by cluster"
    , "\nNormalized expression - read depth, regressed nUMI, brain, lablibrary"
    , "\n"))
dev.off()
################################################################################

### Hierarchical cluster by Seurat cluster mean expression

Hclust_Cluster_Mean_Expression <- function(){
  df <- data.frame(t(noCentExM))
  df$ClusterID <- centSO@ident
  df <- aggregate(.~ClusterID, df, mean)
  row.names(df) <- df$ClusterID
  df <- df[ ,colnames(df) != "ClusterID"]

  d1 <- dist(df, method = "euclidean") # distance matrix
  fit <- hclust(d1)
  ggdendrogram(fit)
}
Hclust_Cluster_Mean_Expression()
ggsave(paste0(outGraph, "PC1-40_hclust.pdf"))
################################################################################

### Output for sharing with highest compression

# Raw counts
raw_counts_mat <- centSO@raw.data[
  ,colnames(centSO@raw.data) %in% names(centSO@ident)]
raw_counts_mat <- as(as.matrix(raw_counts_mat), "sparseMatrix")

# Normalized expression values
norm_mat <- as(as.matrix(noCentExM), "sparseMatrix")

# Metadata
metadata <- centSO@meta.data
metadata <- metadata[ ,colnames(metadata) %in% c("nGene", "nUMI"
  , "CELL", "NEXTERA", "SEQ_RUN", "BRAIN", "REGION", "LIBRARY"
  , "percent.mito", "res.0.54")]
colnames(metadata) <- c("Number_of_Genes_Detected", "Number_of_UMI", "Cell"
  , "Nextera_Index", "Sequencing_Run", "Donor", "Anatomical_Region"
  , "Library_Lab_Batch", "Percentage_Mitochondrial", "Cluster")
# Convert donor to database donor ID
idx <- list("2" = "368", "3" = "370", "4" = "371", "5" = "372")
metadata$Donor <- idx[match(metadata$Donor, names(idx))]

# Cluster annotations
cluster_annot <- c(
  "9" = "vRG"
  , "7" = "oRG"
  , "8" = "Cycling progenitor S phase"
  , "10" = "Cycling progenitor G2/M phase"
  , "2" = "IPC"
  , "0" = "Excitatory Neuron new born migrating"
  , "1" = "Excitatory Neuron"
  , "4" = "Excitatory Neuron (collosal)"
  , "3" = "Deep layer excitatory neuron 1"
  , "13" = "Deep layer excitatory neuron 2"
  , "5" = "Interneuron (SST)"
  , "6" = "Interneuron (CALB2)"
  , "11" = "Oligodendrocyte precursor"
  , "12" = "Endothelial"
  , "14" = "Pericyte"
  , "15" = "Microglia"
)
cluster_annot <- data.frame(
  Cluster = names(cluster_annot), Annotation = cluster_annot
)

save(raw_counts_mat, norm_mat, metadata, cluster_annot
  , file = "/u/flashscratch/d/dpolioud/RNAseq-singleCell-Fetal.Rdata"
  , compression_level = 9
)
################################################################################

# ### Finding differentially expressed genes (cluster biomarkers)
#
# print("### Finding differentially expressed genes (cluster biomarkers)")
#
# # Seurat can help you find markers that define clusters via differential
# # expression. By default, it identifes positive and negative markers of a single
# # cluster (specified in ident.1), compared to all other cells.
# # **FindAllMarkers()** automates this process for all clusters, but you can also
# # test groups of clusters vs. each other, or against all cells.
#
# # The min.pct argument requires a gene to be detected at a minimum percentage in
# # either of the two groups of cells, and the thresh.test argument requires a
# # gene to be differentially expressed (on average) by some amount between the
# # two groups. You can set both of these to 0, but with a dramatic increase in
# # time - since this will test a large number of genes that are unlikely to be
# # highly discriminatory. As another option to speed up these computations,
# # max.cells.per.ident can be set. This will downsample each identity class to
# # have no more cells than whatever this is set to. While there is generally
# # going to be a loss in power, the speed increases can be signficiant and the
# # most highly differentially expressed genes will likely still rise to the top.
#
# # # find all markers of cluster 1
# ldf <- lapply(unique(centSO@ident), function(i) {
#   df <- FindMarkers(centSO, ident.1 = i, only.pos = TRUE
#     , test.use = "negbinom", latent.vars = c("nUMI", "librarylab", "individual")
#     , min.pct = 0.25, thresh.use = 0.25)
#   df$cluster <- i
#   return(df)
# })
# clusterDeDF <- do.call("rbind", ldf)
#
#
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# # clusterDeDF <- FindAllMarkers(centSO, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# # clusterDeDF <- FindAllMarkers(centSO, only.pos = TRUE, test.use = "negbinom"
# #   , latent.vars = c("nUMI", "librarylab", "individual")
# #   , min.pct = 0.25, thresh.use = 0.25)
# # clusterDeDF %>% group_by(cluster) %>% top_n(2, avg_diff)
# clusterDe20DF <- data.frame(clusterDeDF %>% group_by(cluster) %>% top_n(20, avg_diff))
#
# write.table(x = clusterDeDF
#   , file = paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All.txt")
#   , sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(x = clusterDe20DF
#   , file = paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All_Top20.txt")
#   , sep = "\t", quote = FALSE, row.names = FALSE)
#
# save(centSO, noCentExM, metDF, clusterDeDF, file = paste0(outData, "seuratO.Robj"))
#
# # Seurat has four tests for differential expression which can be set with the
# # test.use parameter: ROC test ("roc"), t-test ("t"), LRT test based on
# # zero-inflated data ("bimod", default), LRT test based on tobit-censoring
# # models ("tobit") The ROC test returns the 'classification power' for any
# # individual marker (ranging from 0 - random, to 1 - perfect).
#
# # cluster1.markers <- FindMarkers(centSO, ident.1 = 0, thresh.use = 0.25
# # , test.use = "roc", only.pos = T)
#
# # There are several tools for visualizing marker expression. **VlnPlot()**
# # generates a violin plot which shows the probability density at different
# # expression levels of the gene for each cluster. As seen below, good marker
# # genes will show strongly in a single cluster. It can also be useful to look at
# # gene/gene correlation with **GenePlot()** which returns a plot similar to a
# # 'FACS' plot with cells colored by cluster identity. Also, the
# # **FeaturePlot()** function is useful for viewing the expression of the gene in
# # the context of all the cells and helps validate the specificity of the marker
# # or the quality of the clustering.
#
# # VlnPlot(centSO, c("STMN2","LIMCH1"))
# #
# # #you can plot raw UMI counts as well
# # VlnPlot(centSO, c("STMN2", "LIMCH1"),use.raw = T,y.log = T)
# #
# # FeaturePlot(centSO, c("STMN2", "LIMCH1","SATB2","GPM6A","CSRP2","PLXNA4")
# # , cols.use = c("grey","blue"))
#
# # Heatmaps can also be a good way to examine heterogeneity within/between
# # clusters. The **DoHeatmap()** function will generate a heatmap for given cells
# # and genes. In this case, we are plotting the top 20 markers (or all markers if
# # less than 20) for each cluster.
#
# # Top 10 markers for each cluster
# clusterDeDF %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
#
# # Heatmap of top 10 markers
# # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
# # DoHeatmap(centSO, genes.use = top10$gene, order.by.ident = TRUE
# #   , slim.col.label = TRUE, remove.key = FALSE)
#
# # Centered and scaled
# # Subset to genes, merge to keep duplicate genes (from more than 1 PC)
# ggDF <- merge(top10, centSO@scale.data, by.x = "gene", by.y = "row.names"
#   , all.x = TRUE)
# # For gene ordering by cluster
# idx <- match(rev(top10$gene), ggDF$gene)
# ggDF <- ggDF[idx, ]
# row.names(ggDF) <- paste0(length(ggDF$gene):1, "_", ggDF$gene)
# ggDF <- ggDF[ ,-c(1:6)]
# # Order by clustering
# idx <- match(colnames(ggDF), names(sort(centSO@ident)))
# ggDF <- ggDF[ ,idx]
# # Format for ggplot2
# ggDF <- as.matrix(ggDF)
# ggDF <- melt(ggDF)
# # Add clusters
# idx <- match(ggDF$Var2, names(centSO@ident))
# ggDF$CLUSTERS <- centSO@ident[idx]
# # Set sample order by clustering
# ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# # Set expression limits
# ggDF$value[ggDF$value > 3] <- 3
# ggDF$value[ggDF$value < 0] <- 0
# # # Split by clusters with > or < 1000 cells
# # cl <- names(table(centSO@ident))[table(centSO@ident) < 1000]
# # id <- names(centSO@ident)[centSO@ident %in% cl]
# # gg2DF <- ggDF[ggDF$Var2 %in% id, ]
# # ggplot
# ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile() +
#   facet_grid(~CLUSTERS, scales = "free") +
#   # facet_grid(~CLUSTERS, space = "free", scales = "free") +
#   # scale_fill_gradient2(name = "Normalized\nexpression"
#   #   , high = "#d7191c", low = "white") +
#   scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) +
#   scale_fill_distiller(name = "Normalized\nexpression", type = "div"
#     , palette = 5, direction = -1) +
#   theme_bw() +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks = element_blank()) +
#   theme(text = element_text(size = 12)) +
#   theme(axis.text.y = element_text(size = 12)) +
#   ylab("Genes") +
#   xlab("Cells") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nTop 10 marker genes per cluster"
#     , "\nMean centered, variance scaled, normalized expression"
#     , "\n"))
# ggsave(paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers_CenterScale.png")
#   , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)
#
# # png(paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers_CenterScale.png")
# #   , width = 12, height = 12, units = "in", res = 300)
# # grid.arrange(p1, p2, ncol = 2, top = paste0(graphCodeTitle
# #   , "\n"
# #   , "\nTop 10 marker genes per cluster"
# #   , "\nMean centered, variance scaled, normalized expression"
# #   , "\n"))
# # dev.off()
#
# # No centered and scaled expression
# # Subset to genes, merge to keep duplicate genes (from more than 1 PC)
# ggDF <- merge(top10, noCentExM, by.x = "gene", by.y = "row.names"
#   , all.x = TRUE)
# # For gene ordering by cluster
# idx <- match(rev(top10$gene), ggDF$gene)
# ggDF <- ggDF[idx, ]
# row.names(ggDF) <- paste0(length(ggDF$gene):1, "_", ggDF$gene)
# ggDF <- ggDF[ ,-c(1:6)]
# # Order by clustering
# idx <- match(colnames(ggDF), names(sort(centSO@ident)))
# ggDF <- ggDF[ ,idx]
# # Format for ggplot2
# ggDF <- as.matrix(ggDF)
# ggDF <- melt(ggDF)
# # Add clusters
# idx <- match(ggDF$Var2, names(centSO@ident))
# ggDF$CLUSTERS <- centSO@ident[idx]
# # Set sample order by clustering
# ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# # Set expression limits
# # Set expression limits
# ggDF$value[ggDF$value > 3] <- 3
# ggDF$value[ggDF$value < 0] <- 0
# # ggplot
# ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile() +
#   facet_grid(~CLUSTERS, scales = "free") +
#   # facet_grid(~CLUSTERS, space = "free", scales = "free") +
#   # scale_fill_gradient2(name = "Normalized\nexpression"
#   #   , high = "#d7191c", low = "white") +
#   scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) +
#   scale_fill_distiller(name = "Normalized\nexpression", type = "div"
#     , palette = 5, direction = -1) +
#   theme_bw() +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks = element_blank()) +
#   theme(text = element_text(size = 12)) +
#   theme(axis.text.y = element_text(size = 12)) +
#   ylab("Genes") +
#   xlab("Cells") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nTop 10 marker genes"
#     , "\nNormalized expression"
#     , "\n"))
# ggsave(paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers.png")
#   , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)
#
# # # Violin plots of top 10 markers
# # pdf(paste0(outGraph, "ViolinPlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   VlnPlot(centSO, top10cluster$gene, size.use = 0.5)
# # })
# # dev.off()
# #
# # # Feature plot of top 10 markers
# # pdf(paste0(outGraph, "FeaturePlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   FeaturePlot(centSO, top10cluster$gene, cols.use = c("grey","blue")
# #     , pt.size = 0.7)
# # })
# # dev.off()
# ################################################################################

#
# ###Assigning cell type identity to clusters
# # Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:
#
# #   Cluster ID | Markers       | Cell Type
# # -----------|---------------|----------
# #   0          | IL7R          | CD4 T cells
# # 1          | CD14, LYZ     | CD14+ Monocytes
# # 2          | MS4A1         | B cells
# # 3          | CD8A          | CD8 T cells
# # 4          | FCGR3A, MS4A7 | FCGR3A+ Monocytes
# # 5          | GNLY, NKG7    | NK cells
# # 6          | FCER1A, CST3  | Dendritic Cells
# # 7          | PPBP          | Megakaryocytes
#
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
# new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
# centSO@ident <- plyr::mapvalues(centSO@ident, from = current.cluster.ids, to = new.cluster.ids)
# TSNEPlot(centSO, do.label = T, pt.size = 0.5)
#
# ###Further subdivisions within cell types If you perturb some of our parameter
# #choices above (for example, setting  `resolution=0.8` or changing the number of
# #PCs), you might see the CD4 T cells subdivide into two groups. You can explore
# #this subdivision to find markers separating the two T cell subsets. However,
# #before reclustering (which will overwrite `object@ident`), we can stash our
# #renamed identities to be easily recovered later.
#
# #First lets stash our identities for later
# centSO <- StashIdent(centSO, save.name = "ClusterNames_0.6")
#
# # Note that if you set save.snn=T above, you don't need to recalculate the SNN, and can simply put : centSO=FindClusters(centSO,resolution = 0.8)
# centSO <- FindClusters(centSO, dims.use = 1:10, resolution = 0.8, print.output = F)
#
# # Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
# plot1 <- TSNEPlot(centSO, do.return = T, no.legend = TRUE, do.label = T)
# plot2 <- TSNEPlot(centSO, do.return = T, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = T)
# plot_grid(list(plot1, plot2), cols = 2)
#
# #Find discriminating markers
# tcell.markers <- FindMarkers(centSO, 0, 1)
#
# # Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can
# # see that CCR7 is upregulated in C0, strongly indicating that we can
# # differentiate memory from naive CD4 cells.
#
# # cols.use demarcates the color palette from low to high expression
# FeaturePlot(centSO, c("LIMCH1", "MEF2C"), cols.use = c("green", "blue"))
#
# # The memory/naive split is bit weak, and we would probably benefit from looking
# # at more cells to see if this becomes more convincing (stay tuned!). In the
# # meantime, we can restore our old cluster identities for downstream processing.
#
# centSO <- SetAllIdent(centSO, id = "ClusterNames_0.6")
