# Damon Polioudakis
# 2017-05-16
# WGCNA
################################################################################

rm(list = ls())
sessionInfo()
set.seed(27)

require(Seurat)
require(Matrix)
require(WGCNA)
require(cluster)
require(flashClust)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# disableWGCNAThreads()

# Expression and metadata data stored as Seurat object
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# noCentExM <- ssNoCentExM
# centSO <- ssCentSO

# Variables
outGraph <- "../analysis/graphs/WGCNA/DS2-11/Test_SoftPower_Threshold/WGCNA_Test_SoftPower_Threshold_"
outAnalysis <- "../analysis/analyzed_data/WGCNA/DS2-11/Test_SoftPower_Threshold/WGCNA_Test_SoftPower_Threshold_"
graphCodeTitle <- "WGCNA_Test_SoftPower_Threshold.R DS2-11"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outAnalysis), recursive = TRUE)
################################################################################

### Functions

Soft_Threshold <- function(exDF, powers) {
  exDF <- t(exDF)
  sft = pickSoftThreshold(exDF, powerVector =  powers, verbose =  5
    , corFnc = "bicor")
}

SoftThreshold_Plot <- function(sft, plot.title) {
  # Plot the results:
  par(mfrow = c(1, 2))
  cex1 = 0.9

  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3])*sft$fitIndices[ ,2]
    , xlab = "Soft Threshold (power)"
    , ylab = "Scale Free Topology module Fit,signed R^2",type = "n"
    , main = paste(plot.title))
  text(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3])*sft$fitIndices[ ,2],
    labels = powers,cex = cex1,col = "red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90,col = "red")
  abline(h = 0.80,col = "blue")
  abline(h = 0.70,col = "orange")
  abline(h = 0.60,col = "green")

  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], sft$fitIndices[ ,5]
    , xlab = "Soft Threshold (power)"
    , ylab = "Mean Connectivity"
    , type = "n"
    , main = paste("Mean connectivity"))
  text(sft$fitIndices[ ,1], sft$fitIndices[ ,5], labels = powers, cex = cex1
    , col = "red")
}

# Module identification using hybrid tree cut:
# Function to construct modules
# Args: minimum module size, deep split
Make_Modules <- function (minModSize, deepSplit, geneTree, dissTOM) {
  print("Treecut arguments:")
  # print(c("minModSize"=minModSize,"cut.height"=cutHeightMergeME, "deepSplit"=ds))
  print(c(minModSize, deepSplit))
  tree = cutreeHybrid(dendro = geneTree, pamRespectsDendro = TRUE
    , minClusterSize = minModSize
    # , cutHeight = 0.967
    , deepSplit = deepSplit, distM = as.matrix(dissTOM))
  print("Table of genes per module:")
  print(table(tree$labels))
  tree$labels
}

# Merge modules based on ME function
# Args: Modules colors, Cut height to merge ME
Merge_Modules_ME <- function (exDF, genesModuleColor, cutHeightMergeME) {
  # Call an automatic merging function
  # merged: The merged module colors
  # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
  merged <- mergeCloseModules(exprData = exDF, colors = genesModuleColor,
    cutHeight = cutHeightMergeME, trapErrors = TRUE)
  print("Table of genes per module:")
  print(table(merged$colors))
  labels2colors(merged$colors)
}
################################################################################

### Graphs to determine genes to use

## Gene total counts versus percentile
# Total coutns for each gene
v1 <- rowSums(centSO@raw.data)
ggDF <- data.frame(GENE = names(v1), COUNT_SUM = v1)
ggDF <- ggDF[order(-ggDF$COUNT_SUM), ]
ggDF$RANK <- rank(-ggDF$COUNT_SUM)
# Convert rank to percentile
ggDF$PERCENTILE <- ggDF$RANK/nrow(ggDF)*100
head(ggDF)
# ggplot
ggplot(ggDF, aes(x = PERCENTILE, y = COUNT_SUM)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1000)) +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nGenes ranked by total counts"
  ))
ggsave(paste0(outGraph, "Expression_Percentile.pdf"))
################################################################################

### Filter expression matrix

# Mean expression of each gene to use for filtering
mnEx <- sort(rowMeans(noCentExM), decreasing = TRUE)
# Check top 10% of genes sorted by mean expression
names(mnEx[1:(length(mnEx)*0.1)])

# List of expression matrices after different filters
exLDF <- list(

  # Regression normalization
  # Variable genes + 1000 genes detected threshold
  VarGenes_1000genesDetected = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% centSO@var.genes
    , centSO@meta.data$nGene > 1000]
  # Variable genes + top 10% expressed genes normalized not centered scaled
  # 479
  , NoCenterScale_VarGenes_Top10 = as.data.frame(as.matrix(noCentExM))[
    dimnames(noCentExM)[[1]] %in% centSO@var.genes &
      dimnames(noCentExM)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]

  # Variable genes + top 10% expressed genes
  , VarGenes_Top10 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% centSO@var.genes &
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # Variable genes + top 20% expressed genes
  , VarGenes_Top20 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% centSO@var.genes &
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  # Variable genes + top 30% expressed genes
  , VarGenes_Top30 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% centSO@var.genes &
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.3)]), ]

  # Top 10% expressed genes
  , Top10 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # Top 20% expressed genes
  , Top20 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  # Top 30% expressed genes
  , Top30 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.3)]), ]

  # Top 1000 expressed genes
  , Top1000 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:1000]), ]
  # Top 5000 expressed genes
  , Top5000 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:5000]), ]
  # Top 10000 expressed genes
  , Top10000 = as.data.frame(as.matrix(centSO@scale.data))[
    dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:10000]), ]
)

# Subset for testing
# idx <- sample(1:ncol(exNmVgNg1000DF), 100)
# exNmVgNg1000DF <- exNmVgNg1000DF[ ,idx]
################################################################################

### Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers <- c(seq(2, 30, 2), seq(35, 100, 5))

# Call the network topology analysis function
sftL <- lapply(exLDF, function(df) {
  Soft_Threshold(df, powers)
})

# Plot
pdf(paste0(outGraph, "Power.pdf"), height = 6, width = 12)
lapply(names(sftL), function(name){
  sft <- sftL[[name]]
  SoftThreshold_Plot(sft
    , paste0("Scale independence ", name))
})
dev.off()
################################################################################
