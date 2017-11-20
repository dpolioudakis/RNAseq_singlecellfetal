# Damon Polioudakis
# 2017-02-21
# Test Seurat PCA method and clustering parameters

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
require(gridExtra)
require(ggplot2)
require(cowplot)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# # Log normalized, regressed nUMI and percent mito
# # seuratO
# load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")
# # Log normalized, regressed nUMI and percent mito, mean centered and scaled
# # fetb
# load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

# Seurat
# PC 1-40
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Marker gene lists
deDF <- read.table(
  "../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Luis metaMat results
mmapDF <- read.csv("../source/metaMat/Overlapped-Genes.csv", header = TRUE)

## Variables
graphCodeTitle <- "Human_Specific_Expression.R"
outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Human_Specific_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_Expression_"

## Output Directories
outGraphDir <- dirname(outGraph)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

Subset_To_Specific_Clusters <- function(deDF, cluster, okayClusters, fcHigh, fcLow) {
  # Clusters gene cannot be DE in
  clsNo <- c(0:17)[! c(0:17) %in% okayClusters]
  # Gene is > X FC in cluster
  genes1 <- deDF$GENE[deDF$LOG_FC > fcHigh & deDF$CLUSTER == cluster]
  # Genes in clusters genes cannot be DE in > 0.3
  genes2 <- deDF$GENE[deDF$LOG_FC > fcLow & deDF$CLUSTER %in% clsNo]
  # Check
  print(table(genes1 %in% genes2))
  # Remove genes in clusters genes cannot be DE in
  genes1 <- genes1[! genes1 %in% genes2]
  # Filter DE DF
  utdeDF <- deDF[deDF$GENE %in% genes1, ]
  utdeDF <- utdeDF[utdeDF$CLUSTER == cluster, ]
  return(utdeDF)
}

Combine_DE_and_Expression <- function(deDF, exDF) {
  # Column for setting order of genes
  utdeDF$ORDER <- seq(1, nrow(utdeDF))
  # Merge with expression data frame
  ggDF <- merge(utdeDF[c("GENE", "CLUSTER", "ORDER")], exDF
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
  # Set order
  ggDF <- ggDF[order(-ggDF$ORDER), ]
  # Remove order variable now set
  ggDF <- ggDF[ ,! colnames(ggDF) == "ORDER"]
}
################################################################################

### metaMat cluster DE oRG human specific genes

genes <- mmapDF[mmapDF$X == "Human-specific", "X7"]
genes <- unlist(strsplit(as.character(genes), split = "\\|"))

################################################################################

### Human specific genes DE in oRG cluster

geneGroupDF <- data.frame(GENE = genes, GROUP = "")

ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = noCentExM
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = 0
  , upperLimit = 3
  , geneOrder = NULL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific genes DE in oRG (cluster 7)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
    , "\n")
  )
ggsave(paste0(outGraph, "HumanSpecific_Heatmap_Normalized.png")
  , width = 16, height = 10)

# Heatmap
# Normalized, mean centered and scaled
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = NULL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific genes DE in oRG (cluster 7)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n"))
ggsave(paste0(outGraph, "HumanSpecific_Heatmap_NormalizedCenteredScaled.png")
  , width = 16, height = 10)
################################################################################

### Uniquely expressed in oRG

# 18
length(genes)
# Number of DE genes 7768
nrow(deDF)
# Intersect DE gene lists and TFs, co-factors, chromatin remodelers list: 408
df <- deDF[deDF$GENE %in% genes, ]

# Subset to expressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia
ldf <- lapply(c(7), function(cluster) {
  Subset_To_Specific_Clusters(deDF = df, cluster = cluster
    , okayClusters = c(7, 8, 10, 11, 13, 15, 16), fcHigh = 0.2, fcLow = 0.1)
})
ssDeDF <- do.call("rbind", ldf)

# Heatmap
# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = ssDeDF$GENE, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = NULL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters: > 0.2 log fold change in cluster; < 0.1 for other clusters"
    , "\n")
  )
ggsave(paste0(
  outGraph, "DeUniqueRG_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 10, limitsize = FALSE)

# Heatmap
# Normalized
geneGroupDF <- data.frame(GENE = ssDeDF$GENE, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1
  , upperLimit = 3
  , geneOrder = NULL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters: > 0.2 log fold change in cluster; < 0.1 for other clusters"
    , "\n")
)
ggsave(paste0(
  outGraph, "DeUniqueRG_Heatmap_Normalized.png")
  , width = 12, height = 10, limitsize = FALSE)


## Feature plots

# Collect tSNE values for ggplot
tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

# Normalized
ggL <- FeaturePlot(
  genes = ssDeDF$GENE
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1
  , limHigh = 2
  , geneGrouping = NULL
  , centScale = FALSE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nNormalized expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeUniqueRG_FeaturePlot_Normalized.png")
  , width = 20, height = 35, limitsize = FALSE)

# Normalized centered scaled
ggL <- FeaturePlot(
  genes = ssDeDF$GENE
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nNormalized centered scaled expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeUniqueRG_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 20, height = 35, limitsize = FALSE)
################################################################################
