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

# Set variable to gene of interest

## Input
# Log normalized, regressed nUMI and percent mito: noCentSO, centSO
load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")
load("../DS-005-006/analysis/Seurat_Cluster/Seurat_Cluster_exon_FtMm250_seuratO.Robj")
load("../analysis/Seurat_Cluster_DS-2-3-4-5-6/Seurat_Cluster_DS-2-3-4-5-6_exon_FtMm250_seuratO.Robj")

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_FeaturePlot_MarkersforSingleCell.R"
graphTitle <- "Drop-seq SeqRuns 2,3,4,5,6"
outGraphPfx <- "../analysis/graphs/Seurat_FeaturePlot_MarkersforSingleCell_DS-2-3-4-5-6"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 14))
################################################################################

### Functions

# Calculate mean expression of group of genes for each cell using seurat scaled
# expression values
Mean_Expression <- function(genes, seuratO) {
  genesExDF <- seuratO@scale.data[which(row.names(seuratO@scale.data) %in% genes), ]
  # genesExDF <- as.matrix(seuratO@data)[which(row.names(as.matrix(seuratO@data)) %in% genes), ]
  # Only calculate column means if there are multiple genes
  print("Length genes:")
  print(length(genes))
  if (is.matrix(genesExDF)) {
    mnExDF <- colMeans(genesExDF)
  } else {
    mnExDF <- genesExDF
  }
  # Add to ggplot data frame
  ggDF$EXPRESSION <- mnExDF[match(row.names(ggDF), names(mnExDF))]
  return(ggDF)
}

# Transform data to desired limits for ggplot2
Set_Limits <- function(ggDF) {
  ggDF$EXPRESSION[ggDF$EXPRESSION < 0] <- 0
  ggDF$EXPRESSION[ggDF$EXPRESSION > 2] <- 2
  return(ggDF)  
}

# Color tSNE plot by expression from Mean_Expression()
Feature_Plot <- function(ggDF, grouping) {
  ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = 0.5) +
    scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
      , high = "#e31a1c", limits = c(0, 2)) +
    ggtitle(paste0(graphCodeTitle
      , "\n\n", grouping
      , "\ntSNE plot, each point is a cell"
      , "\nColor indicates genes of interest mean expression"
      , "\n"))
  return(ggFp)
}

# tSNE plot colored by Seurat clustering
TSNE_Plot <- function(seuratO, grouping) {
  # tSNE graph colored by cluster
  ggTsne <- TSNEPlot(seuratO, do.label = TRUE, pt.size = 0.5, do.return = TRUE)
  ggTsne <- ggTsne + ggtitle(paste0(graphCodeTitle
    , "\n\n", graphTitle
    , "\n", grouping
    , "\ntSNE plot, each point is a cell"
    , "\nColor indicates cluster assignment"
    , "\nClusters annotated manually by expression of known marker genes"
    , "\n"))
  return(ggTsne)
}
################################################################################

### tSNE graph colored by cluster and feature plot

# ## Assign annotated cluster names to clusters
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
# new.cluster.ids <- c(
#   "Excitatory Upper Layer Neuron 1"
#   , "Excitatory Neuron"
#   , "Excitatory Upper Layer Neuron 2"
#   , "Excitatory Deep Layer Neuron"
#   , "Intermediate Progenitors"
#   , "Interneuron"
#   , "Mitotic Progenitors"
#   , "oRG"
#   , "Oligodendrocyte Precursor"
#   , "Endothelial")
# seuratO@ident <- plyr::mapvalues(seuratO@ident, from = current.cluster.ids
#   , to = new.cluster.ids)
# # Stash numerical cluster identities if want to use later
# seuratO <- StashIdent(seuratO, save.name = "Cluster_Numbers")

# Collect tSNE values for ggplot
ggDF <- centSO@tsne.rot

## Subset to marker genes of interest for Luis' excel file
# Cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
kmDFL <- split(kmDF, kmDF$Grouping)

# Loop through and plot each group of genes
pdf(paste0(outGraphPfx, ".pdf"), width = 13, height = 6)
lapply(names(kmDFL), function(grouping) {
  print(grouping)
  genes <- kmDF$Gene.Symbol[kmDF$Grouping == grouping]
  ggDF <- Mean_Expression(genes, noCentSO)
  ggDF <- Set_Limits(ggDF)
  ggFp <- Feature_Plot(ggDF, grouping)
  ggTsne <- TSNE_Plot(centSO, grouping)
  MultiPlotList(list(ggTsne, ggFp), cols = 2)
})
dev.off()

## Additional genes of interest not in marker table
# genes <- c("CALB1", "CALB2", "PVALB", "SST")
genes <- c("SLC12A5")
genes[genes %in% row.names(seuratO@scale.data)]

# Loop through and plot each group of genes
pdf(paste0(outGraphPfx, "_SLC12A5.pdf"), width = 13, height = 6)
lapply(genes, function(gene) {
  print(gene)
  ggDF <- Mean_Expression(gene, seuratO)
  str(ggDF)
  ggDF <- Set_Limits(ggDF)
  ggFp <- Feature_Plot(ggDF, gene)
  ggTsne <- TSNE_Plot(seuratO, gene)
  MultiPlotList(list(ggTsne, ggFp), cols = 2)
})
dev.off()

pdf(paste0(outGraphPfx, "_ViolinPlot_SLC12A5.pdf"), width = 13, height = 6)
p <- VlnPlot(seuratO, genes, do.ret = TRUE)
p[[1]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
################################################################################