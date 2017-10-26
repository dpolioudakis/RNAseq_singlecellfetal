# Damon Polioudakis
# 2017-03-07
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
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb
rm(fetb)

# Cell cycle markers from Macosko 2015 Table S2
mkDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_FeaturePlot_MacoskoCellCycle_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Seurat_FeaturePlot_MacoskoCellCycle_DS002_003_exon_FtMm250"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size = 12)))
theme_update(plot.title = element_text(size = 14))
################################################################################

### Functions

# Calculate mean expression of group of genes for each cell using seurat scaled
# expression values
Mean_Expression <- function(genes) {
  genesExDF <- seuratO@scale.data[which(row.names(seuratO@scale.data) %in% genes), ]
  # Only calculate column means if there are multiple genes
  print("Length genes:")
  print(length(genes))
  print("Genes matched in Seurat object:")
  print(table(row.names(seuratO@scale.data) %in% genes))
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
  ggDF$EXPRESSION[ggDF$EXPRESSION > 4] <- 4
  return(ggDF)  
}

# Color tSNE plot by expression from Mean_Expression()
Feature_Plot <- function(ggDF, grouping) {
  ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = 0.5) +
    scale_colour_gradient(name = "Normalized\nExpression", low = "blue"
      , high = "red", limits = c(0, 4)) +
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
    , "\n\n", grouping
    , "\ntSNE plot, each point is a cell"
    , "\nColor indicates cluster assignment"
    , "\nClusters annotated manually by expression of known marker genes"
    , "\n"))
  return(ggTsne)
}

SaveToPNG <- function(...) {
  d = dev.copy(png,...)
  dev.off(d)
}
################################################################################

### tSNE graph colored by cluster and feature plot

## Assign annotated cluster names to clusters
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
seuratO@ident <- plyr::mapvalues(seuratO@ident, from = current.cluster.ids
  , to = new.cluster.ids)
# Stash numerical cluster identities if want to use later
seuratO <- StashIdent(seuratO, save.name = "Cluster_Numbers")

# Collect tSNE values for ggplot
ggDF <- seuratO@tsne.rot

## Subset to genes of interest
# Cleanup marker data frame
mkDF <- mkDF[! mkDF$Gene.Symbol == "", ]
mkDF <- mkDF[! is.na(mkDF$Grouping), ]
mkDF$Grouping <- factor(mkDF$Grouping, levels = unique(mkDF$Grouping))
mkDFL <- split(mkDF, mkDF$Grouping)

# Loop through and plot each group of genes
pdf(paste0(outGraphPfx, ".pdf"), width = 13, height = 6)
lapply(names(mkDF), function(grouping) {
  print(grouping)
  genes <- mkDF[ ,grouping]
  # Remove spaces before gene names
  genes <- gsub(" *", "", genes)
  print(genes)
  ggDF <- Mean_Expression(genes)
  ggDF <- Set_Limits(ggDF)
  ggFp <- Feature_Plot(ggDF, grouping)
  ggTsne <- TSNE_Plot(seuratO, grouping)
  MultiPlotList(list(ggTsne, ggFp), cols = 2)
})
dev.off()
################################################################################