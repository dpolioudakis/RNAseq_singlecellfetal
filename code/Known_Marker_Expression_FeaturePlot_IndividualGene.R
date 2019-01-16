# Damon Polioudakis
# 2018-06-21
# Plot known marker lists as heatmaps, violin plots, feature plots

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
require(viridis)
source("Function_Library.R")

## Command arg to select marker gene group
args <- commandArgs(trailingOnly = TRUE)
print(args)

## Inputs

# Seurat
# PC 1-40
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)


## Variables
graphCodeTitle <- "Known_Marker_Expression.R"
outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Known_Marker_Expression_"
# outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Known_Marker_Expression_"
# outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC_PC1to40/Known_Marker_Expression_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Format

# Cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
################################################################################

### Feature plots

# Feature plot individual expression
# Normalized, no mean centering scaling
# Collect tSNE values for ggplot
km_DFL <- split(kmDF, kmDF$Grouping)
out_tsne_dir <- paste0(dirname(outGraph)
  , "/KnownMarks_FeaturePlot_IndividualGene/")
dir.create(out_tsne_dir)
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
kmDF <- km_DFL[[as.numeric(args[[1]])]]
grouping <- names(km_DFL)[[as.numeric(args[[1]])]]
ggL <- FeaturePlot(
  genes = kmDF$Gene.Symbol
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1
  , limHigh = 3
  , geneGrouping = NULL
  , centScale = FALSE)
Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_heights = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nNormalized expression"
    , "\n")
  )
ggsave(paste0(out_tsne_dir, "Normalized_", grouping, ".png")
  , width = 14, height = length(ggL)*1+4, limitsize = FALSE)

# Feature plot individual expression
# Normalized, mean centered scaled
# Collect tSNE values for ggplot
km_DFL <- split(kmDF, kmDF$Grouping)
out_tsne_dir <- paste0(dirname(outGraph)
  , "/KnownMarks_FeaturePlot_IndividualGene/")
dir.create(out_tsne_dir)
kmDF <- km_DFL[[as.numeric(args[[1]])]]
grouping <- names(km_DFL)[[as.numeric(args[[1]])]]
ggL <- FeaturePlot(
  genes = kmDF$Gene.Symbol
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE)
ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
ggL[1:2] <- lapply(ggL[1:2], function(gg){
    gg + guides(color = guide_legend(
      ncol = 2, title = "Cluster", override.aes = list(size = 3)))
})
Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_heights = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nNormalized centered scaled expression"
    , "\n")
  )
ggsave(paste0(out_tsne_dir, "NormalizedCenteredScaled_", grouping, ".png")
  , width = 14, height = length(ggL)*1+4, limitsize = FALSE)
################################################################################
