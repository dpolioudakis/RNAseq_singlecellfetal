# Damon Polioudakis
# 2017-02-21
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

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)

## Inputs

# Seurat
in_Seurat <- list.files("../analysis/analyzed_data/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/", full.names = TRUE)
load(in_Seurat[as.numeric(args[[1]])])
# load(in_Seurat[9])
centSO <- so
rm(so)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- paste0("MarkerPlots.R\n", centSO@project.name)
outGraph <- paste0( "../analysis/graphs/Seurat_ClusterRound2/MarkerPlots/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/MarkerPlots_"
, gsub(" ", "", centSO@project.name), "_")

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

### Functions

################################################################################

### Format

# Luis markers - cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
################################################################################

### Luis marker genes

Plot_Luis_Markers <- function(){

  # Feature plot individual expression
  # Cluster round 1 normalized, mean centered scaled
  # Collect tSNE values for ggplot
  km_DFL <- split(kmDF, kmDF$Grouping)
  out_tsne_dir <- paste0(dirname(outGraph)
    , "/MarkerPlots_", gsub(" ", "", centSO@project.name)
    , "_FeaturePlot_IndividualGene/")
  dir.create(out_tsne_dir)
  lapply(names(km_DFL), function(grouping){
    kmDF <- km_DFL[[grouping]]
    ggL <- FeaturePlot(
      genes = kmDF$Gene.Symbol
      , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
      , seuratO = centSO
      , exM = centSO@scale.data
      , limLow = -1.5
      , limHigh = 1.5
      , geneGrouping = NULL
      , centScale = TRUE)
    Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_heights = 0.2
      , title = paste0(graphCodeTitle
        , "\n\nExpression of known marker genes"
        , "\nCluster round 1 normalized expression, mean centered and variance scaled"
        , "\n")
      )
    ggsave(paste0(
        out_tsne_dir
        , "Round1NormalizedCenteredScaled_", gsub(" ", "_", grouping), ".png")
      , width = 14, height = length(ggL)*1+4, limitsize = FALSE, dpi = 200)
  })

}
################################################################################

### Run

Main_Function <- function(){
  Plot_Luis_Markers()
}

Main_Function()
################################################################################
