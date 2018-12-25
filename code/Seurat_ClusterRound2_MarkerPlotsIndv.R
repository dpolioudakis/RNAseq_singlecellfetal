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
source("GGplot_Theme.R")

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
# args <- 4
print(args)

## Inputs

# Seurat
in_Seurat <- list.files("../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20181222/", full.names = TRUE)
print(in_Seurat[as.numeric(args[[1]])])
load(in_Seurat[as.numeric(args[[1]])])
# load(in_Seurat[9])
centSO <- so
rm(so)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

## Variables
date <- format(Sys.Date(), "%Y%m%d")
graphCodeTitle <- paste0("MarkerPlots.R\n", centSO@project.name)

## Outputs
outGraph <- paste0( "../analysis/graphs/Seurat_ClusterRound2/MarkerPlots/"
  , date, "/MarkerPlots_"
  , gsub(" ", "", centSO@project.name), "_")

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
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

Plot_Luis_Markers <- function(tsne_slot_name){

  # Feature plot individual expression
  # Cluster round 1 normalized, mean centered scaled
  # Collect tSNE values for ggplot
  km_DFL <- split(kmDF, kmDF$Grouping)
  out_tsne_dir <- paste0(dirname(outGraph)
    , "/MarkerPlots_", gsub(" ", "", centSO@project.name)
    , "_FeaturePlot_IndividualGene/")
  dir.create(out_tsne_dir)
  centSO@dr$tsne <- centSO@dr[[tsne_slot_name]]
  lapply(names(km_DFL), function(grouping){
    kmDF <- km_DFL[[grouping]]
    ggL <- FeaturePlot(
      genes = kmDF$Gene.Symbol
      , tsneDF = as.data.frame(centSO@dr[[tsne_slot_name]]@cell.embeddings)
      , seuratO = centSO
      , exM = centSO@scale.data
      , limLow = -1.5
      , limHigh = 1.5
      , geneGrouping = NULL
      , centScale = TRUE
      , size = 1)
    Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_heights = 0.25
      , title = paste0(graphCodeTitle
        , "\n\nExpression of known marker genes"
        , "\n", tsne_slot_name
        , "\nCluster round 2 normalized expression, mean centered and variance scaled"
        , "\n")
      )
    ggsave(paste0(
        out_tsne_dir
        , "Round2NormalizedCenteredScaled_", tsne_slot_name
        , "_", gsub(" ", "_", grouping), ".png")
      , width = 14, height = length(ggL)*1+4, limitsize = FALSE, dpi = 200)
  })
  # print(outGraph)

}

Plot_Luis_Markers_Wrapper <- function(){
  Plot_Luis_Markers("tsne_pc1to5")
  Plot_Luis_Markers("tsne_pc1to8")
  Plot_Luis_Markers("tsne_pc1to10")
  Plot_Luis_Markers("tsne_pc1to15")
  Plot_Luis_Markers("tsne_pc1to40")
}
################################################################################

### Run

Main_Function <- function(){
  Plot_Luis_Markers_Wrapper()
}

Main_Function()
################################################################################
