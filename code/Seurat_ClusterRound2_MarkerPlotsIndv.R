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
require(tidyverse)
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

Main_Function <- function(){
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to5")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to6")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to8")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to10")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to15")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to40")
  # Plot_Genes_Of_Interest_Violin(cluster_col = "clusters_pc1to5_res_0.7")
  Plot_Genes_Of_Interest_Boxplot(cluster_col = "clusters_pc1to5_res_0.7")
  Plot_Genes_Of_Interest_Boxplot(cluster_col = "clusters_pc1to6_res_0.7")
  Plot_Genes_Of_Interest_Boxplot(cluster_col = "clusters_pc1to8_res_0.7")
  Plot_Genes_Of_Interest_Boxplot(cluster_col = "clusters_pc1to10_res_0.6")

}
################################################################################

### Format

# Luis markers - cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
################################################################################

### Luis marker genes

FeaturePlot_Graph_CentScale <- function(tsneDF, title, limLow, limHigh
  , alpha = 0.5, size = 0.02) {
  ggFp <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = size, alpha = alpha, stroke = 0.5) +
    scale_color_distiller(name = "Normalized\nexpression\nz-score", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    ggtitle(title)
  return(ggFp)
}

Plot_Luis_Markers_tSNE <- function(tsne_slot_name){

  print("Plot_Luis_Markers_tSNE")

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
################################################################################

Plot_Genes_Of_Interest_Violin <- function(cluster_col){

  print("Plot_Genes_Of_Interest_Violin")

  km_DFL <- split(kmDF, kmDF$Grouping)
  out_violin_dir <- paste0(dirname(outGraph)
    , "/MarkerPlots_", gsub(" ", "", centSO@project.name)
    , "_ViolinPlot_IndividualGene/")
  dir.create(out_violin_dir)

  cluster_ids <- factor(centSO@meta.data[[cluster_col]])
  names(cluster_ids) <- rownames(centSO@meta.data)
  centSO@ident <- cluster_ids

  genes <- kmDF$Gene.Symbol[
    kmDF$Grouping %in% c("GABAergic Neuron", "GABAergic interneuron")]

  Gene_Expression_By_Cluster_ViolinPlot(
    genes = genes
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    # , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression"
      , "\n", cluster_col
      , "\n"))
  ggsave(paste0(
      out_violin_dir, "KnownMarks_ViolinPlot_Normalized", cluster_col, ".png")
    , width = 13, height = 26)
  # print(outGraph)

}
################################################################################

plot_expression_by_cluster_boxplot <- function(
  genes
  , expr_m
  , cell_cluster_key
  , cluster_annot_key
  , cluster_order = NULL
  ){
  print("plot_expression_by_cluster_boxplot")

  # make cell_cluster_key a tibble if it's not
  if (any(class(cell_cluster_key) != "tibble")) {
    cell_cluster_key_tb <- cell_cluster_key %>%
      enframe(name = "cell_id", value = "cluster") %>%
      as_tibble
  }

  # gather data
  gg_tb <- expr_m %>%
    as_tibble(rownames = "gene") %>%
    filter(gene %in% genes) %>%
    gather(., key = "cell_id", value = "expression", -gene) %>%
    left_join(., cell_cluster_key_tb)
     # %>%
    # # order clusters
    # mutate(cluster = fct_relevel(cluster, cluster_annot_key)) %>%
    # # order genes
    # mutate(genes = fct_relevel(gene, unique(genes)))

  # plot
  ggplot(gg_tb, aes(x = cluster, y = expression, fill = cluster)) +
    facet_wrap(~gene, ncol = 4, scales = "free") +
    geom_boxplot(outlier.size = 0.1) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      , legend.position = "none")
}

Plot_Genes_Of_Interest_Boxplot <- function(cluster_col){

  print("Plot_Genes_Of_Interest_Boxplot")

  out_boxplot_dir <- paste0(dirname(outGraph)
    , "/MarkerPlots_", gsub(" ", "", centSO@project.name)
    , "_BoxPlot_IndividualGene/")
  dir.create(out_boxplot_dir)

  cluster_ids <- factor(centSO@meta.data[[cluster_col]])
  names(cluster_ids) <- rownames(centSO@meta.data)
  centSO@ident <- cluster_ids

  genes <- as.vector(kmDF$Gene.Symbol[
    kmDF$Grouping %in% c("GABAergic Neuron", "GABAergic interneuron")])

  plot_expression_by_cluster_boxplot(
    genes = genes
    , expr_m = noCentExM
    , cell_cluster_key = cluster_ids
    , cluster_annot_key = sort(unique(cluster_ids))
  ) +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nMarker gene expression by cluster"
    , "\nNormalized expression"
    # , "\n", cluster_col
    , "\n"))
  ggsave(paste0(
      out_boxplot_dir, "KnownMarks_BoxPlot_Normalized", cluster_col, ".png")
    , width = 13, height = 13)
}
################################################################################

### Run

Main_Function()
################################################################################
