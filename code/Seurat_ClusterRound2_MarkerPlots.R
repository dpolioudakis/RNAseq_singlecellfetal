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
args <- 10
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
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv", header = TRUE
  , fill = TRUE)

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccTirosh <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")

# Molyneaux markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

# Lake 2017 cluster DE tables
lake_DE_DF <- read.csv(
  "../source/Lake_2018_TS3_Cluster_DE_Genes.csv"
  , skip = 4, header = TRUE, fill = TRUE
)

# Lake 2017 expression matrix
lake_ex_DF <- read.table(
  "../data/lake_2017/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", header = TRUE
)

# Luis metaMat results
mmapDF <- read.csv("../source/Gene_Lists/Overlapped-Genes.csv", header = TRUE)

# Nowakowski cluster DE table
nowakowski_DE_DF <- read.csv(
  "../data/nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)

## Variables
date <- format(Sys.Date(), "%Y%m%d")
graphCodeTitle <- paste0("MarkerPlots.R\n", centSO@project.name)
outGraph <- paste0( "../analysis/graphs/Seurat_ClusterRound2/MarkerPlots/"
  , date, "/MarkerPlots_"
  , gsub(" ", "", centSO@project.name), "_")

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
################################################################################

### Functions

Main_Function <- function(){
  Plot_Luis_Markers()
  Plot_Cell_Cycle_Genes()
  Plot_Lake2017_Enriched()
  Plot_Nowakowski_Enriched()
}

Plot_Luis_Markers <- function(){
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to5")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to8")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to10")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to15")
  Plot_Luis_Markers_tSNE(tsne_slot_name = "tsne_pc1to40")
  Plot_Luis_Markers_Violin()
  Plot_Luis_Markers_Heatmap()
}

Format_Cell_Cycle_Genes_Table <- function(){
  print("Format_Cell_Cycle_Genes_Table")
  ccDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(ccDF) <- c("Grouping", "Gene.Symbol")
  ccDF$Gene.Symbol <- gsub(" *", "", ccDF$Gene.Symbol)
  ccDF <- ccDF[! ccDF$Gene.Symbol == "", ]
  ccDF <- ccDF[! is.na(ccDF$Grouping), ]
  ccDF$Grouping <- gsub(" *", "", ccDF$Grouping)
  ccDF$Grouping <- factor(ccDF$Grouping, levels = unique(ccDF$Grouping))
  return(ccDF)
}
################################################################################

### Format

# Luis markers - cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))

## Lake
lake_DE_DF <- lake_DE_DF[ ,1:7]
# Subset to clusters in frontal cortex dataset
lake_DE_DF <- lake_DE_DF[
  lake_DE_DF$Cluster %in% unique(gsub("_.*", "", colnames(lake_ex_DF)))
  , ]
# DE filter higher
lake_DE_DF <- lake_DE_DF[
  lake_DE_DF$"Average.Difference..log.fold.change." > 0.25
  , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]

## Nowakowski
# DE filter higher
nowakowski_DE_DF <- nowakowski_DE_DF[
  nowakowski_DE_DF$avg_diff > 0.75, c("gene", "avg_diff", "cluster")]
################################################################################

### Luis marker genes

Plot_Luis_Markers_tSNE <- function(tsne_slot_name = tsne_slot_name){

  print("Plot_Luis_Markers_tSNE")
  print(tsne_slot_name)

  # Select tSNE to use
  centSO@dr$tsne <- centSO@dr[[tsne_slot_name]]

  # Feature plot
  # Normalized, no mean centering scaling
  ggL <- FeaturePlot(
    genes = kmDF$Gene.Symbol
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = 0
    , limHigh = 3
    , geneGrouping = kmDF$Grouping
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nNormalized expression"
    , "\n", tsne_slot_name
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_Normalized_", tsne_slot_name, ".png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = kmDF$Gene.Symbol
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = kmDF$Grouping
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nCluster round 2 normalized expression, mean centered and variance scaled"
    , "\n", tsne_slot_name
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph
      , "KnownMarks_FeaturePlot_Round2NormCenterScale_"
      , tsne_slot_name, ".png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Cluster round 1 normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = kmDF$Gene.Symbol
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = kmDF$Grouping
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nCluster round 1 normalized expression, mean centered and variance scaled"
    , "\n", tsne_slot_name
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(
      outGraph, "KnownMarks_FeaturePlot_Round1NormCenterScale_"
      , tsne_slot_name, ".png")
    , width = 12, height = length(ggL)*1)

}

Plot_Luis_Markers_Violin <- function(){

  print("Plot_Luis_Markers_Violin")

  ## Violin plots of genes by cluster

  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nCluster round 2 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph
    , "KnownMarks_ViolinPlot_Round2NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Cluster round 1 normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol
    , exprM = rd1CentExM
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nCluster round 1 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph, "KnownMarks_ViolinPlot_Round1NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "KnownMarks_ViolinPlot_Normalized.png")
    , width = 13, height = 26)
}

Plot_Luis_Markers_Heatmap <- function(){

  print("Plot_Luis_Markers_Heatmap")

  # Heatmap
  # Normalized, no mean centering scaling
  geneGroupDF <- data.frame(Gene = kmDF$Gene.Symbol, Group = kmDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = noCentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -0.5
    , upperLimit = 2
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of published marker genes by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "KnownMarks_Heatmap_Normalized.png")
    , width = 16, height = 80, limitsize = FALSE)

  # Clustering round 2 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = kmDF$Gene.Symbol, Group = kmDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of published marker genes by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 2 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "KnownMarks_Heatmap_Round2NormalizedCenteredScaled.png")
    , width = 16, height = 80, limitsize = FALSE)

  # Clustering round 1 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = kmDF$Gene.Symbol, Group = kmDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = rd1CentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of published marker genes by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 1 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "KnownMarks_Heatmap_Round1NormalizedCenteredScaled.png")
    , width = 16, height = 80, limitsize = FALSE)

}
################################################################################

### Cell cycle genes

Plot_Cell_Cycle_Genes <- function(){

  print("Plot_Cell_Cycle_Genes")

  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

  # Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  ccDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(ccDF) <- c("Grouping", "Gene.Symbol")
  ccDF$Gene.Symbol <- gsub(" *", "", ccDF$Gene.Symbol)
  ccDF <- ccDF[! ccDF$Gene.Symbol == "", ]
  ccDF <- ccDF[! is.na(ccDF$Grouping), ]
  ccDF$Grouping <- gsub(" *", "", ccDF$Grouping)
  ccDF$Grouping <- factor(ccDF$Grouping, levels = unique(ccDF$Grouping))
  ccDFL <- split(ccDF, ccDF$Grouping)

  # We can segregate this list into markers of G2/M phase and markers of S
  # phase
  ccDFL[["S.Tirosh"]] <- data.frame(
    Grouping = "S.Tirosh", Gene.Symbol = ccTirosh[1:43])
  ccDFL[["G2M.Trisosh"]] <- data.frame(
    Grouping = "G2M.Tirosh", Gene.Symbol = ccTirosh[44:98])

  # Feature plot
  # Normalized, no mean centering scaling
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  ggL <- FeaturePlot(
    genes = ccDF$Gene.Symbol
    , tsneDF = ggDF
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = ccDF$Grouping
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of cell cycle genes from Macosko"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  ggsave(paste0(outGraph, "CellCycle_FeaturePlot_Normalized.png")
    , width = 14, height = length(ggL)*1)

  # Feature plot
  # Cluster round 2 normalized, mean centered scaled
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  ggL <- FeaturePlot(
    genes = ccDF$Gene.Symbol
    , tsneDF = ggDF
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = ccDF$Grouping
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of cell cycle genes from Macosko"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  ggsave(paste0(outGraph
    , "CellCycle_FeaturePlot_Round2NormalizedCenteredScaled.png")
    , width = 14, height = length(ggL)*1)

  # Feature plot
  # Cluster round 1 normalized, mean centered scaled
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  ggL <- FeaturePlot(
    genes = ccDF$Gene.Symbol
    , tsneDF = ggDF
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = ccDF$Grouping
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of cell cycle genes from Macosko"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  ggsave(paste0(outGraph
    , "CellCycle_FeaturePlot_Round1NormalizedCenteredScaled.png")
    , width = 14, height = length(ggL)*1)

  # Heatmap
  # Normalized, no mean centering scaling
  # Round 1 normalized, mean centered and scaled
  geneGroupDF <- data.frame(Gene = ccDF$Gene.Symbol, Group = ccDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = noCentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -0.5
    , upperLimit = 2
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of cell cycle genes (Macosko) by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "CellCycle_Heatmap_Normalized.png")
    , width = 16, height = 44)

  # Round 2 normalized mean centered scaled
  geneGroupDF <- data.frame(Gene = ccDF$Gene.Symbol, Group = ccDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of cell cycle genes (Macosko) by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 2 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "CellCycle_Heatmap_Round2NormalizedCenteredScaled.png")
    , width = 16, height = 44)

  # Round 1 normalized, mean centered and scaled
  geneGroupDF <- data.frame(Gene = ccDF$Gene.Symbol, Group = ccDF$Grouping)
  geneGroupDF <- geneGroupDF[! duplicated(geneGroupDF), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = rd1CentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of cell cycle genes (Macosko) by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 1 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "CellCycle_Heatmap_Round1NormalizedCenteredScaled.png")
    , width = 16, height = 44)

}
################################################################################

### Lake cluster enriched genes

Plot_Lake2017_Enriched <- function(){

  print("Plot_Lake2017_Enriched")

  ## Feature plots

  # Feature plot
  # Normalized, no mean centering scaling
  ggL <- FeaturePlot(
    genes = lake_DE_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = lake_DE_DF$Cluster
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Lake cluster enriched"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "LakeEnriched_FeaturePlot_Normalized.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = lake_DE_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = lake_DE_DF$Cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Lake cluster enriched"
    , "\nCluster round 2 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph
    , "LakeEnriched_FeaturePlot_Round2NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Cluster round 1 normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = lake_DE_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = lake_DE_DF$Cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Lake cluster enriched"
    , "\nCluster round 1 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "LakeEnriched_FeaturePlot_Round1NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)


  ## Violin plots of genes by cluster

  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = lake_DE_DF$Gene
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = lake_DE_DF$Cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Lake cluster enriched"
      , "\nCluster round 2 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph
    , "LakeEnriched_ViolinPlot_Round2NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Cluster round 1 normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = lake_DE_DF$Gene
    , exprM = rd1CentExM
    , clusterIDs = centSO@ident
    , grouping = lake_DE_DF$Cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Lake cluster enriched"
      , "\nCluster round 1 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph, "LakeEnriched_ViolinPlot_Round1NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = lake_DE_DF$Gene
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = lake_DE_DF$Cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Lake cluster enriched"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "LakeEnriched_ViolinPlot_Normalized.png")
    , width = 13, height = 26)


  # Heatmap
  # Normalized, no mean centering scaling
  geneGroupDF <- data.frame(Gene = lake_DE_DF$Gene, Group = lake_DE_DF$Cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = noCentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -0.5
    , upperLimit = 2
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Lake cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "LakeEnriched_Heatmap_Normalized.png")
    , width = 16, height = 40, limitsize = FALSE)

  # Clustering round 2 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = lake_DE_DF$Gene, Group = lake_DE_DF$Cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Lake cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 2 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "LakeEnriched_Heatmap_Round2NormalizedCenteredScaled.png")
    , width = 16, height = 40, limitsize = FALSE)

  # Clustering round 1 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = lake_DE_DF$Gene, Group = lake_DE_DF$Cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = rd1CentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Lake cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 1 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "LakeEnriched_Heatmap_Round1NormalizedCenteredScaled.png")
    , width = 16, height = 40, limitsize = FALSE)

}
################################################################################

### Nowakowski cluster enriched genes

Plot_Nowakowski_Enriched <- function(){

  print("Plot_Nowakowski_Enriched")

  # Feature plot
  # Normalized, no mean centering scaling
  ggL <- FeaturePlot(
    genes = nowakowski_DE_DF$gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = nowakowski_DE_DF$cluster
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Nowakowski cluster enriched"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "NowakowskiEnriched_FeaturePlot_Normalized.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = nowakowski_DE_DF$gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = nowakowski_DE_DF$cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Nowakowski cluster enriched"
    , "\nCluster round 2 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph
    , "NowakowskiEnriched_FeaturePlot_Round2NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Cluster round 1 normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = nowakowski_DE_DF$gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = nowakowski_DE_DF$cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Nowakowski cluster enriched"
    , "\nCluster round 1 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "NowakowskiEnriched_FeaturePlot_Round1NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)


  ## Violin plots of genes by cluster

  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = nowakowski_DE_DF$gene
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = nowakowski_DE_DF$cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Nowakowski cluster enriched"
      , "\nCluster round 2 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph
    , "NowakowskiEnriched_ViolinPlot_Round2NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Cluster round 1 normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = nowakowski_DE_DF$gene
    , exprM = rd1CentExM
    , clusterIDs = centSO@ident
    , grouping = nowakowski_DE_DF$cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Nowakowski cluster enriched"
      , "\nCluster round 1 normalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph, "NowakowskiEnriched_ViolinPlot_Round1NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = nowakowski_DE_DF$gene
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = nowakowski_DE_DF$cluster
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nExpression of Nowakowski cluster enriched"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "NowakowskiEnriched_ViolinPlot_Normalized.png")
    , width = 13, height = 26)


  # Heatmap
  # Normalized, no mean centering scaling
  geneGroupDF <- data.frame(Gene = nowakowski_DE_DF$gene, Group = nowakowski_DE_DF$cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = noCentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -0.5
    , upperLimit = 2
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Nowakowski cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "NowakowskiEnriched_Heatmap_Normalized.png")
    , width = 16, height = 40, limitsize = FALSE)

  # Clustering round 2 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = nowakowski_DE_DF$gene, Group = nowakowski_DE_DF$cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Nowakowski cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 2 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "NowakowskiEnriched_Heatmap_Round2NormalizedCenteredScaled.png")
    , width = 16, height = 40, limitsize = FALSE)

  # Clustering round 1 normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = nowakowski_DE_DF$gene, Group = nowakowski_DE_DF$cluster)
  geneGroup_DFL <- lapply(split(geneGroupDF, geneGroupDF$Group), function(ssGeneGroupDF){
    return(ssGeneGroupDF[1:10, ])
  })
  geneGroupDF <- do.call("rbind", geneGroup_DFL)
  geneGroupDF <- geneGroupDF[! is.na(geneGroupDF$Gene), ]
  cellID_clusterID <- centSO@ident
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = rd1CentExM
    , cellID_clusterID = cellID_clusterID
    , lowerLimit = -1.5
    , upperLimit = 1.5
    , clusterOrder = sort(unique(centSO@ident))
  )
  gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nExpression of Nowakowski cluster top 10 enriched by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nClustering round 1 normalized mean centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "NowakowskiEnriched_Heatmap_Round1NormalizedCenteredScaled.png")
    , width = 16, height = 40, limitsize = FALSE)
}
################################################################################

### Run

Main_Function()
################################################################################
