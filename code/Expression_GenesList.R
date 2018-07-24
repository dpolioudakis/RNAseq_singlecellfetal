# Damon Polioudakis
# 2017-09-13
# TF expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3 +
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

### Set variable to gene of interest
# genes <- c("ARG1", "ARG2")
# genes <- c("TTR", "SLC13A4")
# genes <- c("XIST")
# genes <- c("PLXNA2", "SEMA5A", "SEMA6A")
# genes <- c("RELN", "TBR1", "EMX2", "NECTIN1", "CALB2")
# genes <- c("LDB2", "POU3F1", "POU3F2", "OTX1", "BCL11B", "NEFH", "LMO4", "FEZF2")
# genes <- c("DIAPH3", "IGFBP4", "CRIM1", "FEZF2", "BCL11B", "PCP4", "S100A10", "TBR1")
# genes <- c("SOX2-OT")
# genes <- read.table("/u/project/geschwind/wgpembro/Circadian/CYCLOPS/Seed/Zhang.9.seed.txt")[ ,1]
# genes <- c("CNTNAP2", "NRXN1", "NRXN2", "NRXN3", "KDM5B", "ARID1B", "ASH1L")
# genes <- c("GAPDH", "HPRT1", "UBC")
# genes <- c("FOXP1")
# genes <- c("SCN5A")
genes <- c("DHX57", "MSLT8")

### Set out graphs path (will append type of graph to end, e.g.
# ../analysis/graphs/Expression_GenesList_ARG1-2_ExprHeatmap_NormCentScale.png)
outGraph <- "../analysis/graphs/Expression_GenesList/Expression_GenesList_DHX57_MSLT8_"
################################################################################

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
require(cowplot)
require(ggplot2)
source("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/code/Function_Library.R")

## Inputs

# Log normalized, regressed nUMI and percent mito
# seuratO
load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

# biomaRt gene info
bmDF <- read.csv("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Expression_GenesList.R"

## Output Directories
outGraphDir <- dirname(outGraph)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################


### Convert to gene symbol

# genes <- bmDF$hgnc_symbol[match(genes, bmDF$ensembl_gene_id)]

################################################################################

### Plots

## Violin plots of expression by cluster
Plot_Expression_By_Cluster_Violin <- function(gene){
  gg <- tryCatch(
    {
    # Expression per cluster
    ggDF <- noCentExM
    ggDF <- data.frame(EXPRESSION = ggDF[row.names(ggDF) == gene, ])
    ggDF$CLUSTER <- centSO@ident
    gg <- ggplot(ggDF, aes(x = CLUSTER, y = EXPRESSION)) +
      geom_violin(aes(fill = CLUSTER)) +
      geom_jitter(size = 0.05, height = 0, alpha = 0.1) +
      theme(legend.position = "none") +
      ylab("Normalized expression") +
      xlab("Clusters") +
      ggtitle(gene)
    return(gg)
    },
    error = function(condition){
      message(gene)
      message(condition)
      return(NA)
  })
  return(gg)
}
Plot_Expression_By_Cluster_Violin_Run <- function(){
  ggL <- lapply(genes, function(gene) {
    Plot_Expression_By_Cluster_Violin(gene)
  })
  ggL <- ggL[! is.na(ggL)]
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'l')
  # now add the title
  title <- paste0(graphCodeTitle
    , "\n"
    , "\nExpression of genes of interest by cluster"
    , "\nNormalized expression"
    , "\n")
  # rel_heights values control title margins
  title <- ggdraw() + draw_label(title)
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.1, length(ggL)))
  # Save graph
  ggsave(paste0(outGraph, "violinPlots.png"), width = 14, height = 4+length(ggL))
}
Plot_Expression_By_Cluster_Violin_Run()


## Feature plots

Plot_Expression_FeaturePlotCentScale_Run <- function(){
  # Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- FeaturePlot_CentScale(genes = genes
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , limLow = -1.5
    , limHigh = 1.5
  )
  ggL[-1] <- lapply(ggL[-1], function(gg){
    gg +
    ggplot_set_theme_publication
    # geom_point(size = 0.02, alpha = 0.5)
  })
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\ntSNE plots, each point is a cell"
    , "\nColor indicates normalized expression, mean centered, variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.3, length(ggL)))
  # Save graph
  ggsave(paste0(
    outGraph, "FeaturePlot_NormCentScale.png")
    , width = 15, height = 2.5+length(ggL)*1.25, limitsize = FALSE)
}
Plot_Expression_FeaturePlotCentScale_Run()


Plot_Expression_FeaturePlot_Run <- function(){
  # Feature plot - normalized
  # Loop through and plot each group of genes
  ggL <- FeaturePlot(genes = genes
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -1
    , limHigh = 3
  )
  ggL[-1] <- lapply(ggL[-1], function(gg){
    gg +
    ggplot_set_theme_publication +
    geom_point(size = 0.02, alpha = 0.5)
  })
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\ntSNE plots, each point is a cell"
    , "\nColor indicates normalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.3, length(ggL)))
  # Save graph
  ggsave(paste0(
    outGraph, "FeaturePlot_Norm.png")
    , width = 15, height = 2.5+length(ggL)*1.25, limitsize = FALSE)
}
Plot_Expression_FeaturePlot_Run()

## Heatmaps

Plot_Expression_HeatmapCentScale <- function(){
  print("Plot_Expression_Heatmap")
  # Cluster index
  cluster_annot_idx <- c(
    "vRG" = 9
    , "oRG" = 7
    , "PgS" = 8
    , "PgG2M" = 10
    , "IPC" = 2
    , "ExN" = 0
    , "ExM" = 1
    , "ExCal" = 4
    , "ExDp1" = 3
    , "ExDp2" = 13
    , "InSST" = 5
    , "InCALB2" = 6
    , "OPC" = 11
    , "End" = 12
    , "Per" = 14
    , "Mic" = 15
  )
  names(cluster_annot_idx) <- paste0(
    names(cluster_annot_idx), " = ", cluster_annot_idx)
  cellID_clusterID <- names(cluster_annot_idx)[
    match(centSO@ident, cluster_annot_idx)]
  names(cellID_clusterID) <- names(centSO@ident)
  # Normalized, mean centering scaling
  geneGroupDF <- data.frame(Gene = genes, Group = "")
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
    , clusters = names(cluster_annot_idx)
    , clusterOrder = names(cluster_annot_idx)
  )
  gg + ggtitle(paste0(
    graphCodeTitle
      , "\n\nExpression of genes of interest"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
  )
  ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale.png")
    , width = 12, height = 4 + length(genes)/2, limitsize = FALSE
  )
}
Plot_Expression_HeatmapCentScale()

# Normalized
Plot_Expression_Heatmap <- function(){
  print("Plot_Expression_Heatmap")
  geneGroupDF <- data.frame(Gene = genes, Group = "")
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = as.matrix(noCentExM)
    , cellID_clusterID <- centSO@ident
    , lowerLimit = 0.5
    , upperLimit = 2
  )
  gg <- gg + ggplot_set_theme_publication
  gg + ggtitle(paste0(
    graphCodeTitle
      , "\n\nExpression of genes of interest"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "ExprHeatmap.png")
    , width = 12, height = 4 + length(genes)/2, limitsize = FALSE
  )
}
Plot_Expression_Heatmap()


# # Normalized, mean centering scaling
# geneGroupDF <- data.frame(Gene = genes, Gene = "")
# # Plot
# ggL <- Heatmaps_By_Cluster_Combined(
#   geneGroupDF = geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
#   , lowerLimit = -1.5, upperLimit = 1.5
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
#   , geneOrder = TRUE
# )
# # now add the title
# title <- paste0(graphCodeTitle
#   , "\n\nExpression of genes of interest"
#   , "\nx-axis: Genes"
#   , "\ny-axis: Cells ordered by cluster"
#   , "\nNormalized expression, mean centered, variance scaled"
#   , "\n")
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# # now add the title
# title <- ggdraw() + draw_label(title)
# # rel_heights values control title margins
# pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(1/length(genes), 1))
# # Save graph
# ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale.png")
#   , width = 12, height = 4 + length(genes)/2, limitsize = FALSE)
################################################################################
