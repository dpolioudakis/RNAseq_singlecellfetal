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
genes <- c("PLXNA2", "SEMA5A", "SEMA6A")

### Set out graphs path (will append type of graph to end, e.g.
# ../analysis/graphs/Expression_GenesList_ARG1-2_ExprHeatmap_NormCentScale.png)
outGraph <- "../analysis/graphs/Expression_GenesList/Expression_GenesList_PLXNA2_SEMA5A_SEMA6A_"
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
load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

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

### Plots

## Violin plots of expression by cluster
ggL <- lapply(genes, function(gene) {
  # Expression per cluster
  ggDF <- noCentExM
  ggDF <- data.frame(EXPRESSION = ggDF[row.names(ggDF) == gene, ])
  ggDF$CLUSTER <- centSO@ident
  ggplot(ggDF, aes(x = CLUSTER, y = EXPRESSION)) +
    geom_violin(aes(fill = CLUSTER)) +
    geom_jitter(size = 0.05, height = 0, alpha = 0.1) +
    theme(legend.position = "none") +
    ylab("Normalized expression") +
    xlab("Clusters") +
    ggtitle(gene)
})
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


## Feature plots

# Feature plot - normalized, mean centered scaled on full dataset
# Loop through and plot each group of genes
ggL <- FeaturePlot_CentScale(genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , limLow = -1.5
  , limHigh = 1.5
)
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE plots, each point is a cell"
  , "\nColor indicates normalized expression, mean centered, variance scaled"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1
  , rel_heights = c(length(ggL)*0.25, length(ggL)))
# Save graph
ggsave(paste0(
  outGraph, "FeaturePlot_NormCentScale.png")
  , width = 14, height = 2.5+length(ggL)*1.25, limitsize = FALSE)

# Feature plot - normalized
# Loop through and plot each group of genes
ggL <- FeaturePlot(genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1
  , limHigh = 3
)
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE plots, each point is a cell"
  , "\nColor indicates normalized expression"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1
  , rel_heights = c(length(ggL)*0.25, length(ggL)))
# Save graph
ggsave(paste0(
  outGraph, "FeaturePlot_Norm.png")
  , width = 14, height = 2.5+length(ggL)*1.25, limitsize = FALSE)


## Heatmaps

# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = genes, GROUP = "")
# Plot
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
  , geneOrder = NULL
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of genes of interest"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\n")
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(1/length(genes), 1))
# Save graph
ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale.png")
  , width = 12, height = 4 + length(genes)/2, limitsize = FALSE)
################################################################################