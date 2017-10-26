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
require(viridis)
require(cowplot)
require(monocle)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# Seurat
# PC 1-40
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO_TEST.Robj")

# Monocle
load("../analysis/Monocle/Monocle_monocleO.Robj")

# Luis known markers
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE)

# Molyneaux layer markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

## Variables
graphCodeTitle <- "Neuron_Differentiation.R"
outGraph <- "../analysis/graphs/Neuron_Differentiation/Neuron_Differentiation_DS2-11_PC1-40_"
outTable <- "../analysis/tables/Neuron_Differentiation/Neuron_Differentiation_DS2-11_PC1-40_"

## Output Directories
outGraphDir <- dirname(outGraph)
dir.create(outGraphDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

# Wrapper function for Heatmap_By_Cluster
# Uses plot_grid to combine genes of interest and all other genes expression
# heatmaps to deal with cell number scaling
Heatmaps_By_Clusters_Plots <- function(geneGroupDF, exprM, seuratO
  , clusters1, clusters2, clusters3, lowerLimit, upperLimit, geneOrder) {
  
  gg1 <- Heatmap_By_Cluster(
    geneGroupDF = geneGroupDF
    , exprM = exprM[ ,ids]
    , seuratO = seuratO
    , clusters = 0:17
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneGroupDF$GENE)
  gg1 <- gg1 + theme(legend.position = "none")
  
  gg2 <- Heatmap_By_Cluster(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , seuratO = seuratO
    , clusters = 0:1
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneGroupDF$GENE)
  gg2 <- gg2 + theme(legend.position = "none"
    , axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank())
  
  gg3 <- Heatmap_By_Cluster(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , seuratO = seuratO
    , clusters = c(2:10)
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneGroupDF$GENE)
  gg3 <- gg3 + theme(legend.position = "none"
    , axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank())
  
  gg4 <- Heatmap_By_Cluster(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , seuratO = seuratO
    , clusters = c(11:17)
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneGroupDF$GENE)
  gg4 <- gg4 + theme(legend.position = "none"
    , axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank())
  
  ggL <- list(gg1, gg2, gg3, gg4)
  return(ggL)
}

Subset_Cells_By_Genes_Expression <- function(exprM, keep = NULL, exclude = NULL) {
  # Keep cells with > 1 expression in all these genes
  if (! is.null(keep)) {
    df <- exprM[row.names(exprM) %in% keep, , drop = FALSE] > 1
    keepIdx <- apply(df, 2, function(x) all(x))
    print(table(keepIdx))  
  } else {
    keepIdx <- NULL
  }
  # Filter cells with < 1 expression in all these genes
  if (! is.null(exclude)) {
    df2 <- exprM[row.names(exprM) %in% exclude, , drop = FALSE] < 1
    excludeIdx <- apply(df2, 2, function(x) all(x))
    print(table(excludeIdx))  
  } else {
    excludeIdx <- NULL
  }
  # Subset to cell IDs that match keep and exclude indexes
  if (! is.null(keep) & is.null(exclude)) {
    ids <- colnames(exprM)[keepIdx]
  }
  if (! is.null(exclude) & is.null(keep)) {
    ids <- colnames(exprM)[keepIdx]
  }
  if (! is.null(keep) & ! is.null(exclude)) {
    ids <- colnames(exprM)[keepIdx & excludeIdx]
  }
  # print(exprM[row.names(exprM) %in% c(keep, exclude), ids])
  return(ids)
}

# Barplot of DE values of gene list
DE_Barplot <- function(genes, deDF, title, ylab) {
  
  ggDF <- deDF[row.names(deDF) %in% genes, ]
  colnames(ggDF) <- c("LOG_FC", "PVALUE")
  ggDF$GENE <- row.names(ggDF)
  ggDF$PVALUE <- round(ggDF$PVALUE, 2)
  gg <- ggplot(ggDF, aes(x = GENE, y = LOG_FC)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = PVALUE, x = GENE, y = LOG_FC)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Gene") +
    ylab(ylab) +
    ggtitle(title)
  return(gg)
}

# Bar plots of number of genes in each Seurat cluster based on subset
Plot_Number_In_Subset <- function(exprM, keep = NULL, exclude = NULL) {
  
  # Subset Cell IDs
  ids <- Subset_Cells_By_Genes_Expression(
    noCentExM
    , keep = keep
    , exclude = exclude)
  
  # Number of subset cell IDs in each cluster
  nClusterDF <- data.frame(table(centSO@ident[names(centSO@ident) %in% ids]))
  
  # Barplot of number of cell IDs in each cluster
  gg <- ggplot(nClusterDF, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq, y = Freq), angle = 90) +
    xlab("Seurat cluster") +
    ylab("Number of cells") +
    ggtitle(paste0("Keep: ", paste(keep, collapse = " ")
      , "\nFilter: ", paste(exclude, collapse = " ")))
  
  return(gg)
}
################################################################################

### Neuronal differentiation markers

# EOMES + PAX6 - mutually exclusive in mouse
# PAX6 - RGs
# EOMES + SOX2 - same logic as PAX6
# Neurogenin2 + EOMES
# TUJ1 / TUBB3 - most early neuron marker

## Genes of interest
genes <- c(
  # RG
  "PAX6", "SOX2", "VIM", "HES1"
  # IPC
  , "SSTR2", "EOMES", "PPP1R17", "ELAVL4", "PENK"
  # Neuron
  , "STMN2", "MAP2", "TUBB3")
genes[genes %in% row.names(centSO@scale.data)]
print(genes)

## Heatmaps

# Normalized, no mean centering scaling
p1 <- Heatmap_By_Gene_Expression(genes = genes, exprM = centSO@scale.data
  , limLow = -1.5, limHigh = 1.5)
p1 <- p1 + ylab("Genes") +
  # scale_fill_distiller(name = "Normalized\nexpression\nz-score") +
  ylab("Genes") +
  xlab("Cells ordered by PAX6 expression") +
  ggtitle("Normalized mean centered and scaled expression")

# Normalized
p2 <- Heatmap_By_Gene_Expression(genes = genes, exprM = noCentExM
  , limLow = -1, limHigh = 3)
p2 <- p2 + ylab("Genes") +
  xlab("Cells ordered by PAX6 expression") +
  ggtitle("Normalized expression")

# plot_grid
pg <- plot_grid(p1, p2, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of neuronal differentiation marker genes"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by PAX6 expression"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.25, 1))
ggsave(paste0(outGraph, "Markers_Heatmap.png"), width = 12, height = 7)


## Scatter plots

# scale_color_viridis(name = "EOMES\nexpression", limits = c(limLow, limHigh)) +
#   xlab("PAX6 expression") +
#   ylab("STMN2 expression") +
#   ggtitle("Normalized expression")
p1 <- ScatterPlot_Expression(c("PAX6", "STMN2", "EOMES"), noCentExM)
p2 <- ScatterPlot_Expression(c("PAX6", "MAP2", "EOMES"), noCentExM)
p3 <- ScatterPlot_Expression(c("SOX2", "STMN2", "EOMES"), noCentExM)
p4 <- ScatterPlot_Expression(c("SOX2", "MAP2", "EOMES"), noCentExM)
# plot_grid
pg <- plot_grid(p1, p2, p3, p4, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of neuronal differentiation marker genes"
  , "\nNormalized expression"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "Markers_ScatterPlots.png")
  , width = 9, height = 10)


## Violin plots of genes by cluster
# Normalized, no mean centering scaling
ggDF <- merge(genes, noCentExM
  , by.x = 1, by.y = "row.names", all.x = TRUE)
row.names(ggDF) <- ggDF$x
ggDF <- ggDF[ ,-1]
ggDF <- t(ggDF)
ggDF <- as.data.frame(ggDF)
# Violin plots of genes by cluster
ggDF$CLUSTER <- centSO@ident
ggDF <- melt(ggDF)
png(paste0(outGraph, "ClusterViolinPlot.png")
  , width = 13, height = 6, units = "in", res = 300)
ggplot(ggDF, aes(x = CLUSTER, y = value)) +
  geom_violin(aes(fill = CLUSTER)) +
  geom_jitter(size = 0.05, height = 0, alpha = 0.1) +
  facet_wrap(~variable, scales = "free") +
  theme(legend.position = "none") +
  ylab("Normalized expression") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nMarker gene expression by cluster"
    , "\nNormalized expression - read depth, regressed nUMI, brain, lablibrary"
    , "\n"))
dev.off()


## tSNE feature plot

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "SOX2", "VIM", "HES1", "STMN2", "MAP2")
  , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4", "PENK")
)
# Check cells that PAX6 and STMN2 > 1
noCentExM[row.names(noCentExM) %in%
    c("PAX6", "SOX2", "STMN2", "MAP2", "SSTR2", "PPP1R17", "ELAVL4", "PENK")
  , colnames(noCentExM) %in% ids]
# Add column specifying if cells match ids
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0("\nColor indicates:"
    , "\nPAX6, SOX2, VIM, HES1 > 1 normalized expression"
    , "\nSSTR2, EOMES, PPP1R17, ELAVL4, PENK < 1 normalized expression"
    , "\nSTMN2, MAP2 > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "SOX2", "STMN2", "MAP2"
    , "SSTR2", "EOMES", "PPP1R17", "ELAVL4", "PENK")
)
# Add column specifying if cells match ids
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nPAX6 > 1 normalized expression"
    , "\nSSTR2, EOMES, PPP1R17, ELAVL4, PENK < 1 normalized expression"
    , "\nSTMN2, MAP2 > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
  , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4", "PENK")
)
# Add column specifying if cells match ids
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg3 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nPAX6 > 1 normalized expression"
    , "\nSSTR2, EOMES, PPP1R17, ELAVL4, PENK > 1 normalized expression"
    , "\nSTMN2, MAP2 > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "SOX2", "STMN2", "MAP2")
  , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4", "PENK")
)
# Add column specifying if cells PAX6 > 1 and STMN2 > 1
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg4 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nPAX6, SOX2 > 1 normalized expression"
    , "\nSSTR2, EOMES, PPP1R17, ELAVL4, PENK < 1 normalized expression"
    , "\nSTMN2, MAP2 > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "EOMES")
  # , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4")
)
# Add column specifying if cells PAX6 > 1 and STMN2 > 1
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg5 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nPAX6 > 1 normalized expression"
    , "\nEOMES > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES")
  # , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4")
)
# Add column specifying if cells PAX6 > 1 and STMN2 > 1
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg6 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nEOMES > 1 normalized expression"
    , "\n"))

# plot_grid
pg <- plot_grid(gg1, gg2, gg3, gg4, gg5, gg6, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of RG, IPC, and early neuronal genes"
  , "\nNormalized expression"
  , "\ntSNE plots, each point is a cell colored by passing expression filter"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "tSNE_FeaturePlot.png")
  , width = 10, height = 15)

## Feature Plots of PAX6, STMN2, MAP2 cells and EOMES, STMN2, MAP2 cells

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
  # , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4")
)
# Add column specifying if cells PAX6 > 1 and STMN2 > 1
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nPAX6, STMN2, MAP2 > 1 normalized expression"
    , "\n"))

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES", "STMN2", "MAP2")
  # , exclude = c("SSTR2", "EOMES", "PPP1R17", "ELAVL4")
)
# Add column specifying if cells PAX6 > 1 and STMN2 > 1
ggDF$PASS_EXPRESSION_FILTER <- "False"
ggDF$PASS_EXPRESSION_FILTER[row.names(ggDF) %in% ids] <- "True"
# ggplot
gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PASS_EXPRESSION_FILTER)) +
  geom_point(size = 0.05) +
  geom_point(data = ggDF[ggDF$PASS_EXPRESSION_FILTER == "True", ]
    , aes(x = tSNE_1, y = tSNE_2), size = 0.2) +
  scale_colour_manual(name = "Pass expression\nfilter:"
    , values = c("#a6cee3", "red")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0(
    "\nColor indicates:"
    , "\nEOMES, STMN2, MAP2 > 1 normalized expression"
    , "\n"))

# plot_grid
pg <- plot_grid(gg1, gg2, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of PAX6, STMN2, MAP2 cells and EOMES, STMN2, MAP2 cells"
  , "\nNormalized expression"
  , "\ntSNE plots, each point is a cell colored by passing expression filter"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "tSNE_FeaturePlot_PAX6_EOMES.png")
  , width = 10, height = 6)

## Feature Plot of genes: PAX6, EOMES

# Genes to plot
genes <- c("PAX6", "EOMES", "STMN2", "MAP2")

# FeaturePlot outputs list of ggplots
ggL <- FeaturePlot_CentScale(genes, tsneDF, centSO, -1.5, 1.5)
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nExpression of genes of interest"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nRegress out nUMI, donor, library lab"
  , "\nNormalized expression, mean centered and variance scaled"
  , "\ntSNE PC 1-40"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(outGraph, "FeaturePlot_PAX6_EOMES_NormCentScale.png")
  , width = 14, height = 3+length(ggL)*1.25)

# FeaturePlot outputs list of ggplots
ggL <- FeaturePlot(genes, tsneDF, seuratO = centSO, exM = noCentExM, -1, 3)
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nExpression of genes of interest"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nRegress out nUMI, donor, library lab"
  , "\nNormalized expression"
  , "\ntSNE PC 1-40"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(outGraph, "FeaturePlot_PAX6_EOMES.png")
  , width = 14, height = 3+length(ggL)*1.25)
################################################################################

### Heatmap of cells

ids <- colnames(noCentExM)[
  # RG
  noCentExM[row.names(noCentExM) == "PAX6", ] > 1 &
    noCentExM[row.names(noCentExM) == "SOX2", ] > 1 &
    noCentExM[row.names(noCentExM) == "VIM", ] > 1 &
    noCentExM[row.names(noCentExM) == "HES1", ] > 1 &
    # IPC
    noCentExM[row.names(noCentExM) == "SSTR2", ] < 1 &
    noCentExM[row.names(noCentExM) == "EOMES", ] < 1 &
    noCentExM[row.names(noCentExM) == "PPP1R17", ] < 1 &
    noCentExM[row.names(noCentExM) == "ELAVL4", ] < 1 &
    noCentExM[row.names(noCentExM) == "PENK", ] < 1 &
    # Neuron
    noCentExM[row.names(noCentExM) == "STMN2", ] > 1 &
    noCentExM[row.names(noCentExM) == "MAP2", ] > 1
  # noCentExM[row.names(noCentExM) == "TUBB3", ] > 1
  ]

v1 <- rowMeans(noCentExM[ ,ids])
topGenes <- sort(v1, decreasing = TRUE)[1:500]
geneGroupDF <- data.frame(GENE = factor(names(topGenes), levels = names(topGenes))
  , GROUP = NA)

gg1 <- Heatmap_By_Cluster(
  geneGroupDF = geneGroupDF
  , exprM = noCentExM[ ,ids]
  , seuratO = centSO
  , clusters = 0:17
  , lowerLimit = 0
  , upperLimit = 3
  , geneOrder = geneGroupDF$GENE)
gg1 <- gg1 + theme(legend.position = "none"
  , axis.title.x = element_blank()
      , axis.title.y = element_blank()
      , axis.text.y = element_blank()
      , axis.ticks.y = element_blank())

gg2 <- Heatmap_By_Cluster(
  geneGroupDF = geneGroupDF
  , exprM = noCentExM
  , seuratO = centSO
  , clusters = 0:1
  , lowerLimit = 0
  , upperLimit = 3
  , geneOrder = geneGroupDF$GENE)
gg2 <- gg2 + theme(legend.position = "none")

gg3 <- Heatmap_By_Cluster(
  geneGroupDF = geneGroupDF
  , exprM = noCentExM
  , seuratO = centSO
  , clusters = c(2, 7:10)
  , lowerLimit = 0
  , upperLimit = 3
  , geneOrder = geneGroupDF$GENE)
gg3 <- gg3 + theme(legend.position = "none"
  , axis.title.x = element_blank()
  , axis.title.y = element_blank()
  , axis.text.y = element_blank()
  , axis.ticks.y = element_blank())

# ggL <- Heatmaps_By_Cluster_Combined(
#   geneGroupDF = geneGroupDF
#   , exprM = noCentExM
#   , seuratO = centSO
#   , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17
#   , lowerLimit = 0, upperLimit = 3
#   , geneOrder = geneGroupDF$GENE)
# ggL <- lapply(ggL, function(gg) {
#   gg <- gg + theme(
#     axis.title.x = element_blank()
#     , axis.title.y = element_blank()
#     , axis.text.y = element_blank()
#     , axis.ticks.y = element_blank())
#   return(gg)
# })

# ggL <- append(list(gg), ggL)
# pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
pg <- plot_grid(gg2, gg1, gg3, ncol = 3, align = 'h', axis = 'b')

ggsave(paste0(outGraph, "RG_Neuron_Cells_Heatmap.png"), width = 13)
################################################################################

### Heatmaps of DE genes expression

# Cell IDs
ids <- colnames(noCentExM)[
  # RG
  noCentExM[row.names(noCentExM) == "PAX6", ] > 1 &
    noCentExM[row.names(noCentExM) == "SOX2", ] > 1 &
    noCentExM[row.names(noCentExM) == "VIM", ] > 1 &
    noCentExM[row.names(noCentExM) == "HES1", ] > 1 &
    # IPC
    noCentExM[row.names(noCentExM) == "SSTR2", ] < 1 &
    noCentExM[row.names(noCentExM) == "EOMES", ] < 1 &
    noCentExM[row.names(noCentExM) == "PPP1R17", ] < 1 &
    noCentExM[row.names(noCentExM) == "ELAVL4", ] < 1 &
    noCentExM[row.names(noCentExM) == "PENK", ] < 1 &
    # Neuron
    noCentExM[row.names(noCentExM) == "STMN2", ] > 1 &
    noCentExM[row.names(noCentExM) == "MAP2", ] > 1
  # noCentExM[row.names(noCentExM) == "TUBB3", ] > 1
  ]

## DE

# Filter cells
df <- DE_Filters_ExpMatrix(centSO
  , minPercent = 10, foldChange = 0.2, cellID = ids)

# DE Linear model
termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$id <- FALSE
termsDF$id[colnames(centSO@data) %in% ids] <- TRUE
mod <- "y ~ id+nUMI+librarylab+individual"
deLM <- DE_Linear_Model(exDatDF = df, termsDF = termsDF, mod = mod)

# Format LM output into data frame
# Combine log2 fold changes, p-values
deDF <- data.frame(GENE = row.names(deLM$coefmat)
  , LOG_FC = deLM$coefmat[ ,2]
  , PVALUE = deLM$pvalmat[ ,2])
print(head(deDF))

# Order by log fold change
deDF <- deDF[order(-deDF$LOG_FC), ]


## Heatmaps

## >0.4 log fold change
df <- deDF[deDF$LOG_FC > 0.4, ]
geneGroupDF <- data.frame(GENE = factor(df$GENE, levels = df$GENE)
  , GROUP = NA)

# Heatmaps expression no center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = noCentExM
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = 0, upperLimit = 3)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank()
  , axis.title.y = element_blank()
  , axis.text.y = element_blank()
  , axis.ticks.y = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >0.4 log fold change"
  , "\nNormalized expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE04_Heatmap.png"), height = 10, width = 13)

# Heatmaps expression center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = -1.5, upperLimit = 1.5)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank()
  , axis.title.y = element_blank()
  , axis.text.y = element_blank()
  , axis.ticks.y = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >0.4 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE04_Heatmap_CentScale.png"), height = 10, width = 13)

## >1 log fold change
df <- deDF[deDF$LOG_FC > 1, ]
geneGroupDF <- data.frame(GENE = factor(df$GENE, levels = df$GENE)
  , GROUP = NA)

# Heatmaps expression no center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = noCentExM
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = 0, upperLimit = 3)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank()
  , axis.title.y = element_blank()
  , axis.text.y = element_blank()
  , axis.ticks.y = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >1 log fold change"
  , "\nNormalized expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE1_Heatmap.png"), height = 10, width = 13)

# Heatmaps expression center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = -1.5, upperLimit = 1.5)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank()
  , axis.title.y = element_blank()
  , axis.text.y = element_blank()
  , axis.ticks.y = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >1 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE1_Heatmap_CentScale.png"), height = 10, width = 13)

## >1.5 log fold change
df <- deDF[deDF$LOG_FC > 1.5, ]
geneGroupDF <- data.frame(GENE = factor(df$GENE, levels = df$GENE)
  , GROUP = NA)

# Heatmaps expression no center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = noCentExM
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = 0, upperLimit = 3)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >1.5 log fold change"
  , "\nNormalized expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE15_Heatmap.png"), height = 15, width = 13)

# Heatmaps expression center scale
ggL <- Heatmaps_By_Clusters_Plots(geneGroupDF = geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, geneOrder = geneGroupDF$GENE, lowerLimit = -1.5, upperLimit = 1.5)
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none"
  , axis.title.x = element_blank())})
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes >1.5 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph, "RG_Neuron_Cells_DE15_Heatmap_CentScale.png"), height = 15, width = 13)
################################################################################

### Bar plots of number of genes in each Seurat cluster based on subset

ggL <- list(
  Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("PAX6", "STMN2", "MAP2"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("EOMES", "STMN2", "MAP2"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("EOMES", "PAX6"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("EOMES"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("PAX6"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("STMN2"))
  , Plot_Number_In_Subset(exprM = noCentExM
    , keep = c("MAP2"))
)
# Plot grid
pg <- plot_grid(plotlist = ggL, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nNumbers of cells in each Seurat cluster subset by expression"
  , "\nKeep: Normalized expression > 1"
  , "\nFilter: Normalized expression < 1"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "Number_Cells_Subset_Barplot.png"), width = 9
  , height = 16)
################################################################################

### DE of PAX6, STMN2, "MAP2" cells vs "EOMES", "STMN2", "MAP2" cells

# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids1 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
)
ids2 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES", "STMN2", "MAP2")
)
# Subset expression matrix to cells to keep
idx <- colnames(centSO@data) %in% c(ids1, ids2)
exM <- as.matrix(centSO@data)[ ,idx]
termsDF <- centSO@meta.data[idx, c("CELL", "nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$CELL_FILTER <- "PAX6"
termsDF$CELL_FILTER[termsDF$CELL %in% c(ids2)] <- "EOMES"
# LM
mod <- "y ~ CELL_FILTER+nUMI+librarylab+individual"
# DE via LM
deLM <- DE_Linear_Model(exM, termsDF, mod)
# Format LM output
deDF <- data.frame(LOG_FC_PAX6vsEOMES = deLM$coefmat[ ,c("CELL_FILTERPAX6")]
  , PVALUE = deLM$pvalmat[ ,c("CELL_FILTERPAX6")])
# Save as table
write.csv(deDF, paste0(outTable, "DE_PAX6_STMN2_MAP2_vs_EOMES_SMNT2_MAP2.csv")
  , quote = FALSE)

# Heatmap of expression of DE 0.5 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.5], GROUP = "PAX6 cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.5], GROUP = "EOMES cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of PAX6, STMN2, MAP2 cells vs EOMES, STMN2, MAP2 cells"
  , "\n>0.5 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE05_PAX6_STMN2_MAP2_vs_EOMES_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 6, width = 13)

# Heatmap of expression of DE 0.4 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.4], GROUP = "PAX6 cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.4], GROUP = "EOMES cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of PAX6, STMN2, MAP2 cells vs EOMES, STMN2, MAP2 cells"
  , "\n>0.4 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE04_PAX6_STMN2_MAP2_vs_EOMES_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 8, width = 13)

# Boxplots of marker genes DE

# Excitatory Deep Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Deep Layer Cortical"]
p1 <- DE_Barplot(genes, deDF, "Excitatory Deep Layer Cortical"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Excitatory Upper Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Upper Layer Cortical"]
p2 <- DE_Barplot(genes, deDF, "Excitatory Upper Layer Cortical"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Subplate
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Subplate"]
p3 <- DE_Barplot(genes, deDF, "Subplate"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Plot grid
pg <- plot_grid(p1, p2, p3, ncol = 1) #, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of marker genes"
  , "\nPAX6, STMN2, MAP2 cells vs EOMES, STMN2, MAP2 cells"
  , "\nText indicates p-value"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "DE_PAX6_STMN2_MAP2_vs_EOMES_SMNT2_MAP2_Barplot.pdf")
  , width = 7, height = 16)
###############################################################################

### DE of PAX6, STMN2, MAP2 cells vs STMN2, MAP2, not PAX6 cells

# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids1 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
)
ids2 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("STMN2", "MAP2")
  , exclude = c("PAX6")
)
# Subset expression matrix to cells to keep
idx <- colnames(centSO@data) %in% c(ids1, ids2)
exM <- as.matrix(centSO@data)[ ,idx]
termsDF <- centSO@meta.data[idx, c("CELL", "nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$CELL_FILTER <- "PAX6 STMN2 MAP2"
termsDF$CELL_FILTER[termsDF$CELL %in% c(ids2)] <- "PAX6neg STMN2 MAP2"
# LM
mod <- "y ~ CELL_FILTER+nUMI+librarylab+individual"
# DE via LM
deLM <- DE_Linear_Model(exM, termsDF, mod)
# Format LM output
deDF <- data.frame(LOG_FC_PAX6vsPAX6neg = deLM$coefmat[ ,c("CELL_FILTERPAX6 STMN2 MAP2")]
  , PVALUE = deLM$pvalmat[ ,c("CELL_FILTERPAX6 STMN2 MAP2")])
# Save as table
write.csv(deDF, paste0(outTable, "DE_PAX6_STMN2_MAP2_vs_PAX6neg_SMNT2_MAP2.csv")
  , quote = FALSE)

# Heatmap of expression of DE 0.5 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.5], GROUP = "PAX6+ STMN2+ MAP2+ cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.5], GROUP = "PAX6- STMN2+ MAP2+ cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of PAX6+, STMN2+, MAP2+ cells vs PAX6-, STMN2+, MAP2+ cells"
  , "\n>0.5 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE05_PAX6_STMN2_MAP2_vs_PAX6neg_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 6, width = 13)

# Heatmap of expression of DE 0.4 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.4], GROUP = "PAX6+ STMN2+ MAP2+ cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.4], GROUP = "PAX6- STMN2+ MAP2+ cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of PAX6+, STMN2+, MAP2+ cells vs PAX6-, STMN2+, MAP2+ cells"
  , "\n>0.4 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE04_PAX6_STMN2_MAP2_vs_PAX6neg_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 8, width = 13)

# Boxplots of marker genes DE

# Excitatory Deep Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Deep Layer Cortical"]
p1 <- DE_Barplot(genes, deDF, "Excitatory Deep Layer Cortical"
  , ylab = "Log fold change (PAX6+ cells / PAX6- cells)")
# Excitatory Upper Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Upper Layer Cortical"]
p2 <- DE_Barplot(genes, deDF, "Excitatory Upper Layer Cortical"
  , ylab = "Log fold change (PAX6+ cells / PAX6- cells)")
# Subplate
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Subplate"]
p3 <- DE_Barplot(genes, deDF, "Subplate"
  , ylab = "Log fold change (PAX6+ cells / PAX6- cells)")
# Plot grid
pg <- plot_grid(p1, p2, p3, ncol = 1) #, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of marker genes"
  , "\nPAX6+, STMN2+, MAP2+ cells vs PAX6-, STMN2+, MAP2+ cells"
  , "\nText indicates p-value"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "DE_PAX6_STMN2_MAP2_vs_PAX6neg_SMNT2_MAP2_Barplot.pdf")
  , width = 7, height = 16)
###############################################################################

### DE of CP PAX6, STMN2, MAP2 cells vs CP EOMES, STMN2, MAP2 cells

# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
# And subset to CP
ids1 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
)
ids1 <- ids1[ids1 %in% centSO@meta.data$CELL[centSO@meta.data$REGION == "CP"]]
ids2 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES", "STMN2", "MAP2")
)
ids2 <- ids2[ids2 %in% centSO@meta.data$CELL[centSO@meta.data$REGION == "CP"]]
# Subset expression matrix to cells to keep
idx <- colnames(centSO@data) %in% c(ids1, ids2)
exM <- as.matrix(centSO@data)[ ,idx]
termsDF <- centSO@meta.data[idx, c("CELL", "nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$CELL_FILTER <- "PAX6"
termsDF$CELL_FILTER[termsDF$CELL %in% c(ids2)] <- "EOMES"
# LM
mod <- "y ~ CELL_FILTER+nUMI+librarylab+individual"
# DE via LM
deLM <- DE_Linear_Model(exM, termsDF, mod)
# Format LM output
deDF <- data.frame(LOG_FC_PAX6vsEOMES = deLM$coefmat[ ,c("CELL_FILTERPAX6")]
  , PVALUE = deLM$pvalmat[ ,c("CELL_FILTERPAX6")])
# Save as table
write.csv(deDF, paste0(outTable, "DE_CP_PAX6_STMN2_MAP2_vs_CP_EOMES_SMNT2_MAP2.csv")
  , quote = FALSE)

# Heatmap of expression of DE 0.5 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.5], GROUP = "PAX6 cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.5], GROUP = "EOMES cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of CP PAX6, STMN2, MAP2 cells vs CP EOMES, STMN2, MAP2 cells"
  , "\n>0.5 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE05_CP_PAX6_STMN2_MAP2_vs_CP_EOMES_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 6, width = 13)

# Heatmap of expression of DE 0.4 genes by cluster
geneGroupDF <- rbind(
  data.frame(GENE = row.names(deDF)[deDF$LOG_FC > 0.4], GROUP = "PAX6 cells")
  , data.frame(GENE = row.names(deDF)[deDF$LOG_FC < -0.4], GROUP = "EOMES cells")
)
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
  , clusters1 = 0:1, clusters2 = 2:10, clusters3 = 11:17, lowerLimit = -1.5, upperLimit = 1.5)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of PAX6, STMN2, MAP2 cells vs EOMES, STMN2, MAP2 cells"
  , "\n>0.4 log fold change"
  , "\nNormalized centered scaled expression"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "DE04_CP_PAX6_STMN2_MAP2_vs_CP_EOMES_SMNT2_MAP2_Heatmap_CentScale.png")
  , height = 8, width = 13)

# Boxplots of marker genes DE

# Excitatory Deep Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Deep Layer Cortical"]
p1 <- DE_Barplot(genes, deDF, "Excitatory Deep Layer Cortical"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Excitatory Upper Layer Cortical
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Excitatory Upper Layer Cortical"]
p2 <- DE_Barplot(genes, deDF, "Excitatory Upper Layer Cortical"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Subplate
genes <- kmDF$Gene.Symbol[kmDF$Grouping == "Subplate"]
p3 <- DE_Barplot(genes, deDF, "Subplate"
  , ylab = "Log fold change (PAX6 cells / EOMES cells)")
# Plot grid
pg <- plot_grid(p1, p2, p3, ncol = 1) #, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes of marker genes"
  , "\nPAX6, STMN2, MAP2 cells vs EOMES, STMN2, MAP2 cells"
  , "\nText indicates p-value"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "DE_CP_PAX6_STMN2_MAP2_vs_CP_EOMES_SMNT2_MAP2_Barplot.pdf")
  , width = 7, height = 16)
################################################################################

### Colored on Monocle trajectory

# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
ids1 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
)
ids2 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES", "STMN2", "MAP2")
)

# Add to monocle object for plotting
# PAX6
mo_filtered@phenoData@data$PAX6 <- FALSE
mo_filtered@phenoData@data$PAX6[mo_filtered@phenoData@data$CELL %in% ids1] <- TRUE
# EOMES
mo_filtered@phenoData@data$EOMES <- FALSE
mo_filtered@phenoData@data$EOMES[mo_filtered@phenoData@data$CELL %in% ids2] <- TRUE

# Plot Seurat clusters on trajectory
ggL <- list(
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "PAX6", cell_size = 0.01) +
    scale_colour_manual(name = "Cluster:", values = c("#a6cee3", "red")) +
    ggtitle("PAX6+, STMN2+, MAP2+") +
    theme(legend.position = "right")
  , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "EOMES", cell_size = 0.01) +
    scale_colour_manual(name = "Cluster:", values = c("#a6cee3", "red")) +
    ggtitle("EOMES+, STMN2+, MAP2+") +
    theme(legend.position = "right")
)
# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 7)))
})
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 2, align = 'h', axis = 't')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.2))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by cells passing expression filters"
  , "\nComponent 1 versus 2"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
# Save
ggsave(paste0(outGraph, "Monocle_PAX6_EOMES.png"), width = 13, height = 7)

## CP PAX6+ STMN2+ MAP2+ and CP EOMES+ STMN2+ MAP2+

# IDs of cells to keep / exclude (> 1 or < 1 normalized expression)
# And subset to CP
ids1 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("PAX6", "STMN2", "MAP2")
)
ids1 <- ids1[ids1 %in% centSO@meta.data$CELL[centSO@meta.data$REGION == "CP"]]
ids2 <- Subset_Cells_By_Genes_Expression(noCentExM
  , keep = c("EOMES", "STMN2", "MAP2")
)
ids2 <- ids2[ids2 %in% centSO@meta.data$CELL[centSO@meta.data$REGION == "CP"]]

# Add to monocle object for plotting
# PAX6
mo_filtered@phenoData@data$PAX6 <- FALSE
mo_filtered@phenoData@data$PAX6[mo_filtered@phenoData@data$CELL %in% ids1] <- TRUE
# EOMES
mo_filtered@phenoData@data$EOMES <- FALSE
mo_filtered@phenoData@data$EOMES[mo_filtered@phenoData@data$CELL %in% ids2] <- TRUE

# Plot Seurat clusters on trajectory
ggL <- list(
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "PAX6", cell_size = 0.01) +
    scale_colour_manual(name = "Cluster:", values = c("#a6cee3", "red")) +
    ggtitle("CP PAX6+, STMN2+, MAP2+") +
    theme(legend.position = "right")
  , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "EOMES", cell_size = 0.01) +
    scale_colour_manual(name = "Cluster:", values = c("#a6cee3", "red")) +
    ggtitle("CP EOMES+, STMN2+, MAP2+") +
    theme(legend.position = "right")
)
# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 7)))
})
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 2, align = 'h', axis = 't')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.2))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by cells passing expression filters"
  , "\nComponent 1 versus 2"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
# Save
ggsave(paste0(outGraph, "Monocle_CP_PAX6_EOMES.png"), width = 13, height = 7)