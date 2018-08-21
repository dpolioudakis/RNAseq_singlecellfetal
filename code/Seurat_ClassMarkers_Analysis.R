# Damon Polioudakis
# 2017-08-28
# Pairwise DE of each cluster to find cluster exclusively expressed genes

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
# require(irlba)
require(gridExtra)
require(fdrtool)
require(cowplot)

source("Function_Library.R")

## Inputs
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM
# Biomart to add ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClassMarkers.R"
# Output paths
# Sub path
out_sub_path <- paste0(
  "Seurat_ClassMarkers/DS2-11/"
  ,"FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/"
  , "res054/"
  , "Seurat_ClassMarkers_"
)
outGraph <- paste0("../analysis/graphs/", out_sub_path)
outTable <- paste0("../analysis/tables/", out_sub_path)
outData <- paste0("../analysis/analyzed_data/", out_sub_path)

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 16)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

# Heatmap of mean expression per cluster
DE_Mean_Heatmap <- function(
  markers_DF, exDF, clusterIDs, ggtitle, upLim, lowLim) {

  # Subset to DE genes
  ggDF <- exDF[match(markers_DF$GENE, row.names(exDF)), ]

  # Change column names from Cell IDs to Cluster ID
  colnames(ggDF) <- clusterIDs
  # Melt
  ggDF <- melt(ggDF)
  # Mean expression per cluster
  ggDF <- aggregate(value~Var2+Var1, mean, data = ggDF)
  # Order genes by fold change
  ggDF$Var1 <- factor(ggDF$Var1, rev(levels(ggDF$Var1)))
  # Set expression limits
  ggDF$value[ggDF$value > upLim] <- upLim
  ggDF$value[ggDF$value < lowLim] <- lowLim
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    # scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(lowLim, upLim)) +
    scale_x_continuous(breaks = unique(ggDF$Var2), labels = unique(ggDF$Var2)) +
    theme_bw() +
    # theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    # theme(axis.text.y = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(ggtitle)
}
################################################################################

### Format and compile

## DE tables into one table

inTables <- list.files(dirname(outTable), full.names = TRUE)
inTables <- inTables[grepl("Markers_Class", inTables)]

# Loop through DE text files and compile into one table
dfl <- lapply(inTables, function(in_table) {
  df <- read.csv(in_table, header = TRUE)
  df <- df[! duplicated(df$Gene), c("Gene", "Class1", "Ensembl")]
  # Add ensembl IDs
  # df$Ensembl <- Convert_Mixed_GeneSym_EnsID_To_EnsID(as.character(df$Gene))
  return(df)
})
markers_DF <- do.call("rbind", dfl)
colnames(markers_DF) <- c("Gene", "Class", "Ensembl")

# Order by class
class_cluster_idx <- c(
  "Radial glia" = 7
  , "Radial glia" = 9
  , "Cycling progenitor" = 8
  , "Cycling progenitor" = 10
  , "Intermediate progenitor" = 2
  , "Excitatory Neuron" = 0
  , "Excitatory Neuron" = 1
  , "Excitatory Neuron" = 4
  , "Excitatory Neuron" = 3
  , "Excitatory Neuron" = 13
  , "Interneuron" = 5
  , "Interneuron" = 6
  , "Oligodendrocyte precursor" = 11
  , "Endothelial" = 12
  , "Pericyte" = 14
  , "Microglia" = 15
)
class_order <- unique(names(class_cluster_idx))
markers_DF$Class <- factor(markers_DF$Class, levels = class_order)
markers_DF <- markers_DF[order(markers_DF$Class), ]

# Check
head(markers_DF)
tail(markers_DF)
dim(markers_DF)
length(grep("ENSG", markers_DF$Ensembl))
head(markers_DF[order(markers_DF$Class, markers_DF$Gene), ], 20)
# Write out as tab delimited
write.csv(x = markers_DF
  , file = paste0(outTable, "Markers_Class.csv")
  , quote = FALSE, row.names = FALSE
)
################################################################################

### Plots

# Plot
geneGroupDF <- data.frame(
  Gene = rev(markers_DF$Gene)
  , Group = ""
)
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of cell type specific genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "HeatmapSetColWid_NormCenterScale.png"
  )
  , width = 7, height = 16
)
################################################################################

## Plot pvalues
Plot_DE_Pval_Hist <- function(){
  print("Plot_DE_Pval_Hist")
  deDFL <- split(markers_DF, markers_DF$Cluster)
  pgL <- lapply(names(deDFL), function(name){
    deDF <- deDFL[[name]]
    # ggplot histogram p-values
    p1 <- ggplot(deDF, aes(x = Pvalue)) +
      geom_histogram() +
      xlab("P-value") +
      ylab("Count") +
      ggtitle("P-value") +
      theme(plot.title = element_text(hjust = 0.5))
    # ggplot histogram FDR
    p2 <- ggplot(deDF, aes(x = FDR)) +
      geom_histogram() +
      xlab("FDR corrected p-value") +
      ylab("Count") +
      ggtitle("Benjamini Hochberg corrected p-value") +
      theme(plot.title = element_text(hjust = 0.5))
    # plot_grid
    pg <- plot_grid(p1, p2, ncol = 2)
    # now add the title
    title <- ggdraw() + draw_label(paste0("\nCluster: ", name))
    # rel_heights values control title margins
    plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  })
  Plot_Grid(pgL, rel_height = 0.05, ncol = 1
    , title = paste0(graphCodeTitle
      , "\n\nHistograms of p-values before and after Benjamini Hochberg correction")
  )
}
Plot_DE_Pval_Hist()
ggsave(paste0(outGraph, "Pvalue_Histogram_Cluster.pdf")
  , width = 9, height = 45)


# ## Expression heatmaps of DE genes
#
# # Centered scaled
# geneGroupDF <- data.frame(GENE = deDF$Gene, GROUP = "")
# # Heatmaps
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
#   , seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
#   , clusters1 = c(0:16), geneOrder = geneGroupDF$GENE
# )
# # Remove y-axis labels
# ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
#   , "\nMean centered, variance scaled, normalized expression"
#   , "\nCells sorted by cluster (columns)"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale_Cluster", clusterID, ".png")
#   , width = 12, height = 8)
#
# # Not centered scaled
# geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
# # Heatmaps
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = noCentExM
#   , seuratO = centSO, lowerLimit = 0, upperLimit = 3
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
# )
# # Remove y-axis labels
# ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
#   , "\nNormalized expression"
#   , "\nCells sorted by cluster (columns)"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# ggsave(paste0(outGraph, "ExprHeatmap_Norm_Cluster", clusterID, ".png")
#   , width = 12, height = 8)
################################################################################

### DE Heatmaps, top 20 DE table, and save seurat object

# Top 20 DE genes
clusterDe20DF <- data.frame(
  markers_DF %>% group_by(Cluster) %>% top_n(20, Log2_Fold_Change))
# Write out as tab delimited
write.table(x = clusterDe20DF
  , file = paste0(outTable, "ClusterX_Vs_All_Clusters_Top20.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)

# Save seurat object with cluster DE data frame
# save(centSO, noCentExM, metDF, markers_DF, file = paste0(outData, "seuratO.Robj"))

# Violin plots of fold changes by cluster
ggDF <- markers_DF
ggDF$Cluster <- as.factor(ggDF$Cluster)
ggplot(ggDF, aes(x = Cluster, y = Log2_Fold_Change)) +
  geom_violin(aes(fill = Cluster)) +
  geom_jitter(size = 0.05) +
  ylab("Log fold change") +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nFold changes for DE genes by cluster"
    , "\n"))
ggsave(paste0(outGraph, "DE_Violin.png"), width = 12, height = 8)

# Histograms of fold changes by cluster
ggDF <- markers_DF
ggDF$Cluster <- as.factor(ggDF$Cluster)
ggplot(ggDF, aes(x = Log2_Fold_Change)) +
  geom_histogram() +
  facet_wrap(~Cluster) +
  ylab("Counts") +
  xlab("Log fold change") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nHistogram of fold changes for DE genes by cluster"
    , "\n"))
ggsave(paste0(outGraph, "DE_Histogram.png"), width = 12, height = 8)

# MA plots for each cluster
ggL <- lapply(unique(markers_DF$Cluster), function(clusterID) {
  print(clusterID)
  # Subset expression matrix to cluster
  cdf <- noCentExM[ ,centSO@ident == clusterID, drop = FALSE]
  # Subset expression matrix to all other cells
  ndf <- noCentExM[ , ! centSO@ident == clusterID]
  # Fold change
  lfc <- rowMeans(cdf) - rowMeans(ndf)
  # Mean expression
  mn <- rowMeans(noCentExM)
  # Combine in DF for ggplot
  ggDF <- data.frame(LOG_FOLD_CHANGE = lfc, MEAN_EXPRESSION = mn
    , GENE = names(lfc))
  ggDF$DE_GENE <- row.names(ggDF) %in% markers_DF$Gene[
    markers_DF$Cluster == clusterID &
    markers_DF$FDR < 0.05 &
    markers_DF$Log2_Fold_Change > 0.2
    ]
  print(head(ggDF))
  # ggplot
  gg <- ggplot(ggDF, aes(x = MEAN_EXPRESSION, y = LOG_FOLD_CHANGE)) +
    geom_point(size = 0.1, alpha = 0.25, aes(color = DE_GENE)) +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "none") +
    xlab("Mean normalized expression") +
    ylab("Log fold change") +
    ggtitle(clusterID)
  return(gg)
})
# plot_grid combine cluster heatmaps
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMA plots: Cells in cluster Vs all other cells"
  , "\nColor indicates > 0.2 log fold change and FDR < 0.05"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "MAplot.png"), width = 12, height = 4+length(ggL))

# Number of DE genes per cluster barplots
ggDF1 <- data.frame(table(markers_DF$Cluster[
  markers_DF$FDR < 0.05]))
ggDF2 <- data.frame(table(markers_DF$Cluster[
  markers_DF$Log2_Fold_Change > 0.2 &
  markers_DF$FDR < 0.05]))
ggDF3 <- data.frame(table(markers_DF$Cluster[
  markers_DF$Log2_Fold_Change > 0.4 &
  markers_DF$FDR < 0.05]))
ggDF4 <- data.frame(table(markers_DF$Cluster[
  markers_DF$Log2_Fold_Change > 0.7 &
  markers_DF$FDR < 0.05]))
p1 <- ggplot(ggDF1, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("FDR < 0.05")
p2 <- ggplot(ggDF2, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.2 log fold change; FDR < 0.05")
p3 <- ggplot(ggDF3, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.4 log fold change; FDR < 0.05")
p4 <- ggplot(ggDF4, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.7 log fold change; FDR < 0.05")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, p4, ncol = 4, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_Number_Barplot.pdf"), width = 17, height = 6)

# Number of DE genes per cluster versus cluster size
ggDF1 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$FDR < 0.05])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
ggDF2 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.2 &
    markers_DF$FDR < 0.05])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
ggDF3 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.4 &
    markers_DF$FDR < 0.05])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
ggDF4 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.7 &
    markers_DF$FDR < 0.05])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
p1 <- ggplot(ggDF1, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("FDR < 0.05")
p2 <- ggplot(ggDF1, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.2 log fold change; FDR < 0.05")
p3 <- ggplot(ggDF2, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.4 log fold change; FDR < 0.05")
p4 <- ggplot(ggDF3, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.7 log fold change; FDR < 0.05")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, p4, ncol = 4, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster versus number of cells per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_NumberVsCells_ScatterPlot.pdf")
  , width = 17, height = 6)

# Number of DE genes per cluster versus mean genes detected per cluster
ggDF1 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$FDR < 0.05])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggDF2 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.2 &
    markers_DF$FDR < 0.05])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggDF3 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.2 &
    markers_DF$FDR < 0.05])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggDF4 <- data.frame(
  table(markers_DF$Cluster[
    markers_DF$Log2_Fold_Change > 0.2 &
    markers_DF$FDR < 0.05])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
p1 <- ggplot(ggDF1, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("FDR < 0.05")
p2 <- ggplot(ggDF1, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.2 log fold change; FDR < 0.05")
p3 <- ggplot(ggDF2, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.4 log fold change; FDR < 0.05")
p4 <- ggplot(ggDF3, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.7 log fold change; FDR < 0.05")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, p4, ncol = 4, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster versus number of cells per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_NumberVsnGene_ScatterPlot.pdf")
  , width = 17, height = 6)
################################################################################

### Hierarchical cluster by Seurat cluster mean expression of DE genes

# Top DE genes
clusterDeTopDF <- data.frame(
  markers_DF %>% group_by(Cluster) %>% top_n(20, Log2_Fold_Change))

# Subset to clusters of interest
clusterDeTopDF <- clusterDeTopDF[clusterDeTopDF$Cluster != 16, ]

# Subset expression matrix
df <- data.frame(t(centSO@scale.data[
  row.names(centSO@scale.data) %in% markers_DF$Gene, ]))
df$ClusterID <- centSO@ident
# Subset to clusters of interest
df <- df[df$ClusterID %in% clusterDeTopDF$Cluster, ]
# Mean
df <- aggregate(.~ClusterID, df, mean)
row.names(df) <- df$ClusterID
df <- df[ ,colnames(df) != "ClusterID"]

# Obtain the dendrogram
dend <- as.dendrogram(hclust(
  d = dist(df, method = "euclidean")
  , method = "ward.D2"))
# Extract dend data
dend_data <- dendro_data(dend)
# Exract the hcluster order
hclust_order <- as.numeric(as.character(dend_data$labels$label))

# Plot dendrogram
# ggplot
ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
    hjust = 1, angle = 90, size = 3)
ggsave(paste0(outGraph, "hclust_dend.pdf"))

# Heatmap plot
# Use the dendrogram label data to position the labels
# Order genes by dendro order
geneGroupDF <- data.frame(GENE = clusterDeTopDF$Gene, GROUP = clusterDeTopDF$Cluster)
geneGroupDF$GROUP <- factor(geneGroupDF$GROUP, levels = hclust_order)
geneGroupDF <- geneGroupDF[order(geneGroupDF$GROUP), ]
# Plot heatmap
ggL <- lapply(hclust_order, function(cluster){
  # tryCatch(
    Heatmap_By_Cluster(
      geneGroupDF = geneGroupDF
      , exprM = as.matrix(centSO@scale.data)
      , seuratO = centSO
      , clusters = cluster
      , lowerLimit = -1.5
      , upperLimit = 1.5
      , geneOrder = TRUE
    )
    # , error = function(e) NULL)
})
# Remove nulls from ggplot list
ggL <- ggL[! sapply(ggL, is.null)]
# Extract legend
legend <- get_legend(ggL[[1]])
# Duplicate one plot to save labels
ggLabels <- ggL[[1]]
# Remove axis labels
ggL[1:length(ggL)] <- lapply(ggL[1:length(ggL)], function(gg) {
  # gg <- gg + facet_grid(~SEURAT_ClusterS, space = "free", scales = "free")
  gg + theme(
    strip.text.y = element_blank()
    , axis.text.y = element_blank()
    , legend.position = "none"
    , axis.title.y = element_blank()
    , axis.ticks.y = element_blank()
    , axis.title.x = element_blank()
    # margin: top, right, bottom, and left
    , plot.margin = unit(c(1, 0.02, 1, 0.02), "cm")
  )
})
# Add genes as labels by plotting and reducing margin to cutoff heatmap
ggLabels <- lbak
ggLabels <- ggLabels + theme(
  strip.text.y = element_blank()
  # , axis.text.y = element_blank()
  , legend.position = "none"
  , axis.title.y = element_blank()
  , axis.ticks.y = element_blank()
  , axis.title.x = element_blank()
  # margin: top, right, bottom, and left
  , plot.margin = unit(c(1, -50, 1, 1), "cm")
)

# Combine individual heatmaps and dendrogram
# # Determine relative widths by scaling to log 10 cluster size
# df <- data.frame(table(centSO@ident))
# df$rel_widths <- as.vector(log((df$Freq + 1), 10)) + 1
# df <- df[match(hclust_order, df$Var1), ]
# rel_widths <- c(30, df$rel_widths, 1)
rel_widths <- c(30, rep(3, length(hclust_order)), 1)
# Combine
pg <- plot_grid(plotlist = append(list(ggLabels), ggL), ncol = 19
  , rel_widths = rel_widths, align = 'h', axis = 't')
# Add legend
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.1))
# Add title
title <- paste0(graphCodeTitle
  , "\n\nTop 10 DE genes for each cluster"
  , "\nMean centered, variance scaled normalized expression"
  , "\nCells sorted by hierarchical clustering of Seurat clusters (columns)")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))

ggsave(paste0(outGraph, "hclust_heatmap.png"), width = 12, height = 56
  , limitsize = FALSE)

## Markers in heatmap

Intersect_Top_DE_Genes_And_Markers <- function(clusters, groupings){
  genes <- intersect(
    clusterDeTopDF$Gene[clusterDeTopDF$Cluster %in% clusters]
    , kmDF$Gene.Symbol[kmDF$Grouping %in% groupings]
  )
  return(genes)
}
# "RGS5"
Intersect_Top_DE_Genes_And_Markers(
  clusters = 14
  , groupings = c("Pericyte")
)
# "CCL3" "AIF1"
Intersect_Top_DE_Genes_And_Markers(
  clusters = 15
  , groupings = c("Microglia")
)
# "ITM2A" "CLDN5" "ESAM"
Intersect_Top_DE_Genes_And_Markers(
  clusters = 12
  , groupings = c("Endothelial Cell")
)
# "SOX5"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(3,13)
  , groupings = c("Excitatory Deep Layer Cortical", "Neuron")
)
# "DLX1" "DLX2" "LHX6" "DLX5"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(5,6)
  , groupings = c("GABAergic interneuron")
)
# "STMN2"   "NEUROD6" "SATB2"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(4)
  , groupings = c("Neuron", "Excitatory Upper Layer Cortical")
)
# [1] "PPP1R17" "SSTR2"   "EOMES"   "PENK"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(2)
  , groupings = c("IP")
)
# [1] "SATB2"   "STMN2"   "NEUROD6"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(1)
  , groupings = c("Neuron", "Excitatory Upper Layer Cortical")
)
# "NEUROD6" "POU3F2"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(0)
  , groupings = c("Neuron", "Excitatory Upper Layer Cortical")
)
# [1] "HMGB2" "SOX2"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(8,10)
  , groupings = c("RG", "IP")
)
# [1] "PTPRZ1" "OLIG1"  "PDGFRA"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(11)
  , groupings = c("OPC")
)
# "VIM"    "PTPRZ1" "SOX2"   "SLC1A3" "HES1"   "HOPX"
Intersect_Top_DE_Genes_And_Markers(
  clusters = c(7,9)
  , groupings = c("RG", "oRG", "vRG")
)
################################################################################

### Heatmaps

# Top 10 markers for each cluster
markers_DF %>% group_by(Cluster) %>% top_n(10, Log2_Fold_Change) -> top10
geneGroupDF <- data.frame(GENE = as.character(top10$Gene)
  , GROUP = top10$Cluster)

Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, clusters = c(1:15), lowerLimit = -1.5, upperLimit = 1.5
  # , geneOrder = TRUE
  , clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
ggsave(paste0(outGraph, "Top10DE_ExprHeatmap_CentScale.png")
  , width = 13, height = 24)
################################################################################
#
# ### DE genes expression across Kang et al. cortex stages 1-8
#
# ## Plot DE genes expression across Kang et al. cortex stages 1-8
#
# ggL <- lapply(sort(unique(markers_DF$Cluster)), function(cluster){
#   # Subset DE data frame to cluster
#   specificClusterDeDF <- markers_DF[markers_DF$Cluster == cluster, ]
#   # Subset to genes > 0.4 log fold change
#   specificClusterDeDF <- specificClusterDeDF[specificClusterDeDF$Log2_Fold_Change > 0.4, ]
#   # Subset Kang expression matrix to DE genes
#   exM <- kangExM[row.names(kangExM) %in% specificClusterDeDF$Ensembl, ]
#   # Format for ggplot2 and add stage
#   df <- melt(exM)
#   df$Stage <- kangMtDF$Stage[match(df$Var2, kangMtDF$X)]
#   # Plot
#   gg <- ggplot(df, aes(x = Stage, y = value)) +
#     geom_jitter(size = 0.01, alpha = 0.2) +
#     stat_summary(geom = "pointrange", fun.data = mean_cl_normal,
#       fun.args = list(conf.int = 0.95), color = "red", fatten = 0.25) +
#     # stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point") +
#     # geom_smooth(method = "loess") +
#     ylab("Normalized expression") +
#     ggtitle(paste0("Cluster ", cluster))
#   return(gg)
# })
# Plot_Grid(ggPlotsL = ggL, ncol = 4, rel_height = 0.15
#   , title = paste0(graphCodeTitle
#     , "\n\nSeurat cluster DE genes expression in Kang cortex stage 1-8"
#     , "\nDE > 0.4 log fold change"
#     , "\nRed points = mean"
#     , "\nRed bar = 95% confidence intervals"
#     , "\nBlack points = expression of each gene in each sample"))
# ggsave(paste0(outGraph, "KangExpr.png"), width = 13, height = 18)
#
#
# ## Correlation to Kang stages
#
# # Convert hgnc symbols to ensembl and leave ensembl IDs unchanged
# # (Gene IDs are a mix of hgnc symbols and ensembl IDs for those genes that have
# # no hgnc symbol)
#
#
# Convert_Mixed_GeneSym_EnsID_To_EnsID <- function(ids){
#   idx <- match(ids, bmDF$hgnc_symbol)
#   ens <- bmDF$ensembl_gene_id[idx]
#   ids[! is.na(ens)] <- as.character(ens[! is.na(ens)])
#   return(ids)
# }
#
# Cluster_Correlation_To_Kang_Stages <- function (topGenes) {
#   ll <- lapply(sort(unique(markers_DF$Cluster)), function(cluster){
#     # Subset expression matrix to cells in cluster
#     exM <- noCentExM[ ,centSO@ident %in% cluster]
#     # Convert hgnc symbols to ensembl and leave ensembl IDs unchanged
#     # (Gene IDs are a mix of hgnc symbols and ensembl IDs for those genes that have
#     # no hgnc symbol)
#     row.names(exM) <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(exM))
#     # Subset to genes above mean expression threshold
#     rMns <- rowMeans(exM)
#     rMns <- sort(rMns, decreasing = TRUE)
#     ids <- names(rMns)[1:topGenes]
#     ssKangExM <- kangExM[row.names(kangExM) %in% ids, ]
#     df <- melt(ssKangExM)
#     df$Stage <- kangMtDF$Stage[match(df$Var2, kangMtDF$X)]
#     df <- aggregate(value~Stage+Var1, df, mean)
#     df <- dcast(df, Var1~Stage)
#     df$Cluster_Mean_Expression <- rMns[match(df$Var1, names(rMns))]
#     scor <- apply(df[ ,-1], 2, function(col){
#       cor(df$Cluster_Mean_Expression, col, method = "spearman")
#     })
#     return(scor)
#   })
#   df <- do.call("cbind", ll)
#   df <- df[-9, ]
#   colnames(df) <- sort(unique(markers_DF$Cluster))
#   df <- melt(df)
#   return(df)
# }
#
# Plot_Cluster_Correlation_To_Kang_Stages <- function (ggDF, topGenes) {
#   ggplot(ggDF, aes(x = Var1, y = Var2, fill = value)) +
#     geom_tile() +
#     scale_fill_gradient(low = "white", high = "red", space = "Lab"
#       , name = "Spearman") +
#     geom_text(aes(label = round(value, 2))) +
#     scale_x_continuous(breaks = unique(ggDF$Var1)) +
#     scale_y_continuous(breaks = unique(ggDF$Var2)) +
#     xlab("Stage") +
#     ylab("Cluster") +
#     ggtitle(paste0("Used top ", topGenes, " expressed genes in cluster"))
# }
#
# df <- Cluster_Correlation_To_Kang_Stages(topGenes = 2500)
# gg1 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 2500)
#
# df <- Cluster_Correlation_To_Kang_Stages(topGenes = 5000)
# gg2 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 5000)
#
# df <- Cluster_Correlation_To_Kang_Stages(topGenes = 10000)
# gg3 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 10000)
#
# Plot_Grid(ggPlotsL = list(gg1, gg2, gg3), ncol = 3, rel_height = 0.15
#   , title = paste0(graphCodeTitle
#     , "\n\nCorrelation of Seurat cluster mean expression profile to"
#     , "\nKang stage expression")
# )
# ggsave(paste0(outGraph, "KangCorrelation.png"), width = 16, height = 9)
# ################################################################################
#
# ### DE genes across Miller zones
#
# ## Expression across Miller zones
#
# ggL <- lapply(sort(unique(markers_DF$Cluster)), function(cluster){
#   # Subset DE data frame to cluster
#   specificClusterDeDF <- markers_DF[markers_DF$Cluster == cluster, ]
#   # Subset to genes > 0.4 log fold change
#   specificClusterDeDF <- specificClusterDeDF[specificClusterDeDF$Log2_Fold_Change > 0.4, ]
#   # Subset Kang expression matrix to DE genes
#   exM <- millerExDF[row.names(millerExDF) %in% specificClusterDeDF$Ensembl, ]
#   # Format for ggplot2 and add stage
#   df <- melt(exM)
#   df$Zone <- millerMtDF$Zone[match(df$variable, millerMtDF$well_id)]
#   df$Zone <- factor(df$Zone, levels = c("VZ", "SZi", "SZo", "IZ", "SP", "CPi", "CPo", "MZ", "SG"))
#   df <- df[! is.na(df$Zone), ]
#   # Plot
#   gg <- ggplot(df, aes(x = Zone, y = value)) +
#     geom_jitter(size = 0.01, alpha = 0.2) +
#     stat_summary(geom = "pointrange", fun.data = mean_cl_normal,
#       fun.args = list(conf.int = 0.95), color = "lightcoral", alpha = 0.75) +
#     stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point") +
#     # geom_smooth(method = "loess") +
#     ylab("Normalized expression") +
#     ggtitle(paste0("Cluster ", cluster))
#   return(gg)
# })
# Plot_Grid(ggPlotsL = ggL, ncol = 4, rel_height = 0.15
#   , title = paste0(graphCodeTitle
#     , "\n\nSeurat cluster DE genes expression in Miller zones"
#     , "\nDE > 0.4 log fold change"
#     , "\nRed points = mean"
#     , "\nBlue bar = 95% confidence intervals"
#     , "\nBlack points = expression of each gene in each sample"))
# ggsave(paste0(outGraph, "MillerExpr.png"), width = 13, height = 18)
# ################################################################################
