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

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)

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

# Filter expression matrix by:
# percent of cells in cluster 1 gene is expressed in
# fold change of gene in cluster 1 versus cluster 2
DE_Filters_ExpMatrix <- function(
  minPercent, foldChange, clusterID1, clusterID2) {

  # Expressed > 0 counts in > X% of cells in cluster
  # Subset expression matrix to cluster
  cdf <- as.matrix(centSO@data)[ ,centSO@ident %in% clusterID1]
  # Expressed > 0 counts in > X% of cells in cluster
  idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (minPercent / 100)
  print(paste0("Genes expressed in > ", minPercent, "% of cells in cluster"))
  print(table(idxp))

  # Fold change > Y of gene in cluster versus all other cells
  # Subset expression matrix to cluster1
  c1df <- noCentExM[ ,centSO@ident %in% clusterID1]
  # Subset expression matrix to cluster2
  c2df <- noCentExM[ ,centSO@ident %in% clusterID2]
  # Fold change
  v1 <- rowMeans(c1df) - rowMeans(c2df)
  idxf <- v1 > foldChange
  print(paste0("Genes > ", foldChange, " fold change in cluster 1 versus cluster 2"))
  print(table(idxf))

  # Filter exDF
  exDF <- as.matrix(centSO@data[idxp & idxf, ])
  return(exDF)
}

# Format output of linear model into data frame
Format_DE <- function (deLM, classID1, classID2, clusterID1, clusterID2) {
  # Combine log2 fold changes, p-values
  deDF <- data.frame(Gene = row.names(deLM$coefmat)
    , Log2_Fold_Change = deLM$coefmat[ ,2]
    , Pvalue = deLM$pvalmat[ ,2])
  # Order by pvalue
  deDF <- deDF[order(deDF$Pvalue), ]
  # Add cluster ID
  deDF$Class1 <- classID1
  deDF$Class2 <- classID2
  # Percent of cells in cluster expressing gene > 0 counts
  cdf <- as.matrix(centSO@data)[
    row.names(centSO@data) %in% deDF$Gene, centSO@ident %in% clusterID1]
  idx <- match(deDF$Gene, rownames(cdf))
  cdf <- cdf[idx, ]
  deDF$Percent_Cluster <- (rowSums(cdf > 0) / ncol(cdf)) * 100
  deDF$Percent_Cluster <- round(deDF$Percent_Cluster, 1)
  # Percent of all cells expressing gene > 0 counts
  percent_all <- (
    rowSums(as.matrix(centSO@data)[row.names(centSO@data) %in% deDF$Gene, ] > 0)
    / ncol(centSO@data)) * 100
  names(percent_all) <- row.names(centSO@data)[
    row.names(centSO@data) %in% deDF$Gene]
  idx <- match(deDF$Gene, names(percent_all))
  percent_all <- percent_all[idx]
  deDF$Percent_All <- percent_all
  deDF$Percent_All <- round(deDF$Percent_All, 1)
  return(deDF)
}

DE_Heatmap <- function(clusterDeDF, exDF, clusterIDs, ggtitle, upLim, lowLim) {
  # Subset to genes, merge to keep duplicate genes (from more than one cluster)
  ggDF <- merge(clusterDeDF, exDF, by.x = "Gene", by.y = "row.names"
    , all.x = TRUE)
  # For gene ordering by cluster
  idx <- match(rev(clusterDeDF$Gene), ggDF$Gene)
  ggDF <- ggDF[idx, ]
  row.names(ggDF) <- paste0(length(ggDF$Gene):1, "_", ggDF$Gene)
  ggDF <- ggDF[ ,-c(1:ncol(clusterDeDF))]
  # Order by clustering
  idx <- match(colnames(ggDF), names(sort(clusterIDs)))
  ggDF <- ggDF[ ,idx]
  # Format for ggplot2
  ggDF <- as.matrix(ggDF)
  ggDF <- melt(ggDF)
  # Add clusters
  idx <- match(ggDF$Var2, names(clusterIDs))
  ggDF$ClusterS <- clusterIDs[idx]
  # Set sample order by clustering
  ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(clusterIDs)))
  # Set expression limits
  ggDF$value[ggDF$value > upLim] <- upLim
  ggDF$value[ggDF$value < lowLim] <- lowLim
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(~ClusterS, scales = "free") +
    # facet_grid(~ClusterS, space = "free", scales = "free") +
    # scale_fill_gradient2(name = "Normalized\nexpression"
    #   , high = "#d7191c", low = "white") +
    scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    # theme(axis.text.y = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(ggtitle)
}
################################################################################

### Finding differentially expressed genes (cluster biomarkers)

idx <- match(names(centSO@ident), colnames(noCentExM))
noCentExM <- noCentExM[ ,idx]

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
class_cellID <- centSO@ident
class_cellID <- names(class_cluster_idx)[match(class_cellID, class_cluster_idx)]
names(class_cellID) <- names(centSO@ident)

## Find DE genes of cluster versus each other cluster
# List of IDs to loop through minus cluster ID to calculate DE versus
class_ids <- as.character(unique(class_cellID))
class_ids <- class_ids[! is.na(class_ids)]
classID1 <- class_ids[as.numeric(args[1])]
class_ids <- class_ids[! class_ids %in% classID1]
print(paste0("### Finding differentially expressed genes for: ", classID1))

# Loop through list of IDs and calculate DE
# Return list of data frames, each data frame is DE clusterA vs clusterB
ldf <- lapply(class_ids, function(classID2) {
  clusterID1 <- class_cluster_idx[names(class_cluster_idx) == classID1]
  clusterID2 <- class_cluster_idx[names(class_cluster_idx) == classID2]
  print(clusterID1)
  print(clusterID2)

  # Filter cells
  print("Filtering cells")
  exDF <- DE_Filters_ExpMatrix(minPercent = 10, foldChange = 0
    , clusterID1 = clusterID1, clusterID2 = clusterID2)

  # DE Linear model
  print("DE via linear model")
  termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual")]
  # Add term cluster 1 / cluster 2
  termsDF$cluster <- NA
  termsDF$cluster[centSO@ident %in% clusterID1] <- "cluster 1"
  termsDF$cluster[centSO@ident %in% clusterID2] <- "cluster 2"
  mod <- "y ~ cluster+nUMI+librarylab+individual"
  # Subset to cluster 1 and cluster 2 cells
  v1 <- names(centSO@ident[centSO@ident %in% c(clusterID1, clusterID2)])
  termsDF <- termsDF[row.names(termsDF) %in% v1, ]
  exDF <- exDF[ ,colnames(exDF) %in% v1]
  # Set factor levels of cluster 1/2
  termsDF$cluster <- factor(termsDF$cluster, levels = c("cluster 2", "cluster 1"))
  print(head(termsDF))
  print(str(exDF))
  deLM <- DE_Linear_Model(exDatDF = exDF, termsDF = termsDF, mod = mod)

  # Format LM output into data frame
  deDF <- Format_DE(deLM, classID1, classID2, clusterID1, clusterID2)
  head(deDF)

  # Add ensembl
  deDF$Ensembl <- Convert_Mixed_GeneSym_EnsID_To_EnsID(as.character(deDF$Gene))

  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method = "BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)

  return(deDF)
})
de_PW_DF <- do.call("rbind", ldf)

# Write to tab delimited table
print("Writing DE pairwise table to text file")
write.csv(x = de_PW_DF
  , file = paste0(outTable, "Class_", gsub(" ", "", classID1), ".txt")
  , quote = FALSE, row.names = FALSE)

# Subset to genes with pvalue < 0.05
de_PW_DF <- de_PW_DF[de_PW_DF$FDR < 0.05, ]
de_PW_DF <- de_PW_DF[de_PW_DF$Log2_Fold_Change > 0, ]

# Subset to genes DE for the cluster vs each other cluster
de_PW_DFL <- split(de_PW_DF, de_PW_DF$Gene)
idx <- as.vector(sapply(de_PW_DFL, nrow) == length(unique(de_PW_DF$Class2)))
de_PW_DFL <- de_PW_DFL[idx]
de_PW_DF <- do.call("rbind", de_PW_DFL)

# Writing to tab delimited table
print("Writing marker table to text file")
write.csv(x = de_PW_DF
  , file = paste0(
    outTable, "Markers_Class_", gsub(" ", "", classID1), ".txt")
  , quote = FALSE, row.names = FALSE)
################################################################################
#
# ### Plots
#
# # Expression heatmap of marker genes
# # Centered scaled
# p1 <- DE_Heatmap(clusterDeDF = mkDF
#   , exDF = centSO@scale.data
#   , clusterIDs = centSO@ident
#   , upLim = 1.5
#   , lowLim = -1.5
#   , ggtitle = paste0(
#     "\nMean centered, variance scaled, normalized expression"
#     , "\n")
# )
# # Not centered scaled
# p2 <- DE_Heatmap(clusterDeDF = mkDF
#   , exDF = noCentExM
#   , clusterIDs = centSO@ident
#   , upLim = 3
#   , lowLim = 0
#   , ggtitle = paste0(
#     "\nNormalized expression"
#     , "\n")
# )
# # plot_grid
# pg <- plot_grid(p1, p2, ncol = 2)
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nMarker genes for cluster ", clusterID1, " (pairwise DE)"
#   , "\nMarker = gene DE < 0.05 p-value for cluster versus each other cluster"
#   , "\nCells sorted by cluster (columns)"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# ggsave(paste0(outGraph, "ExprHeatmap_Cluster", clusterID1, ".png")
#   , width = 12, height = 8)
# ################################################################################
#
# ### DE Heatmaps, top 20 DE table, and save seurat object
#
# print("### DE Heatmaps, top 20 DE table, and save seurat object")
#
# # If cluster DE is done for all clusters compile DE tables and make DE heatmaps
# if (clusterID1 == max(as.numeric(as.character(unique(centSO@ident))))) {
#
#   # Cluster pairwise DE
#   ldf <- lapply(sort(unique(centSO@ident)), function(clid) {
#     read.table(
#       paste0(outTable, "PC1-40_Cluster", clid, "_Vs_All_Clusters.txt")
#       , header = TRUE)
#   })
#   clusterDeDF <- do.call("rbind", ldf)
#
#   # Write out as tab delimited
#   write.table(x = clusterDeDF
#     , file = paste0(outTable, "PC1-40_ClusterAll_Vs_All_Clusters.txt")
#     , sep = "\t", quote = FALSE, row.names = FALSE)
#
#
#   # # Cluster markers
#   # l1 <- split(clusterDeDF, clusterDeDF$Cluster1)
#   # ldf <- lapply(names(l1), function(clusterID1) {
#   #   df1 <- l1[[clusterID1]]
#   #   # Marker genes
#   #   # Intersection of DE genes across all pairwise comparisons for cluster
#   #   l <- split(df1, df1$Cluster2)
#   #   l <- lapply(l, function(df) {df$Gene})
#   #   if (! identical(Reduce(intersect, l), character(0))) {
#   #     mkDF <- data.frame(Gene = Reduce(intersect, l), Cluster = clusterID1)
#   #     return(mkDF)
#   #   }
#   # })
#   # mkDF <- do.call("rbind", ldf)
#
#   # Cluster markers
#   ldf <- lapply(sort(unique(centSO@ident)), function(clid) {
#     path <- paste0(outTable, "PC1-40_Markers_Cluster", clid, ".txt")
#     if (file.exists(path)) {
#       read.table(path, header = TRUE)
#     }
#   })
#   mkDF <- do.call("rbind", ldf)
#
#   # # Top 20 DE genes
#   # clusterDe20DF <- data.frame(
#   #   clusterDeDF %>% group_by(Cluster) %>% top_n(20, Log2_Fold_Change))
#   # # Write out as tab delimited
#   # write.table(x = clusterDe20DF
#   #   , file = paste0(outTable, "PC1-40_ClusterAll_Vs_All_Clusters_Top20.txt")
#   #   , sep = "\t", quote = FALSE, row.names = FALSE)
#
#   # Heatmaps - format data
#   ggDF <- merge(mkDF, exDF, by.x = "Gene", by.y = "row.names"
#     , all.x = TRUE)
#   # For gene ordering by cluster
#   idx <- match(rev(mkDF$Gene), ggDF$Gene)
#   ggDF <- ggDF[idx, ]
#   row.names(ggDF) <- paste0(length(ggDF$Gene):1, "_", ggDF$Gene)
#   ggDF <- ggDF[ ,-c(1:ncol(mkDF))]
#   # Order by clustering
#   idx <- match(colnames(ggDF), names(sort(clusterIDs)))
#   ggDF <- ggDF[ ,idx]
#   # Format for ggplot2
#   ggDF <- as.matrix(ggDF)
#   ggDF <- melt(ggDF)
#   # Add clusters
#   idx <- match(ggDF$Var2, names(clusterIDs))
#   ggDF$ClusterS <- clusterIDs[idx]
#   # Set sample order by clustering
#   ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(clusterIDs)))
#
#   # Heatmaps - centered scaled
#   # # Top 10 markers for each cluster
#   # clusterDeDF %>% group_by(Cluster) %>% top_n(10, Log2_Fold_Change) -> top10
#   # Split by cluster
#   ldf <- split(mkDF, mkDF$Cluster)
#   # Heatmap for each cluster
#   ggL <- lapply(names(ldf), function(cl) {
#     deDF <- ldf[[cl]]
#     gg <- DE_Heatmap(clusterDeDF = deDF
#       , exDF = centSO@scale.data
#       , clusterIDs = centSO@ident
#       , upLim = 1.5
#       , lowLim = -1.5
#       , ggtitle = paste0("Cluster: ", cl)
#     )
#     gg <- gg + theme(axis.text.y = element_text(size = 10)) +
#       return(gg)
#   })
#   # extract the legend from one of the plots
#   legend <- get_legend(ggL[[1]])
#   # Remove legends from plots
#   ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = ggL, ncol = 2)
#   # add the legend to the row we made earlier. Give it one-third of the width
#   # of one plot (via rel_widths).
#   pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nMarker genes for each cluster determined by pairwise DE (log FC > 0.4)"
#     , "\nMean centered, variance scaled normalized expression"
#     , "\nCells sorted by cluster (columns)"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
#   ggsave(paste0(outGraph, "ExprHeatmap_CentScale.png")
#     , width = 13, height = length(ggL)*2.5)
#
#   # Heatmaps - no center scale
#   # Top 10 markers for each cluster
#   clusterDeDF %>% group_by(Cluster) %>% top_n(10, Log2_Fold_Change) -> top10
#   # Split by cluster
#   ldf <- split(top10, top10$Cluster)
#   # Heatmap for each cluster
#   ggL <- lapply(names(ldf), function(cl) {
#     deDF <- ldf[[cl]]
#     gg <- DE_Heatmap(clusterDeDF = deDF
#       , exDF = noCentExM
#       , clusterIDs = centSO@ident
#       , upLim = 3
#       , lowLim = 0
#       , ggtitle = paste0("Cluster: ", cl)
#     )
#     gg <- gg + theme(axis.text.y = element_text(size = 10)) +
#       return(gg)
#   })
#   # extract the legend from one of the plots
#   legend <- get_legend(ggL[[1]])
#   # Remove legends from plots
#   ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = ggL, ncol = 2)
#   # add the legend to the row we made earlier. Give it one-third of the width
#   # of one plot (via rel_widths).
#   pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nTop 10 DE genes for each cluster"
#     , "\nNormalized expression"
#     , "\nCells sorted by cluster (columns)"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
#   ggsave(paste0(outGraph, "ExprHeatmap_NoCentScale.png")
#     , width = 13, height = length(ggL)*2.5)
#
# }
#
# # # Violin plots of top 10 markers
# # pdf(paste0(outGraph, "ViolinPlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   VlnPlot(centSO, top10cluster$gene, size.use = 0.5)
# # })
# # dev.off()
# #
# # # Feature plot of top 10 markers
# # pdf(paste0(outGraph, "FeaturePlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   FeaturePlot(centSO, top10cluster$gene, cols.use = c("grey","blue")
# #     , pt.size = 0.7)
# # })
# # dev.off()
# ################################################################################
