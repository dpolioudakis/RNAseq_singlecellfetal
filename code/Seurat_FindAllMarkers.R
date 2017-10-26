# Damon Polioudakis
# 2017-08-28
# Run Seurat FindMarkers on all clusters

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
require(irlba)
require(gridExtra)
# require(xlsx)

args <- commandArgs(trailingOnly = TRUE)
print(args)

load("../analysis/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11_exon_FtMm250_200-3sdgd_Mt5_RegNumiPMtLibBrain_PC1to40_seuratO.Robj")

## Variables
clusterID <- (as.numeric(args[1]) - 1)
print(paste0("Cluster ID: ", clusterID))

graphCodeTitle <- "Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11_exon_FtMm250_200-3sdgd_Mt5_RegNumiPMtLibBrain_"
outTable <- "../analysis/tables/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11_exon_FtMm250_200-3sdgd_Mt5_RegNumiPMtLibBrain_"
outData <- "../analysis/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11/Seurat_Cluster_DS-2-3-4-5-6-7-8-9-11_exon_FtMm250_200-3sdgd_Mt5_RegNumiPMtLibBrain_PC1to40_test_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outData)
dir.create(outRdatDir, recursive = TRUE)

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

DE_Heatmap <- function(clusterDeDF, exDF, clusterIDs, ggtitle) {
  # Subset to genes, merge to keep duplicate genes (from more than one cluster)
  ggDF <- merge(clusterDeDF, exDF, by.x = "gene", by.y = "row.names"
    , all.x = TRUE)
  # For gene ordering by cluster
  idx <- match(rev(clusterDeDF$gene), ggDF$gene)
  ggDF <- ggDF[idx, ]
  row.names(ggDF) <- paste0(length(ggDF$gene):1, "_", ggDF$gene)
  ggDF <- ggDF[ ,-c(1:6)]
  # Order by clustering
  idx <- match(colnames(ggDF), names(sort(clusterIDs)))
  ggDF <- ggDF[ ,idx]
  # Format for ggplot2
  ggDF <- as.matrix(ggDF)
  ggDF <- melt(ggDF)
  # Add clusters
  idx <- match(ggDF$Var2, names(clusterIDs))
  ggDF$CLUSTERS <- clusterIDs[idx]
  # # Duplicate CLUSTERS variable to facet genes by cluster as well
  # ggDF$CLUSTERS2 <- ggDF$CLUSTERS
  # Set sample order by clustering
  ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(clusterIDs)))
  # Set expression limits
  ggDF$value[ggDF$value > 3] <- 3
  ggDF$value[ggDF$value < 0] <- 0
  # # Split by clusters with > or < 1000 cells
  # cl <- names(table(clusterIDs))[table(clusterIDs) < 1000]
  # id <- names(clusterIDs)[clusterIDs %in% cl]
  # gg2DF <- ggDF[ggDF$Var2 %in% id, ]
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(~CLUSTERS, scales = "free") +
    # facet_grid(~CLUSTERS, space = "free", scales = "free") +
    # scale_fill_gradient2(name = "Normalized\nexpression"
    #   , high = "#d7191c", low = "white") +
    scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(ggtitle)
}
################################################################################

### Finding differentially expressed genes (cluster biomarkers)

print("### Finding differentially expressed genes for cluster ", clusterID)

# Seurat can help you find markers that define clusters via differential
# expression. By default, it identifes positive and negative markers of a single
# cluster (specified in ident.1), compared to all other cells. 
# **FindAllMarkers()** automates this process for all clusters, but you can also
# test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a gene to be detected at a minimum percentage in
# either of the two groups of cells, and the thresh.test argument requires a
# gene to be differentially expressed (on average) by some amount between the
# two groups. You can set both of these to 0, but with a dramatic increase in
# time - since this will test a large number of genes that are unlikely to be
# highly discriminatory. As another option to speed up these computations,
# max.cells.per.ident can be set. This will downsample each identity class to
# have no more cells than whatever this is set to. While there is generally
# going to be a loss in power, the speed increases can be signficiant and the
# most highly differentially expressed genes will likely still rise to the top.

# Find DE genes of cluster
clusterDeDF <- FindMarkers(centSO, ident.1 = clusterID, only.pos = TRUE
  , test.use = "negbinom", latent.vars = c("nUMI", "librarylab", "individual")
  , min.pct = 0.25, thresh.use = 0.25)
clusterDeDF$gene <- row.names(clusterDeDF)
clusterDeDF$cluster <- clusterID

# Write to tab delimited table
# If 1st cluster write new file, otherwise append to file
if (clusterID == min(as.numeric(unique(centSO@ident)))) {
  write.table(x = clusterDeDF
    , file = paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All.txt")
    , sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(x = clusterDeDF
    , file = paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All.txt")
    , sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
}
################################################################################

### DE Heatmaps, top 20 DE table, and save seurat object

print("### DE Heatmaps, top 20 DE table, and save seurat object")
  
# If cluster DE is done for all clusters compile DE tables and make DE heatmaps
if (clusterID == max(as.numeric(unique(centSO@ident)))) {
  
  clusterDeDF <- read.table(
    paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All.txt"), header = TRUE)
  
  # Top 20 DE genes
  clusterDe20DF <- data.frame(
    clusterDeDF %>% group_by(cluster) %>% top_n(20, avg_diff))
  # Write out as tab delimited
  write.table(x = clusterDe20DF
    , file = paste0(outTable, "PC1-40_Marker_Genes_Clusters_Vs_All_Top20.txt")
    , sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save seurat object with cluster DE data frame
  save(centSO, noCentExM, metDF, clusterDeDF, file = paste0(outData, "seuratO.Robj"))
  
  # Heatmaps
  # Top 10 markers for each cluster
  clusterDeDF %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
  # Normalized expression - centered and scaled
  p1 <- DE_Heatmap(clusterDeDF = top10
    , exDF = centSO@scale.data
    , clusterIDs = centSO@ident
    , ggtitle = paste0(graphCodeTitle
      , "\n"
      , "\nTop 10 marker genes per cluster"
      , "\nMean centered, variance scaled, normalized expression"
      , "\n")
    )
  ggsave(p1, paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers_CenterScale.png")
    , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)
  # Normalized expression
  p2 <- DE_Heatmap(clusterDeDF = top10
    , exDF = noCentExM
    , clusterIDs = centSO@ident
    , ggtitle = paste0(graphCodeTitle
      , "\n"
      , "\nTop 10 marker genes per cluster"
      , "\nNormalized expression"
      , "\n")
  )
  ggsave(p2, paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers.png")
    , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)
  
  # Plot centered scaled and not centered scaled heatmaps next to each other
  png(paste0(outGraph, "PC1-40_DoHeatmap_Top10Markers_CenterScale.png")
    , width = 12, height = 4 + 0.15*nrow(top10), units = "in", res = 300)
  grid.arrange(p1, p2, ncol = 2, top = paste0(graphCodeTitle
    , "\n"
    , "\nTop 10 marker genes per cluster"
    , "\n"))
  dev.off()
}

# # Violin plots of top 10 markers
# pdf(paste0(outGraph, "ViolinPlot_Top10Markers.pdf"), width = 10)
# top10L <- split(top10, top10$cluster)
# lapply(top10L, function(top10cluster) {
#   VlnPlot(centSO, top10cluster$gene, size.use = 0.5)
# })
# dev.off()
# 
# # Feature plot of top 10 markers
# pdf(paste0(outGraph, "FeaturePlot_Top10Markers.pdf"), width = 10)
# top10L <- split(top10, top10$cluster)
# lapply(top10L, function(top10cluster) {
#   FeaturePlot(centSO, top10cluster$gene, cols.use = c("grey","blue")
#     , pt.size = 0.7)
# })
# dev.off()
################################################################################