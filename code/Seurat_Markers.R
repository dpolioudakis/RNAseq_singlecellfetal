# Damon Polioudakis
# 2017-02-13
# Clustering of Drop-seq cells by digital gene expression

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
# require(xlsx)

args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(args)

## Input data

# Seurat object
load(args[1])

## Variables
graphCodeTitle <- "Cluster_Seurat.R"
outGraphPfx <- args[2]
outTablePfx <- args[3]
outRdatPfx <- args[4]
# outGraphPfx <- "../analysis/graphs/Seurat_Markers_exon_FtMm250_"
# outTablePfx <- "../analysis/tables/Seurat_Markers_exon_FtMm250_"
# outRdatPfx <- "../analysis/Seurat_Markers_exon_FtMm250_"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
outTableDir <- dirname(outTablePfx)
outRdatDir <- dirname(outRdatPfx)
dir.create(outGraphDir, recursive = TRUE)
dir.create(outTableDir, recursive = TRUE)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Finding differentially expressed genes (cluster biomarkers)

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

# Find markers for every cluster compared to all remaining cells, report only
# the positive ones
markers <- FindAllMarkers(seuratO, only.pos = TRUE, min.pct = 0.25
  , thresh.use = 0.25, return.thresh = 1)
# markers <- FindAllMarkers(seuratO, only.pos = TRUE, min.pct = 0
#   , min.diff.pct = 0, thresh.use = 0, return.thresh = 1)
markersDF <- data.frame(
  markers %>% group_by(cluster) %>% top_n(20, avg_diff))
# cluster1.markers <- FindMarkers(seuratO, ident.1 = 1, min.pct = 0
# , thresh.use = 0, min.diff.pct = 0)

write.table(x = markers
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = markersDF
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All_Top20.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)
################################################################################