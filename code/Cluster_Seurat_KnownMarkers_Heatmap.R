# Damon Polioudakis
# 2016-12-06
# Clustering of Drop-seq cells by digital gene expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell.csv", header = TRUE)

# Cluster identities
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

## Variables
graphCodeTitle <- "Cluster_Seurat_KnownMarkers_Heatmap.R"
outGraphPfx <- "../analysis/graphs/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_"
################################################################################

### Heatmap of known markers genes

# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
kmDF <- kmDF[! is.na(kmDF$Gene.Symbol), ]
# Filter for gene symbols in data
kmDF <- kmDF[kmDF$Gene.Symbol %in% fetb@data@Dimnames[1][[1]], ]
pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf")
  , width = 8, height = 4 + 0.2*nrow(kmDF))
pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"))
my_palette <- colorRampPalette(c("black", "yellow"))(n = 100)
DoHeatmap(fetb, genes.use = as.character(kmDF$Gene.Symbol)
  , order.by.ident = TRUE
  , slim.col.label = TRUE
  # Turn off mean centering
  , do.scale = FALSE
  # Turn off histogram trace on key
  # , density.info = "none"
  , remove.key = FALSE
  , col.use = my_palette
  , main = paste0(graphCodeTitle
    , "\nKnown Markers"))
dev.off()

# pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"), width = 10, height = 4 + 0.2*nrow(kmDF))
# kmDFL <- split(kmDF, kmDF$Population.Marked)
# sapply(names(kmDFL), function(kmLName) {
#   # print(kmLName)
#   df <- kmDFL[[kmLName]]
#   # print(df)
#   print(as.character(df$Gene.Symbol))
#   DoHeatmap(fetb, genes.use = as.character(df$Gene.Symbol)
#     , order.by.ident = TRUE
#     , slim.col.label = TRUE
#     , remove.key = FALSE
#     , main = (kmLName)
#   )
# })
# dev.off()