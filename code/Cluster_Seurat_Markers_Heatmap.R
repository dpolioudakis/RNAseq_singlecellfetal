# Damon Polioudakis
# 2016-12-08
# Heatmap of all markers from clusters output by Cluster_Seurat.R

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())

require(Seurat)

# Cluster identities
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
################################################################################

# Find markers
fetb.markers <- FindAllMarkers(fetb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)


# Heatmap of markers
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
pdf("..analysis/graphs/Cluster_Seurat_Markers_Heatmap_exon_FtMm250.pdf"
  , width = 10, height = 12)
my_palette <- colorRampPalette(c("black", "yellow"))(n = 100)
DoHeatmap(fetb, genes.use = fetb.markers$gene
  , order.by.ident = TRUE
  , slim.col.label = TRUE
  # Turn off mean centering
  , do.scale = FALSE
  # Turn off histogram trace on key
  , density.info = "none"
  , remove.key = FALSE
  , col.use = my_palette
  , main = "Known Markers")
dev.off()