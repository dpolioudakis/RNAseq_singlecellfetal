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

## Input
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb

## Variables
pcs <- c(1:10)
res <- 3.0

graphCodeTitle <- "Seurat_tSNE_Cluster_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Seurat_tSNE_Cluster_DS002_003_exon_FtMm250_200gd_2500dg_Mt5_Pc110_Res3_"
outTablePfx <- "../analysis/tables/Seurat_tSNE_Cluster_DS002_003_exon_FtMm250_200gd_2500dg_Mt5_Pc110_Res3_"
outRdat <- "../analysis/Seurat_tSNE_Cluster_DS002_003_exon_FtMm250_200gd_2500dg_Mt5_Pc110_Res3.Robj"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
dir.create(outGraphDir, recursive = TRUE)
outTableDir <- dirname(outTablePfx)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdat)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Cluster, run tSNE, plot, and save Seurat object

# Cluster
seuratO <- FindClusters(seuratO, pc.use = pcs, resolution = res, print.output = 0
  , save.SNN = T)
# Run tSNE
seuratO <- RunTSNE(seuratO, dims.use = pcs, do.fast = T)
# tSNE plot
pdf(paste0(outGraphPfx, "tSNE.pdf"))
ggO <- TSNEPlot(seuratO, pt.size = 0.75, do.return = TRUE)
ggO + ggtitle(paste0("PCs: 1-", pcs[length(pcs)], " Resolution: ", res))
dev.off()
# Save seurat object
save(seuratO, file = outRdat)
################################################################################