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
graphCodeTitle <- "Seurat_PCA_Test_Function_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Seurat_PCA_Test_Function_DS002_003_exon_FtMm250_"
outTablePfx <- "../analysis/tables/Seurat_PCA_Test_Function_DS002_003_exon_FtMm250_"
outRdat <- "../analysis/Seurat_PCA_Test_Function_DS002_003_exon_FtMm250.Robj"

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

### PCA of expression

# Run prcomp PCA
# Center
pcaPcmpCtrDF <- prcomp(t(seuratO@scale.data), center = TRUE, scale = FALSE)
# No centering
pcaPcmpNoCtrDF <- prcomp(t(seuratO@scale.data), center = FALSE, scale = FALSE)

# Seurat PCA
sPcaO <- PCA(seuratO, pc.genes = unlist(seuratO@data@Dimnames[1])
  , do.print = TRUE, pcs.print = 5, genes.print = 5)

# Run PCA with the IRLBA package (iteratively computes the top dimensions
# , dramatic increase in speed since we are throwing away most PCs anyway)
sFastPcaO <- PCAFast(seuratO, pc.genes = unlist(seuratO@data@Dimnames[1])
  , pcs.compute = 40, pcs.print = 10)

# Save PCA as robj to compare prcomp to PCAFast functions
save(pcaPcmpCtrDF, pcaPcmpNoCtrDF, sPcaO, sFastPcaO, file = outRdat)
################################################################################