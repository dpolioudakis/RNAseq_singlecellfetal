# Damon Polioudakis
# 2017-09-13
# Subset seurat object for testing

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)

# Input Seurat object
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")

# Subsetting for testing
idx <- sample(1:ncol(centSO@scale.data), 2000)

exDF <- centSO@raw.data
exDF <- exDF[, colnames(exDF) %in% names(centSO@ident)]
exDF <- exDF[ ,idx]
centSO@raw.data <- exDF

centSO@meta.data <- centSO@meta.data[idx, ]

centSO@ident <- centSO@ident[idx]

noCentExM <- noCentExM[ ,idx]

centSO@data <- centSO@data[ ,idx]

centSO@scale.data <- centSO@scale.data[ ,idx]

centSO@dr$tsne@cell.embeddings <- centSO@dr$tsne@cell.embeddings[idx, ]

save(centSO, noCentExM
  , file = "../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO_TEST.Robj")
################################################################################