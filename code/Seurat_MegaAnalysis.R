# Damon Polioudakis
# 2017-09-29
# Clustering of Drop-seq cells by digital gene expression

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
require(cowplot)
# require(xlsx)
source("Function_Library.R")

# load data
seqwell.data <- read.table(file = paste0("~/Downloads/IntegratedAnalysis_ExpressionMatrices/", 
  "pbmc_SeqWell.expressionMatrix.txt"))
tenx.data <- read.table(file = paste0("~/Downloads/IntegratedAnalysis_ExpressionMatrices/", 
  "pbmc_10X.expressionMatrix.txt"))

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

seqwell <- CreateSeuratObject(raw.data = seqwell.data)
seqwell <- NormalizeData(object = seqwell)
seqwell <- ScaleData(object = seqwell)
seqwell <- FindVariableGenes(object = seqwell, do.plot = FALSE)

tenx <- CreateSeuratObject(raw.data = tenx.data)
tenx <- NormalizeData(object = tenx)
tenx <- ScaleData(object = tenx)
tenx <- FindVariableGenes(object = tenx, do.plot = FALSE)

# we will take the union of the top 2k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
hvg.seqwell <- rownames(x = head(x = seqwell@hvg.info, n = 2000))
hvg.tenx <- rownames(x = head(x = tenx@hvg.info, n = 2000))
hvg.union <- union(x = hvg.seqwell, y = hvg.tenx)

# lastly, we set the 'protocol' in each dataset for easy identification
# later it will be transferred to the merged object in RunCCA
tenx@meta.data[, "protocol"] <- "10X"
seqwell@meta.data[, "protocol"] <- "SeqWell"