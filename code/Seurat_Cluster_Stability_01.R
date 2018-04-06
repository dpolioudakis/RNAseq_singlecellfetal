# Damon Polioudakis
# 2017-03-29
# Clustering stability using bootstrapping + Jaccard method from Christian
# Hennig “Cluster-wise assessment of cluster stability”

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(Seurat)
require(Matrix)
require(fpc)
require(methods)
require(dplyr)

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
i <- as.numeric(args[1])
print(paste0("Seed: ", i))
# i <- 3

## Input
# DS-002-003-004-005-006-007-008-009-011: exDF and metDF
load("../analysis/analyzed_data/Expression_Matrix_Compile/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7-8-9-11.Rdata")
# Subset for testing
# idx <- sample(1:30000, 1000)
# exDF <- exDF[ ,idx]
# metDF <- metDF[idx, ]

## Variables
graphCodeTitle <- "Known_Marker_Expression.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"
outTable <- "../analysis/tables/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"
outTable_iterations <- "../analysis/tables/Seurat_Cluster_Stability/DS2-11/Iterations/Seurat_Cluster_Stability_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outTable_iterations), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Sample with replacement then reprocess and cluster

Sample_CellIDs_With_Replacement <- function(){
  print("Sample_CellIDs_With_Replacement")
  # Set seed
  set.seed(i)
  # Sample with replacement
  sampled_cellIDs <- colnames(exDF)[sample(ncol(exDF), replace = TRUE)]
  return(sampled_cellIDs)
}

Sample_ExprMatrix_With_Replacement <- function(sampled_cellIDs){
  print("Sample_ExprMatrix_With_Replacement")
  idx <- match(sampled_cellIDs, colnames(exDF))
  exDF <- exDF[ ,idx]
  return(exDF)
}

Sample_Metadata_With_Replacement <- function(sampled_cellIDs){
  print("Sample_Metadata_With_Replacement")
  idx <- match(sampled_cellIDs, metDF$CELL)
  metDF <- metDF[idx, ]
  return(metDF)
}

Setup_Seurat_Object <- function(){
  print("Setup_Seurat_Object")
  # Initialize the Seurat object
  centSO <- CreateSeuratObject(raw.data = exDF, min.cells = 3, min.genes = 200
    , normalization.method = "LogNormalize", scale.factor = 10000
    , project = "DS-2-3-4-5-6-7-8-9-11", do.scale = FALSE, do.center = FALSE)
  # Add metadata
  metDF$CELL <- colnames(exDF)
  row.names(metDF) <- metDF$CELL
  centSO <- AddMetaData(centSO, metadata = metDF)
  mito.genes <- grep("^MT-", rownames(centSO@data), value = TRUE)
  percent.mito <- apply(
    expm1(centSO@data[mito.genes, ]), 2, sum)/apply(expm1(centSO@data), 2, sum)
  centSO <- AddMetaData(centSO, percent.mito, "percent.mito")
  return(centSO)
}

Filter_Seurat_Object <- function(){
  print("Filter_Seurat_Object")
  high.thresholds <- round(
    mean(centSO@meta.data$nGene) + 3*sd(centSO@meta.data$nGene), 0)
  centSO <- FilterCells(centSO, subset.name = "nGene"
    , high.thresholds = high.thresholds)
  centSO <- FilterCells(centSO, subset.name = "percent.mito"
    , high.thresholds = 0.05)
  return(centSO)
}

Add_Covariates_To_Seurat_Object <- function(){
  print("Add_Covariates_To_Seurat_Object")

  # Individual for regressing out
  individual <- rep("Brain 2", length(centSO@data@Dimnames[[2]]))
  # Brain 3
  br3Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "3"])
  individual[centSO@data@Dimnames[[2]] %in% br3Cells] <- "Brain 3"
  # Brain 4
  br4Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "4"])
  individual[centSO@data@Dimnames[[2]] %in% br4Cells] <- "Brain 4"
  # Brain 5
  br5Cells <- c(as.character(metDF$CELL)[metDF$BRAIN == "5"])
  individual[centSO@data@Dimnames[[2]] %in% br5Cells] <- "Brain 5"
  # Add to Seurat object
  names(individual) <- centSO@data@Dimnames[[2]]
  centSO <- AddMetaData(centSO, individual, "individual")

  # Lab library prep for regressing out
  librarylab <- rep("Plath", length(centSO@data@Dimnames[[2]]))
  dhgCells <- c(as.character(metDF$CELL)[metDF$LIBRARY == "Geschwind"])
  librarylab[centSO@data@Dimnames[[2]] %in% dhgCells] <- "Geschwind"
  names(librarylab) <- centSO@data@Dimnames[[2]]
  centSO <- AddMetaData(centSO, librarylab, "librarylab")

  return(centSO)
}

## Run
sampled_cellIDs <- Sample_CellIDs_With_Replacement()
exDF <- Sample_ExprMatrix_With_Replacement(sampled_cellIDs = sampled_cellIDs)
metDF <- Sample_Metadata_With_Replacement(sampled_cellIDs = sampled_cellIDs)
centSO <- Setup_Seurat_Object()
centSO <- Filter_Seurat_Object()
centSO <- Add_Covariates_To_Seurat_Object()
centSO <- ScaleData(centSO
  , vars.to.regress = c("nUMI", "librarylab", "individual")
)
centSO <- FindVariableGenes(centSO, mean.function = ExpMean
  , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3
  , y.cutoff = 0.5)
centSO <- RunPCA(object = centSO, pc.genes = centSO@var.genes, pcs.compute = 50
  , pcs.print = 1:5, genes.print = 5, maxit = 500, weight.by.var = FALSE)
centSO <- FindClusters(centSO, dims.use = 1:40, resolution = 0.54
  , print.output = 0, save.SNN = TRUE)
cell_cluster_IDs_DF <- data.frame(
  Cell_ID = names(centSO@ident)
  , Cluster = centSO@ident
)
# Write out
write.csv(cell_cluster_IDs_DF
  , file = paste0(outTable_iterations, "ClusterIDs_", i, ".csv")
  , row.names = FALSE
)
################################################################################
