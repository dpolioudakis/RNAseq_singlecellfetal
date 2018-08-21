# Damon Polioudakis
# 2018-06-14
# Plot known marker lists as heatmaps, violin plots, feature plots

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
require(reshape2)
require(gridExtra)
require(ggplot2)
require(cowplot)
source("Function_Library.R")

## Inputs

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

## Variables
graphCodeTitle <- "Seurat_Cluster_Expressed.R"
outTable <- "../analysis/tables/Seurat_Cluster_Expressed/Seurat_Cluster_Expressed_"

## Output Directories
dir.create(dirname(outTable), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Table of genes expressed in cluster >10%

Rank_Genes_By_Expression_For_Cluster <- function(
  cluster, cluster_annot, seuratO){
    print(paste0("Rank_Genes_By_Expression_For_Cluster cluster: ", cluster))

    # Subset to cells in cluster
    ss_cell_IDs <- names(seuratO@ident)[seuratO@ident %in% cluster]
    ex_M <- as.matrix(seuratO@data)[ ,colnames(seuratO@data) %in% ss_cell_IDs]

    # Remove MALAT1
    ex_M <- ex_M[! rownames(ex_M) %in% "MALAT1", ]

    # Start table with mean expression for each gene
    express_rank_DF <- data.frame(rev(sort(rowMeans(ex_M))))
    # Add cluster annotation
    express_rank_DF$Cluster <- cluster_annot[names(cluster_annot) %in% cluster]
    # Add gene and rank
    express_rank_DF$Gene <- rownames(express_rank_DF)
    express_rank_DF$Rank <- c(1:nrow(express_rank_DF))

    # Percent of cells gene detected in for each gene
    expressed_M <- as.matrix(seuratO@data)[
      ,colnames(seuratO@data) %in% ss_cell_IDs] > 0
    express_percent <- (rowSums(expressed_M) / ncol(expressed_M)) * 100
    express_percent <- round(express_percent, 1)
    idx <- match(express_rank_DF$Gene, names(express_percent))
    express_rank_DF$Percent_of_Cells <- express_percent[idx]

    # Subset to genes expressed in >= 10% of cells
    express_rank_DF <- express_rank_DF[express_rank_DF$Percent_of_Cells >= 10, ]

    # Cleanup and format
    express_rank_DF <- express_rank_DF[ ,-1]

    return(express_rank_DF)
}

Main_Function <- function(){

  # Cluster annotations
  cluster_annot <- c(
    "9" = "vRG"
    , "7" = "oRG"
    , "8" = "Cycling progenitor S phase"
    , "10" = "Cycling progenitor G2/M phase"
    , "2" = "IPC"
    , "0" = "Excitatory neuron new born migrating"
    , "1" = "Excitatory neuron"
    , "4" = "Excitatory neuron (callosal)"
    , "3" = "Deep layer excitatory neuron 1"
    , "13" = "Deep layer excitatory neuron 2"
    , "5" = "Interneuron (SST)"
    , "6" = "Interneuron (CALB2)"
    , "11" = "Oligodendrocyte precursor"
    , "12" = "Endothelial"
    , "14" = "Pericyte"
    , "15" = "Microglia"
  )

  # Determine expressed genes and rank them for each cluster
  express_rank_DFL <- lapply(c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  , function(cluster){
    express_rank_DF <- Rank_Genes_By_Expression_For_Cluster(
      cluster = cluster
      , cluster_annot = cluster_annot
      , seuratO = centSO
    )
    # Add empty column for formatting
    # express_rank_DF <- cbind(express_rank_DF, "")
  })

  # Combine into one matrix
  express_rank_M <- matrix("", ncol = 5*length(express_rank_DFL)
    , nrow = max(sapply(express_rank_DFL, nrow))
  )
  for(i in 1:length(express_rank_DFL)){
    mt1 <- as.matrix(express_rank_DFL[[i]])
    col_idx <- ((i-1)*5)
    express_rank_M[1:nrow(mt1), 1:ncol(mt1)+col_idx] <-
      mt1[1:nrow(mt1), 1:ncol(mt1)]
  }

  # Add colnames
  colnames(express_rank_M) <- rep(
    c(colnames(express_rank_DFL[[1]]), "")
    , length(express_rank_DFL)
  )

  # Out
  write.csv(express_rank_M, file = paste0(outTable, "10percent.csv")
    , quote = FALSE, row.names = FALSE)
}

Main_Function()
################################################################################
