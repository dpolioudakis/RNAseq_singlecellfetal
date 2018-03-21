# Damon Polioudakis
# 2017-08-28
# DE genes for each cluster

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
# require(irlba)
require(gridExtra)
require(fdrtool)
require(cowplot)
source("Function_Library.R")

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
clusterID <- (as.numeric(args[1]) - 1)
print(paste0("Cluster ID: ", clusterID))

## Inputs
# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM
# Biomart to add ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterDE.R"
# Output paths
# Sub path
out_sub_path <- paste0(
  "Seurat_ClusterDE_DS2-11/"
  ,"FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/"
  , "res054/"
  , "Seurat_ClusterDE_DS2-11_"
)
outGraph <- paste0("../analysis/graphs/", out_sub_path)
outTable <- paste0("../analysis/tables/", out_sub_path)
outData <- paste0("../analysis/processed_data/", out_sub_path)

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 14)))
theme_update(plot.title = element_text(size = 14))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

### Differentially expressed genes for each cluster versus all other cells

print(paste0(
  "### Finding differentially expressed genes for cluster ", clusterID)
)

# # Filter cells
# df <- DE_Filters_ExpMatrix(centSO
#   , minPercent = 10, foldChange = 0.2, clusterID = clusterID)

# DE Linear model
termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$cluster <- FALSE
termsDF$cluster[centSO@ident == clusterID] <- TRUE
mod <- "y ~ cluster+nUMI+librarylab+individual"
deLM <- DE_Linear_Model(exDatDF = centSO@data, termsDF = termsDF, mod = mod)

# Format LM output into data frame
deDF <- Format_DE(deLM, centSO, clusterID)

# FDR correct
# NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
deDF$FDR <- p.adjust(deDF$Pvalue, method = "BH")
# Check
table(deDF$Pvalue < 0.05)
table(deDF$FDR < 0.05)

# Write to tab delimited table
print("Writing DE table to text file")
write.table(x = deDF
  , file = paste0(outTable, "Cluster", clusterID, "_Vs_All_Clusters.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)
################################################################################
