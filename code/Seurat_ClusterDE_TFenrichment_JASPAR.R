# Damon Polioudakis
# 2017-10-03
# Submit list of genes to Vijays TF binding site enrichment pipeline JASPAR
################################################################################

rm(list = ls())

require(RegFacEnc)

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
clusterID <- as.character((as.numeric(args[1]) - 1))
print(paste0("Cluster ID: ", clusterID))

setwd("/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/code")

## Inputs
# Seurat cluster DE
deDF <- read.table(
  "../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

## Output Directories
dir.create("../analysis/Seurat_ClusterDE_TFenrichment/JASPAR", recursive = TRUE)
################################################################################

### Format DE file for TF enrichment

# Split DE table by cluster
deLDF <- split(deDF, deDF$CLUSTER)

# Check number of genes above log FC in each cluster
v1 <- sapply(deLDF, function(df) table(df$LOG_FC > 0.4))
# Format into DF
nFiltDF <- data.frame(t(data.frame(v1)))
nFiltDF$CLUSTER <- row.names(nFiltDF)
nFiltDF$NUMBER_DE_GENES <- nFiltDF$TRUE.
nFiltDF <- nFiltDF[, -c(1:2)]
write.csv(nFiltDF, "../analysis/tables/Seurat_ClusterDE_TFenrichment_logFC04.csv"
  , quote = FALSE)

# Subset to cluster
df <- deLDF[[clusterID]]

# Filter by fold change
df <- df[df$LOG_FC > 0.4, ]

# Format and save as csv
df <- data.frame(Genename = df$GENE, EnsemblID = df$ENSEMBL)
df$Sno <- c(1:nrow(df))
df <- df[c("Sno", "Genename", "EnsemblID")]

# Two output directories, one to run JASPAR, one TRANSFAC because Vijay's tool
# doesn't allow to specify output directory
write.csv(df, file = paste0(
  "../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/Seurat_ClusterDE_TFenrichment_logFC04_Cluster"
  , clusterID, ".csv"), row.names = FALSE)
################################################################################

### Run TF enrichment

# JASPAR

print("Running for JASPAR...")

setwd("/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/analysis/Seurat_ClusterDE_TFenrichment/JASPAR")

inFile <- paste0("Seurat_ClusterDE_TFenrichment_logFC04_Cluster", clusterID, ".csv")

TFBSenrich(user.file = inFile
  , TF.db = system.file("data/JASPAR_CLOVER", package = "RegFacEnc")
  , TF.nome = system.file(
    "data/JASPAR_NOMENCLATURE_TABLE", package = "RegFacEnc")
  , species = "Human", TF_motifs = "JASPAR", BF_Type = "Protein Coding"
  , pval = "0.05", option = "-t")

print("Done running for JASPAR...")
################################################################################

