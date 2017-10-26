
# Damon Polioudakis
# 2017-08-15
# Run monocle
################################################################################

# Error occured after including "ERC1", "ERCC1"
# No issue with "ERC1"


rm(list = ls())

require(Seurat)
require(biomaRt)
source("2016-11-03-runMetaNeighbor.R")
# source("2017-08-28-runMN-US.R")
# options(stringsAsFactors = FALSE)

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
taskID <- as.numeric(args[1])

# load("dev/MetaNeighbor_sample_data.Rdata")

load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# biomart
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv", header = TRUE)

## Variables
graphCodeTitle <- "MetaNeighbor_DS2-11.R"
outGraph <- "../analysis/graphs/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_"
outTable <- "../analysis/tables/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_"
outRdat <- "../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdat)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Test MetaNeighbor

# load("dev/MetaNeighbor_sample_data.Rdata")
# # rnames <- row.names(data)
# # cnames <- colnames(data)
# # data <- exM[1:3157, 1:1051]
# # # row.names(data) <- rnames
# # colnames(data) <- cnames
# # data[1:5, 1:5]
# 
# cell.lab <- cbind(cell.lab, cell.lab)
# 
# # Run MetaNeighbor
# AUROC.scores = run_MetaNeighbor(data = data
#   , experiment_labels = exp.lab
#   , celltype_labels = cell.lab
#   , genesets = genesets
#   , file_ext = paste0(outRdat, "Results.rdat"))
# 
# load("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_Results.rdat.IDscore.list.Rdata")
# load("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_Results.rdat.IDscore.matrix.Rdata")
# ROCs
# nv.mat
################################################################################

### MetaNeighbor

exM <- as.matrix(centSO@data)
mns <- rowMeans(exM)
idx <- sort(mns, decreasing = TRUE)[1:10000]
exM <- exM[row.names(exM) %in% names(idx), ]

# Convert hgnc symbols to ensembl
idx <- match(row.names(exM), bmDF$hgnc_symbol)
ens <- bmDF$ensembl_gene_id[idx]
row.names(exM)[! is.na(ens)] <- as.character(ens[! is.na(ens)])

# Subset raw counts to filtered cells and genes
# exM <- centSO@raw.data[row.names(centSO@raw.data) %in% row.names(centSO@scale.data)
#   , colnames(centSO@raw.data) %in% colnames(centSO@scale.data)]

# # Subset to genes that are in biomart
# genes <- row.names(exM)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# genesBmDF <- getBM(attributes = c("hgnc_symbol")
#   , filters = "hgnc_symbol",
#   values = genes, mart = mart)
# exM <- exM[row.names(exM) %in% genesBmDF$hgnc_symbol, ]

# Genes to look up GO for
genes <- row.names(exM)

# GO ontologies from biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
goDF <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id")
  , filters = "ensembl_gene_id",
  values = genes, mart = mart)
# Format GO for MetaNeighbor
ldf <- split(goDF, goDF$go_id)
genesGoL <- lapply(ldf, function(goDF) goDF$ensembl_gene_id)
# Filter ontologies that only have one gene in dataset (breaks MetaNeighbor function)
idx <- sapply(genesGoL, function(go) {
  sum(row.names(exM) %in% unlist(go)) > 1
})
genesGoL <- genesGoL[idx]
# Subset ontologies to those with > 5 genes
v1 <- sapply(genesGoL, length)
idx <- v1 > 10
genesGoL <- genesGoL[idx]

# Format cluster IDs for MetaNeighbor
celltype_labels <- matrix( , nrow = length(names(centSO@ident))
  , ncol = length(unique(centSO@ident)))
for(column in 1:ncol(celltype_labels)) {
  celltype_labels[ ,column] <- 0
  clusterID <- sort(unique(centSO@ident))[column]
  celltype_labels[centSO@ident == clusterID, column] <- 1
}

# celltype_labels <- as.matrix(as.numeric(as.character(centSO@ident)))
row.names(celltype_labels) <- names(centSO@ident)

# Run MetaNeighbor
AUROC.scores = run_MetaNeighbor(data = exM
  , experiment_labels = as.character(centSO@meta.data$BRAIN)
  , celltype_labels = celltype_labels
  , genesets = genesGoL[taskID]
  , file_ext = paste0(outRdat, "Results", taskID, ".rdat"))

save(AUROC.scores, file = paste0(outRdat, "AUROCscores", taskID, ".rdat"))
# # 
# load("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_ResultsNA.rdat.IDscore.list.Rdata")
# load("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/MetaNeighbor_DS2-11_ResultsNA.rdat.IDscore.matrix.Rdata")
# print(ROCs)
# print(nv.mat)
################################################################################

# inROCs <- list.files("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/"
#   , full.names = TRUE)
# inROCs <- inROCs[grep("list", inROCs)]
# lapply(inROCs, function(inROC){
#   load(inROC)
#   print(ROCs)
# })
# 
# inNVMs <- list.files("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/"
#   , full.names = TRUE)
# inNVMs <- inNVMs[grep("matrix", inNVMs)]
# lm <- lapply(inNVMs, function(inNVM){
#   load(inNVM)
#   print(nv.mat)
# })
# m1 <- do.call("rbind", lm)
# 
# inAUROCs <- list.files("../analysis/MetaNeighbor_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40/"
#   , full.names = TRUE)
# inAUROCs <- inAUROCs[grep("AUROC", inAUROCs)]
# lapply(inAUROCs, function(inAUROC){
#   load(inAUROC)
#   print(AUROC.scores)
# })
# 
# 
# 
# 
