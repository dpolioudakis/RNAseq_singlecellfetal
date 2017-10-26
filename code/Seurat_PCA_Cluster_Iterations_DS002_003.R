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
graphCodeTitle <- "Seurat_PCA_Cluster_Iterations_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Seurat_PCA_Cluster_Iterations_DS002_003_exon_FtMm250_"
outTablePfx <- "../analysis/tables/Seurat_PCA_Cluster_Iterations_DS002_003_exon_FtMm250_"
outRdat <- "../analysis/Seurat_PCA_Cluster_Iterations_DS002_003_exon_FtMm250.Robj"

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

## PCA all genes

# Run PCA with the IRLBA package (iteratively computes the top dimensions
# , dramatic increase in speed since we are throwing away most PCs anyway)
pcAllGenesSO <- PCAFast(seuratO, pc.genes = unlist(seuratO@data@Dimnames[1])
  , pcs.compute = 40, pcs.print = 10)

## PCA Plots

## Plot PCA graph colored by VZ or CP samples
# Collect PCA values
ggDF <- pcAllGenesSO@pca.rot
# Format sample names for ggplot2
ggDF$SAMPLE <- gsub(".*_", "", rownames(ggDF))
ggDF$SAMPLE[ggDF$SAMPLE == "CPH1" | ggDF$SAMPLE == "CPH2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "CPS1" | ggDF$SAMPLE == "CPS2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "VZH1" | ggDF$SAMPLE == "VZH2"] <- "Ventricular zone"
ggDF$SAMPLE[ggDF$SAMPLE == "VZS1" | ggDF$SAMPLE == "VZS2"] <- "Ventricular zone"
# ggplot
ggplot(data = ggDF, aes(x = PC1, y = PC2, col = SAMPLE)) +
  geom_point(shape = 21, alpha = 0.5)
ggsave(paste0(outGraphPfx, "PCAplots_Region_AllGenes.pdf"), width = 6, height = 5)

## Plot PCA graph colored by VZ or CP samples
# Collect PCA values
ggDF <- seuratO@pca.rot
# Format sample names for ggplot2
ggDF$SAMPLE <- gsub(".*_", "", rownames(ggDF))
ggDF$SAMPLE[ggDF$SAMPLE == "CPH1" | ggDF$SAMPLE == "CPH2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "CPS1" | ggDF$SAMPLE == "CPS2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "VZH1" | ggDF$SAMPLE == "VZH2"] <- "Ventricular zone"
ggDF$SAMPLE[ggDF$SAMPLE == "VZS1" | ggDF$SAMPLE == "VZS2"] <- "Ventricular zone"
# ggplot
ggplot(data = ggDF, aes(x = PC1, y = PC2, col = SAMPLE)) +
  geom_point(shape = 21, alpha = 0.5)
ggsave(paste0(outGraphPfx, "PCAplots_Region_VarGenes.pdf"), width = 6
  , height = 5)

## Plot PCA graph colored by Seurat clustering - variable genes
ggO <- PCAPlot(seuratO, 1, 2, pt.size = 0.5, do.return = TRUE)
ggsave(paste0(outGraphPfx, "PCAplots_Cluster_VarGenes.pdf"), plot = ggO
  , width = 5, height = 5)

## Plot PCA graph colored by Seurat clustering - all genes
ggO <- PCAPlot(pcAllGenesSO, 1, 2, pt.size = 0.5, do.return = TRUE)
ggsave(paste0(outGraphPfx, "PCAplots_Cluster_AllGenes.pdf"), plot = ggO
  , width = 5, height = 5)

# PC elbow plot - all genes
pdf(paste0(outGraphPfx, "PCElbowPlot_AllGenes.pdf"))
PCElbowPlot(pcAllGenesSO, num.pc = 40)
dev.off()

# PC elbow plot - variable genes
pdf(paste0(outGraphPfx, "PCElbowPlot_VarGenes.pdf"))
PCElbowPlot(seuratO, num.pc = 40)
dev.off()

# # Jackstraw plot - all genes
# pdf(paste0(outGraphPfx, "JackStraw_AllGenes.pdf"))
# seuratO <- JackStraw(seuratO, num.replicate = 100, do.print = FALSE)
# JackStrawPlot(pcAllGenesSO, PCs = 1:20)
# dev.off()
# 
# # Jackstraw plot - variable genes
# pdf(paste0(outGraphPfx, "JackStraw_VarGenes.pdf"))
# seuratO <- JackStraw(seuratO, num.replicate = 100, do.print = FALSE)
# JackStrawPlot(seuratO, PCs = 1:20)
# dev.off()
################################################################################

### Cluster and tSNE

# Test clustering using resolutions 0.1 to 4.0 on PC1-10
resIdx <- seq(0.1, 4.0, 0.1)
# Variable genes
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_PCA1to10_VarGenes_tSNE_ResTest.pdf"))
lapply(resIdx, function(i) {
  print(i)
  res <- i
  seuratO <- FindClusters(seuratO, pc.use = 1:10, resolution = res, print.output = 0
    , save.SNN = T)
  seuratO <- RunTSNE(seuratO, dims.use = 1:10, do.fast = T)
  ggO <- TSNEPlot(seuratO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
  ggO + ggtitle(paste0("Resolution: ", i))
})
dev.off()
# Test clustering using resolutions 0.1 to 4.0 on PC1-15 variable genes
resIdx <- seq(0.1, 4.0, 0.1)
# Variable genes
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_PCA1to10_VarGenes_tSNE_ResTest.pdf"))
lapply(resIdx, function(i) {
  print(i)
  res <- i
  seuratO <- FindClusters(seuratO, pc.use = 1:15, resolution = res
    , print.output = 0, save.SNN = T)
  seuratO <- RunTSNE(seuratO, dims.use = 1:10, do.fast = T)
  ggO <- TSNEPlot(seuratO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
  ggO + ggtitle(paste0("Resolution: ", i))
})
dev.off()
# Test clustering using resolutions 0.1 to 4.0 on PC1-15 all genes
resIdx <- seq(0.1, 4.0, 0.1)
# Variable genes
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_PCA1to10_AllGenes_tSNE_ResTest.pdf"))
lapply(resIdx, function(i) {
  print(i)
  res <- i
  pcAllGenesSO <- FindClusters(pcAllGenesSO, pc.use = 1:15, resolution = res
    , print.output = 0, save.SNN = T)
  pcAllGenesSO <- RunTSNE(pcAllGenesSO, dims.use = 1:15, do.fast = T)
  ggO <- TSNEPlot(pcAllGenesSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
  ggO + ggtitle(paste0("Resolution: ", i))
})
dev.off()

# Test using PCs 1:2 to 1:20
pcIdx <- seq(2,20)
# Variable genes
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_VarGenes_Res0.6_tSNE_PCAtest.pdf"))
lapply(pcIdx, function(i) {
  print(i)
  pcs <- c(1:i)
  seuratO <- FindClusters(seuratO, pc.use = pcs, resolution = 0.6, print.output = 0
    , save.SNN = T)
  seuratO <- RunTSNE(seuratO, dims.use = pcs, do.fast = T)
  ggO <- TSNEPlot(seuratO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
  ggO + ggtitle(paste0("PCs: 1_", i[length(i)]))
})
dev.off()
# All genes
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_AllGenes_Res0.6_tSNE_PCAtest.pdf"))
lapply(pcIdx, function(i) {
  print(i)
  pcs <- c(1:i)
  pcAllGenesSO <- FindClusters(pcAllGenesSO, pc.use = pcs, resolution = 0.6, print.output = 0
    , save.SNN = T)
  pcAllGenesSO <- RunTSNE(pcAllGenesSO, dims.use = pcs, do.fast = T)
  ggO <- TSNEPlot(pcAllGenesSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
  ggO + ggtitle(paste0("PCs: 1_", i[length(i)]))
})
dev.off()


# Jackstraw plot - all genes
pdf(paste0(outGraphPfx, "JackStraw_AllGenes.pdf"))
seuratO <- JackStraw(seuratO, num.replicate = 100, do.print = FALSE)
JackStrawPlot(pcAllGenesSO, PCs = 1:20)
dev.off()

# Jackstraw plot - variable genes
pdf(paste0(outGraphPfx, "JackStraw_VarGenes.pdf"))
seuratO <- JackStraw(seuratO, num.replicate = 100, do.print = FALSE)
JackStrawPlot(seuratO, PCs = 1:20)
dev.off()


# Test tSNE on scaled expression values, not PCs
pdf(paste0(outGraphPfx, "200gd_2500dg_Mt5_AllGenesNoPC_Res0.6_tSNE.pdf"))
seuratO <- RunTSNE(seuratO, do.fast = T, genes.use = seuratO@scale.data)
TSNEPlot(pcAllGenesSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
dev.off()
################################################################################