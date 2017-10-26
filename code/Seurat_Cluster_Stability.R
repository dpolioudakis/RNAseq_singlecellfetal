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

## Input
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb
rm(fetb)

## Variables
graphCodeTitle <- "Seurat_Cluster_Stability.R"
outGraphPfx <- "../analysis/graphs/Seurat_Cluster_Stability_"
outTablePfx <- "../analysis/tables/Seurat_Cluster_Stability_"
outRdat <- "../analysis/Seurat_Cluster_Stability.Rdat"

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

# # Subset to less cells for testing
# tSO <- seuratO@raw.data[ ,1:500]
# tSO <- new("seurat", raw.data = tSO)
# tSO <- Setup(tSO, min.cells = 3, min.genes = 200, do.logNormalize = T
#   , total.expr = 1e4, project = "10X_PBMC", do.scale = FALSE, do.center = FALSE)
# mito.genes <- grep("^MT-", rownames(tSO@data), value = T)
# percent.mito <- colSums(expm1(tSO@data[mito.genes, ]))/colSums(expm1(tSO@data))
# tSO <- AddMetaData(tSO, percent.mito, "percent.mito")
# tSO <- RegressOut(tSO, latent.vars = c("nUMI", "percent.mito"))
# tSO <- MeanVarPlot(tSO ,fxn.x = expMean, fxn.y = logVarDivMean
#   , x.low.cutoff = -4, x.high.cutoff = 10, y.cutoff = 0, do.contour = F)
# # Seurat PCA
# tSO <- PCA(tSO, pc.genes = tSO@var.genes, do.print = FALSE)
# # Seurat Cluster
# tSO <- FindClusters(tSO, pc.use = 1:10, resolution = 0.6
#   , print.output = 1, save.SNN = TRUE)
# tSO@ident

# Number of bootstraps
nBt <- 20
# Number of clusters in original
nOrCl <- length(unique(seuratO@ident))
# Make matrix to fill in bootstrap jaccard results
jcrdM <- matrix(data = NA, nrow = nOrCl, ncol = nBt)
row.names(jcrdM) <- paste("Cluster", seq(1, nOrCl))
colnames(jcrdM) <- paste("Boot", seq(1, nBt))

# Bootstrap, recluster, + jaccard to test cluster stability
for (i in 1:nBt) {
  print("Bootstrap:", nBt)
  
  # Set seed
  set.seed(i)
  
  # Sample with replacement
  exM <- seuratO@scale.data[ ,sample(ncol(seuratO@scale.data)
    , replace = TRUE)]
  # Unique ID to colnames
  colnames(exM) <- paste0(colnames(exM)
    , "_", 1:length(colnames(exM)))
  
  # Remake seurat object after sampling
  btSO <- new("seurat", raw.data = exM)
  btSO <- Setup(btSO, min.cells = 3, min.genes = 200, do.logNormalize = T
    , total.expr = 1e4, project = "10X_PBMC", do.scale = FALSE, do.center = FALSE)
  # Overwrite expression matrix with bootstrap expression matrix
  btSO@scale.data <- exM
  # Use same variable genes as original clustering
  btSO@var.genes <- seuratO@var.genes
  
  # Seurat PCA
  btSO <- PCA(btSO, pc.genes = btSO@var.genes, do.print = FALSE)
  # Seurat Cluster
  btSO <- FindClusters(btSO, pc.use = 1:10, resolution = 0.6
    , print.output = 0, save.SNN = TRUE)
  btSO@ident
  
  # Intersection of original and bootstrap samples
  int <- intersect(names(seuratO@ident), gsub("*_[0-9]+", "", names(btSO@ident)))
  # Subset original to intersection samples
  clsOr <- seuratO@ident[names(seuratO@ident) %in% int]
  # Remove unique integer IDs added to bootstrap samples
  clsSb <- btSO@ident
  names(clsSb) <- gsub("*_[0-9]+", "", names(btSO@ident))
  # Subset bootstrap to intersection samples
  clsSb <- clsSb[names(clsSb) %in% int]
  # Remove duplicates
  clsSb <- clsSb[unique(names(clsSb))]
  # Merge original and bootstrap to have same order
  clDF <- merge(data.frame(clsOr), data.frame(clsSb), by = "row.names")
  # Make TRUE FALSE list of cluster membership for clujaccard() (fpc package)
  # Original
  origCl <- lapply(unique(clDF$clsOr), function(cl) {clDF$clsOr == cl})
  # Bootstrap
  btClL <- lapply(unique(clDF$clsSb), function(cl) {clDF$clsSb == cl})
  # Max jaccard index for each original cluster and bootstrap clusters
  jcrd <- sapply(origCl, function(oCl) {
    max(sapply(btClL, function(btCl) {clujaccard(oCl, btCl)}))  
  })
  # Add as column of matrix (rows are clusters)
  jcrdM[ ,i] <- jcrd
}
# Reset seed
set.seed(27)

save(jcrdM, file = outRdat)

print(jcrdM)

# Mean jaccard index for each cluster
jcrdMn <- rowMeans(jcrdM)

# Save matrix of jaccard indecies for bootstrap and mean jaccard index
save(jcrdM, jcrdMn, file = outRdat)

# Annotate clusters
an <- c(
  "1 Excitatory Upper Layer Neuron"
  , "2 Excitatory Neuron"
  , "3 Excitatory Upper Layer Neuron"
  , "4 Excitatory Deep Layer Neuron"
  , "5 Intermediate Progenitors"
  , "6 Interneuron"
  , "7 Mitotic Progenitors"
  , "8 oRG"
  , "9 Oligodendrocyte Precursor"
  , "10 Endothelial")
names(jcrdMn) <- an
row.names(jcrdM) <- an

# Format in data frame and output as csv
# Jaccard Coefficient (mean across bootstraps)
jcrdMnDF <- data.frame(CLUSTER = names(jcrdMn), JACCARD_MEAN = jcrdMn) 
jcrdMnM <- as.matrix(jcrdMnDF)
write.csv(jcrdMnM, paste0(outTablePfx, "JaccardMean.csv"), quote = FALSE
  , row.names = FALSE)
# Jaccard index for each boostrap
write.csv(jcrdM, paste0(outTablePfx, "JaccardBoot.csv"), quote = FALSE
  , row.names = TRUE)
################################################################################