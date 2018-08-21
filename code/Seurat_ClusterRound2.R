# Damon Polioudakis
# 2017-05-28
# 2nd iteration of Seurat clustering

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
require(viridis)
source("Function_Library.R")
# require(xlsx)

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)

## Inputs

# Seurat clustering object
# PC 1-40
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Cell cycle markers from Macosko 2015 Table S2 to remove from variable gene
# list used for clustering
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterRound2.R"
outGraph <- "../analysis/graphs/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_"
outTable <- "../analysis/tables/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_"
outData <- "../analysis/analyzed_data/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_"
# outGraph <- "../analysis/graphs/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"
# outTable <- "../analysis/tables/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"
# outData <- "../analysis/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"

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

# load(paste0(outData, "seuratO.Robj"))
################################################################################

### Split by cluster

## Select Seurat cluster to re-cluster
# List of cluster IDs
# clustersL <- append(list(c(0, 1, 4), c(0, 1), c(3, 13), c(5, 6), c(7, 9))
#   , sort(as.numeric(as.character(unique(centSO@ident)))))
clustersL <- sort(as.numeric(as.character(unique(centSO@ident))))
clusterID <- clustersL[[as.numeric(args[1])]]
print(paste0("Cluster ID: ", clusterID))

Subset_Seurat_Raw_Counts_By_Cluster <- function(clusterID, seuratO){
  print("Subset_Seurat_Raw_Counts_By_Cluster")
  exDF <- centSO@raw.data
  clusterIDs_cellIDs <- centSO@ident
  clusterIDs_cellIDs <- clusterIDs_cellIDs[clusterIDs_cellIDs %in% clusterID]
  exDF <- exDF[, colnames(exDF) %in% names(clusterIDs_cellIDs)]
  return(exDF)
}

Create_Subset_Seurat_Object <- function(seuratO, clusterID){
  print("Create_Subset_Seurat_Object")
  exDF <- Subset_Seurat_Raw_Counts_By_Cluster(
    clusterID = clusterID, seuratO = seuratO
  )
  so <- CreateSeuratObject(raw.data = exDF
    , min.cells = 0, min.genes = 0
    , normalization.method = "LogNormalize", scale.factor = 10000
    , project = paste0("Cluster ", clusterID) , do.scale = FALSE, do.center = FALSE)
  # Add metadata
  metDF <- centSO@meta.data
  ssMetDF <- metDF[row.names(metDF) %in% colnames(so@raw.data), ]
  so <- AddMetaData(so, metadata = ssMetDF)
  return(so)
}

# Initialize Seurat object for each cluster
so <- Create_Subset_Seurat_Object(seuratO = centSO, clusterID = clusterID)

# Save centered scaled expression from Seurat round 1
rd1CentExM <- centSO@scale.data[
  , colnames(centSO@scale.data) %in% colnames(so@raw.data)]
idx <- match(colnames(so@raw.data), colnames(rd1CentExM))
rd1CentExM <- rd1CentExM[ ,idx]
################################################################################

### Detection of variable genes across the single cells

print("### Detection of variable genes across the single cells")

so <- FindVariableGenes(so, mean.function = ExpMean
    , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3
    , y.cutoff = 0.5
    )

# Check length of var.genes
print(paste0("Cluster: ", clusterID, " length var.genes: "
 , length(so@var.genes)
))
################################################################################

### Regress out unwanted sources of variation

print("### Regress out unwanted sources of variation")

# Your single cell dataset likely contains ‘uninteresting’ sources of variation.
# This could include not only technical noise, but batch effects, or even
# biological sources of variation (cell cycle stage). As suggested in Buettner
# et al, NBT, 2015, regressing these signals out of the analysis can improve
# downstream dimensionality reduction and clustering. To mitigate the effect of
# these signals, Seurat constructs linear models to predict gene expression
# based on user-defined variables. The scaled z-scored residuals of these models
# are stored in the scale.data slot, and are used for dimensionality reduction
# and clustering.

# note that this overwrites so@scale.data. Therefore, if you intend to use
# ScaleData, you can set do.scale=F and do.center=F in the original object to
# save some time.

## Regress covariates and center scale

so <- ScaleData(so, vars.to.regress = c("nUMI", "librarylab", "individual"))

## No center or scale

# Make dataframe of covariates to regress out
covDF <- data.frame(nUMI = so@meta.data$nUMI
  , librarylab = so@meta.data$librarylab
  , individual = so@meta.data$individual)
# Regress out confounding variables
RegressCovariates <- function (exM, covDF) {
  exRegCovM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
  rownames(exRegCovM) <- rownames(exM)
  colnames(exRegCovM) <- colnames(exM)
  # ncol(covDF)+1 when condition has 2 levels
  coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
  for (i in 1:nrow(exM)) {
    if (i%%1000 == 0) {print(paste("Done for", i, "genes..."))}
    mod <- lm(as.numeric(exM[i, ]) ~ ., data = covDF)
    # The full data - the undesired covariates
    exRegCovM[i,] <- coef(mod)[1] + mod$residuals
    # lmmod1 <- lm(as.numeric(exM[i, ]) ~ condition + age + sex + pmi, data = covDF)
  }
  return(exRegCovM)
}

noCentExM <- RegressCovariates(so@data, covDF)
################################################################################

### Perform linear dimensional reduction

print("### Perform linear dimensional reduction")

# Perform PCA on the scaled data. By default, the genes in object\@var.genes are
# used as input, but can be defined using pc.genes. We have typically found that
# running dimensionality reduction on genes with high-dispersion can improve
# performance. However, with UMI data - particularly after using ScaleData, we
# often see that PCA returns similar (albeit slower) results when run on much
# larger subsets of genes, including the whole transcriptome.
# r2SO <- PCA(r2SO, pc.genes = r2SO@var.genes, do.print = TRUE
#   , pcs.print = 5, genes.print = 5)

# Run PCA with the IRLBA package (iteratively computes the top dimensions,
# dramatic increase in speed since we only use a fraction of the PCs anyways) if
# you see the warning "did not converge–results might be invalid!; try
# increasing maxit or fastpath=FALSE", try increasing maxit
so <- RunPCA(object = so, pc.genes = so@var.genes, pcs.compute = 50
  , pcs.print = 1:2, genes.print = 5, maxit = 500, weight.by.var = FALSE)

# ProjectPCA scores each gene in the dataset (including genes not included in
# the PCA) based on their correlation with the calculated components. Though we
# don't use this further here, it can be used to identify markers that are
# strongly correlated with cellular heterogeneity, but may not have passed
# through variable gene selection. The results of the projected PCA can be
# explored by setting use.full=T in the functions below.
# Saves in slot @pca.x.full
so <- ProjectPCA(so, pcs.print = 0)

# Plot genes with highest PC loadings
ggL <- lapply(1:8, function(pc) {
  df <- rbind(data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc])[1:10])
  , data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc], decreasing = TRUE)[1:10])
    )
  df$GENE <- factor(row.names(df), levels = row.names(df))
  gg <- ggplot(df, aes(x = PC, y = GENE)) +
    geom_point() +
    xlab("Loading") +
    ylab("Gene") +
    ggtitle(paste0("Cluster: ", clusterID, "\nPC ", pc))
  return(gg)
})
Plot_Grid(ggL, ncol = 4, align = 'v', axis = 'r', rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nGenes with highest PC loadings"
    , "\nCluster: ", clusterID)
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_PCAplots.pdf"), width = 13
  , height = 14, limitsize = FALSE)
################################################################################

### Determine statistically significant principal components

print("### Determine statistically significant principal components")

# PC variance plot
gg <- DimElbowPlot(so, reduction.type = "pca", dims.plot = 50)
gg <- gg + ggtitle(
    paste0(graphCodeTitle
      , "\n\nPC variance for each cluster"
      , "\nCluster: ", clusterID)
    ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0, max(so@dr$pca@sdev)))
ggsave(paste0(outGraph, "Cluster", clusterID, "_PCElbowPlot.pdf")
  , width = 5, height = 5)

save(so, noCentExM, rd1CentExM
  , file = paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))
################################################################################

### Run Non-linear dimensional reduction (tSNE)

print("### Run Non-linear dimensional reduction (tSNE)")

## PCs test
Seurat_PCs_Test <- function(seuratO, resolution, dims_use, title){
  print("Seurat_PCs_Test")
  ggL <- lapply(dims_use, function(dims_use) {
    tryCatch({
      print(paste0("Dimensions: ", dims_use[1], "-", dims_use[length(dims_use)]))
      seuratO <- RunTSNE(seuratO, dims.use = dims_use, do.fast = TRUE)
      seuratO <- FindClusters(
        seuratO, dims.use = dims_use, resolution = resolution
        , print.output = 0, save.SNN = TRUE)
      gg <- TSNEPlot(seuratO, pt.size = 0.01, do.return = TRUE, do.label = TRUE)
      gg <- gg + ggplot_set_theme_publication
      gg <- gg + ggtitle(paste0(
        "PCs: ", dims_use[1], "-", dims_use[length(dims_use)]
      ))
      return(gg)
      } ,
      error = function(cond) {
        print(paste0("Error for dimensions: ", dims_use[1], "-", dims_use[length(dims_use)]))
        message(cond)
        # Choose a return value in case of error
        return(NA)
      }
    )

  })
  return(ggL)
}
# Plot
ggL <- Seurat_PCs_Test(
  seuratO = so
  , title = name
  , dims_use = list(c(1:5), c(1:10), c(1:15), c(1:20), c(1:30), c(1:40))
  # , dims_use = list(c(1:10), c(1:20))
  , resolution = 0.4
)
ggL <- ggL[! is.na(ggL)]
# Combine tSNE graphs
pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different PCs used"
    , "\nCluster: ", clusterID)
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PCsTest.png")
  , width = 12, height = 7)


## Resolution test
Seurat_Resolution_Test <- function(seuratO, resolutions, dims_use){
  print("Seurat_Resolution_Test")
  so <- tryCatch({
    RunTSNE(seuratO, dims.use = dims_use, do.fast = TRUE)
    } ,
    error = function(cond) {
      print(paste0("Error for dimensions: ", dims_use))
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
  ggL <- lapply(resolutions, function(resolution) {
    tryCatch({
      print(paste0("Resolution: ", resolution))
      so <- FindClusters(so, dims.use = dims_use, resolution = resolution
        , print.output = 0, save.SNN = TRUE)
      gg <- TSNEPlot(so, pt.size = 0.01, do.return = TRUE, do.label = TRUE)
      gg <- gg + ggtitle(paste0("Resolution: ", resolution))
      gg <- gg + ggplot_set_theme_publication
      return(gg)
    } ,
    error = function(cond) {
      print(paste0("Error for resolution: ", resolution))
      message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
  })
  return(ggL)
}
# Plot using PC1-40
ggL <- Seurat_Resolution_Test(
  seuratO = so
  , dims_use = 1:40
  # , resolutions = c(0.1, 0.2)
  , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
)
ggL <- ggL[! is.na(ggL)]
# Combine tSNE graphs
pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different resolutions"
    , "\nCluster: ", clusterID
    , "\nPC 1-40"
  )
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-40_ResolutionTest.png")
  , width = 12, height = 7)
# Plot using PC1-10
ggL <- Seurat_Resolution_Test(
  seuratO = so
  , dims_use = 1:10
  # , resolutions = c(0.1, 0.2)
  , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
)
ggL <- ggL[! is.na(ggL)]
# Combine tSNE graphs
pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different resolutions"
    , "\nCluster: ", clusterID
    , "\nPC 1-10"
  )
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-10_ResolutionTest.png")
  , width = 12, height = 7)

# Plot using PC1-15
ggL <- Seurat_Resolution_Test(
  seuratO = so
  , dims_use = 1:15
  # , resolutions = c(0.1, 0.2)
  , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
)
ggL <- ggL[! is.na(ggL)]
# Combine tSNE graphs
pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different resolutions"
    , "\nCluster: ", clusterID
    , "\nPC 1-10"
  )
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-15_ResolutionTest.png")
  , width = 12, height = 7)

# Plot using PC1-5
ggL <- Seurat_Resolution_Test(
  seuratO = so
  , dims_use = 1:5
  # , resolutions = c(0.1, 0.2)
  , resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
)
ggL <- ggL[! is.na(ggL)]
# Combine tSNE graphs
pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE with different resolutions"
    , "\nCluster: ", clusterID
    , "\nPC 1-5"
  )
)
ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-5_ResolutionTest.png")
  , width = 12, height = 7)

## Clustering parameters to save
so <- RunTSNE(so, dims.use = 1:10, do.fast = TRUE)
so <- FindClusters(so, dims.use = 1:10, resolution = 0.5
  , print.output = 0, save.SNN = TRUE)

save(so, noCentExM, rd1CentExM
  , file = paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))
################################################################################

print("### Feature plot of covariates")

Feature_Plot_Of_Covariate <- function(seuratO, covariate){
  ("Feature_Plot_Of_Covariate")
  # Collect tSNE values
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  # Add cluster identity
  ggDF$CLUSTER <- so@ident
  # Add metadata
  ggDF <- data.frame(ggDF, so@meta.data[match(row.names(ggDF), so@meta.data$CELL), ])
  ggDF$BRAIN <- as.factor(ggDF$BRAIN)
  # Plot
  gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF[[covariate]])) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggplot_set_theme_publication +
    ggtitle(covariate)
  if(class(ggDF[[covariate]]) == "numeric"){gg <- gg + scale_color_viridis()}
  return(gg)
}

## PC 1-10
# tSNE colored by clustering
ggTsne <- TSNE_Plot(so) +
  theme(legend.position = "none") +
  ggplot_set_theme_publication
# tSNE colored by covariate
gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
# plot grid
pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
  , align = "v", axis = "r", ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat clustering round 2 and covariates"
    , "\nCluster: ", clusterID
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-10 resolution 0.5"
    , "\n")
  )
ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-10_covariates.png")
  , width = 14, height = 7, limitsize = FALSE)

## PC 1-15
# Re do tSNE and clustering
so <- RunTSNE(so, dims.use = 1:15, do.fast = TRUE)
so <- FindClusters(so, dims.use = 1:15, resolution = 0.5
, print.output = 0, save.SNN = TRUE)
# tSNE colored by clustering
ggTsne <- TSNE_Plot(so) +
  theme(legend.position = "none") +
  ggplot_set_theme_publication
# tSNE colored by covariate
gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
# plot grid
pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
  , align = "v", axis = "r", ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat clustering round 2 and covariates"
    , "\nCluster: ", clusterID
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-15 resolution 0.5"
    , "\n")
  )
ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-15_covariates.png")
  , width = 14, height = 7, limitsize = FALSE)

## PC 1-5
# Re do tSNE and clustering
so <- RunTSNE(so, dims.use = 1:5, do.fast = TRUE)
so <- FindClusters(so, dims.use = 1:5, resolution = 0.5
, print.output = 0, save.SNN = TRUE)
# tSNE colored by clustering
ggTsne <- TSNE_Plot(so) +
  theme(legend.position = "none") +
  ggplot_set_theme_publication
# tSNE colored by covariate
gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
# plot grid
pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
  , align = "v", axis = "r", ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat clustering round 2 and covariates"
    , "\nCluster: ", clusterID
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-5 resolution 0.5"
    , "\n")
  )
ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-5_covariates.png")
  , width = 14, height = 7, limitsize = FALSE)

## PC 1-40
# Re do tSNE and clustering
so <- RunTSNE(so, dims.use = 1:40, do.fast = TRUE)
so <- FindClusters(so, dims.use = 1:40, resolution = 0.5
, print.output = 0, save.SNN = TRUE)
# tSNE colored by clustering
ggTsne <- TSNE_Plot(so) +
  theme(legend.position = "none") +
  ggplot_set_theme_publication
# tSNE colored by covariate
gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
# plot grid
pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
  , align = "v", axis = "r", ncol = 3, rel_height = 0.3
  , title = paste0(graphCodeTitle
    , "\n\nSeurat clustering round 2 and covariates"
    , "\nCluster: ", clusterID
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-40 resolution 0.5"
    , "\n")
  )
ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-40_covariates.png")
  , width = 14, height = 7, limitsize = FALSE)
################################################################################

### Color 40k cell tSNE by sub-clustering

# Collect tSNE values
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# Add round 1 cluster identity
ggDF$Cluster <- centSO@ident
# Add round 2 cluster identity
idx <- match(rownames(ggDF), names(so@ident))
ggDF$Cluster_Round2 <- so@ident[idx]
# Plot
gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster)) +
  geom_point(size = 0.02, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2
    , title = "Cluster")) +
  ggplot_set_theme_publication +
  ggtitle("Clustering round 1")
gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster_Round2)) +
  geom_point(size = 0.02, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 3)
    , title = "Cluster")) +
  ggplot_set_theme_publication +
  ggtitle("Clustering round 2")
Plot_Grid(list(gg1, gg2), ncol = 2, rel_height = 0.4, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nSeurat clustering round 2 and covariates"
    , "\nCluster: ", clusterID
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-40 resolution 0.5"
    , "\n")
)
ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_ClusterRound1.png")
  , width = 10, height = 5, limitsize = FALSE)
################################################################################
