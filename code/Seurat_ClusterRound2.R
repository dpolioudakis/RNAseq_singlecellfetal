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
source("Function_Library.R")
# require(xlsx)

args <- commandArgs(trailingOnly = TRUE)
print(args)

## Inputs

# Seurat clustering object
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
centSO <- ssCentSO
noCentExM <- ssNoCentExM

# Cell cycle markers from Macosko 2015 Table S2 to remove from variable gene
# list used for clustering
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterRound2.R"
# outGraph <- "../analysis/graphs/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/RegNumiLibBrain/PC1-30/Seurat_ClusterRound2_DS2-11_"
# outTable <- "../analysis/tables/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/RegNumiLibBrain/PC1-30/Seurat_ClusterRound2_DS2-11_"
# outData <- "../analysis/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/RegNumiLibBrain/PC1-30/Seurat_ClusterRound2_DS2-11_"
outGraph <- "../analysis/graphs/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"
outTable <- "../analysis/tables/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"
outData <- "../analysis/Seurat_ClusterRound2_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/AllGenes/PC1-30/Seurat_ClusterRound2_DS2-11_"

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
# 
# ### Functions
# 
# # Calculate mean expression of group of genes for each cell using seurat scaled
# # expression values
# Mean_Expression <- function(ggDF, genes, exM) {
#   genesExDF <- exM[which(row.names(exM) %in% genes), ]
#   # genesExDF <- as.matrix(seuratO@data)[which(row.names(as.matrix(seuratO@data)) %in% genes), ]
#   # Only calculate column means if there are multiple genes
#   print("Length genes:")
#   print(length(genes))
#   if (is.matrix(genesExDF)) {
#     mnExDF <- colMeans(genesExDF)
#   } else {
#     mnExDF <- genesExDF
#   }
#   # Add to ggplot data frame
#   ggDF$EXPRESSION <- mnExDF[match(row.names(ggDF), names(mnExDF))]
#   return(ggDF)
# }
# 
# # Transform data to desired limits for ggplot2
# Set_Limits <- function(ggDF, limHigh, limLow) {
#   ggDF$EXPRESSION[ggDF$EXPRESSION < limLow] <- limLow
#   ggDF$EXPRESSION[ggDF$EXPRESSION > limHigh] <- limHigh
#   return(ggDF)  
# }
# 
# # Color tSNE plot by expression from Mean_Expression()
# Feature_Plot <- function(ggDF, title, limLow, limHigh) {
#   ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
#     geom_point(size = 0.1) +
#     # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
#     #   , high = "#e31a1c", limits = c(0, 2)) +
#     # scale_color_distiller(name = "Normalized\nexpression", type = "div"
#     #   , palette = 5, direction = -1) +
#     scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
#       , high = "red", limits = c(limLow, limHigh)) +
#     ggtitle(title)
#   return(ggFp)
# }
# 
# # Color tSNE plot by expression from Mean_Expression()
# Feature_Plot_CentScale <- function(ggDF, title, limLow, limHigh) {
#   ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
#     geom_point(size = 0.1) +
#     # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
#     #   , high = "#e31a1c", limits = c(0, 2)) +
#     scale_color_distiller(name = "Normalized\nexpression\nz-score", type = "div"
#       , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
#     # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
#     #   , high = "red", limits = c(limLow, limHigh)) +
#     ggtitle(title)
#   return(ggFp)
# }
# 
# # tSNE plot colored by Seurat clustering
# TSNE_Plot <- function(seuratO) {
#   # tSNE graph colored by cluster
#   ggTsne <- TSNEPlot(seuratO, do.label = TRUE, pt.size = 0.1, do.return = TRUE
#     , no.legend = FALSE)
#   ggTsne <- ggTsne + ggtitle(paste0(
#     "tSNE plot, each point is a cell"
#     , "\nColor indicates cluster assignment"
#   ))
#   ggTsne <- ggTsne + theme_set(theme_bw()) +
#     theme_set(theme_get() + theme(text = element_text(size = 16))) +
#     theme_update(plot.title = element_text(size = 12)) +
#     theme_update(axis.line = element_line(colour = "black")
#       , plot.background = element_blank() 
#       , panel.border = element_blank()
#     )
#   return(ggTsne)
# }
################################################################################

### Split by cluster

exDF <- centSO@raw.data
exDF <- exDF[, colnames(exDF) %in% names(centSO@ident)]
metDF <- centSO@meta.data
clids <- centSO@ident
# # # Subsetting for testing
# idx <- sample(1:ncol(exDF), 3000)
# exDF <- exDF[ ,idx]
# metDF <- metDF[idx, ]
# clids <- clids[idx]
# bak <- noCentExM
# noCentExM <- noCentExM[ ,idx]

# Initialize Seurat object for each cluster
# If cluster is less than 100 cells skip by returning NULL
clustersL <- append(list(c(0, 1, 4, 12), c(0, 1), c(3, 14), c(5, 6), c(7, 9))
  , sort(as.numeric(as.character(unique(centSO@ident)))))
lso <- lapply(clustersL, function(clid){
  print(paste0("Initializing Seurat object for Cluster ID: "
    , paste(clid, collapse = " ")))
  if (sum(clids %in% clid) < 100) {
    print("Less than 100 cells in cluster, skipping...")
    NULL
  } else {
    # Initialize Seurat object
    so <- CreateSeuratObject(raw.data = exDF[ ,clids %in% clid]
      , min.cells = 3, min.genes = 200
      , normalization.method = "LogNormalize", scale.factor = 10000
      , project = "Cluster 0", do.scale = FALSE, do.center = FALSE)
    # Add metadata
    so <- AddMetaData(so, metadata = metDF[clids %in% clid, ])
    # Replace @data with covariates regressed not centered so that mean centered data uses regressed data
    so@data <- noCentExM[row.names(noCentExM) %in% row.names(so@raw.data)
      , clids %in% clid]
    # Mean center, variance scale
    so <- ScaleData(so)
  }
})
# Remove NULLs and change list names of seurat objects to cluster IDs
idx <- sapply(lso, is.null)
lso <- lso[! idx]
names(lso) <- sapply(clustersL, function(cl) paste(cl, collapse = " "))[! idx]
################################################################################

### Detection of variable genes across the single cells

print("### Detection of variable genes across the single cells")

## Using all genes for clustering (set as variable genes)
lso <- lapply(lso, function(so) {
  so@var.genes <- row.names(so@scale.data)
  return(so)
})

# # Seurat calculates highly variable genes and focuses on these for downstream
# # analysis. FindVariableGenes calculates the average expression and dispersion
# # for each gene, places these genes into bins, and then calculates a z-score for
# # dispersion within each bin. This helps control for the relationship between
# # variability and average expression. This function is unchanged from (Macosko
# # et al.), but new methods for variable gene expression identification are
# # coming soon. We suggest that users set these parameters to mark visual
# # outliers on the dispersion plot, but the exact parameter settings may vary
# # based on the data type, heterogeneity in the sample, and normalization
# # strategy. The parameters here identify ~2,000 variable genes, and represent
# # typical parameter settings for UMI data that is normalized to a total of 1e4
# # molecules.
# 
# # Save names
# n <- names(lso)
# # Calculate variable genes
# lso <- lapply(names(lso), function(cl) {
#   print(cl)
#   so <- lso[[cl]]
#   png(paste0(outGraph, "MeanVarPlot_Cluster", cl, ".png"), width = 7, height = 5
#     , units = "in", res = 300)
#   # Was returning kde2d() error which is used for plotting, turned off plotting
#   so <- FindVariableGenes(so, mean.function = ExpMean
#     , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3
#     , y.cutoff = 0.5, num.bin = 20, do.plot = FALSE)
#   # dev.off()
#   return(so)
# })
# print("Number of variable genes used for clustering:")
# sapply(lso, function(so) {length(so@var.genes)})
# # Add names back to list
# names(lso) <- n

# ## Remove CC genes from variable genes
# # Number of variable genes before
# sapply(lso, function(so) {length(so@var.genes)})
# # Remove CC genes
# lso <- lapply(lso, function(so) {
#   ccDF <- data.frame(lapply(ccDF, as.character), stringsAsFactors = FALSE)
#   cc <- c(unlist(ccDF))
#   cc <- gsub(" *", "", cc)
#   idx <- so@var.genes %in% cc
#   so@var.genes <- so@var.genes[! idx]
#   return(so)
# })
# # Number of variable genes after
# sapply(lso, function(so) {length(so@var.genes)})
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

# note that this overwrites centSO@scale.data. Therefore, if you intend to use
# ScaleData, you can set do.scale=F and do.center=F in the original object to
# save some time.

# # Regress covariates and center scale
# lso <- lapply(lso, function(so) {
#   # Regress covariates and center scale
#   so <- ScaleData(so
#     , vars.to.regress = c("nUMI", "librarylab", "individual"))
#   return(so)
# })
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
lso <- lapply(lso, function(so) {
  pcs.compute <- ncol(so@scale.data) - 1
  if (pcs.compute > 500) {pcs.compute <- 500}
  print(dim(so@scale.data))
  print(length(so@var.genes))
  RunPCA(object = so, pc.genes = so@var.genes, pcs.print = 1:2
    , genes.print = 5, pcs.compute = pcs.compute)
})

# ProjectPCA scores each gene in the dataset (including genes not included in
# the PCA) based on their correlation with the calculated components. Though we
# don't use this further here, it can be used to identify markers that are
# strongly correlated with cellular heterogeneity, but may not have passed
# through variable gene selection. The results of the projected PCA can be
# explored by setting use.full=T in the functions below.
# Saves in slot @pca.x.full
lso <- lapply(lso, function(so) {ProjectPCA(so)})

# # Seurat provides several useful ways of visualizing both cells and genes that
# # define the PCA, including **PrintPCA()**, **VizPCA()**, **PCAPlot()**, and
# # **PCHeatmap()** Examine  and visualize PCA results a few different ways
# pdf(paste0(outGraph, "PCAplots.pdf"), width = 6, height = 8)
# lapply(lso, function(so) {
#   # Returns the top genes ranked by the score's absolute values
#   # From @pca.obj[[1]]$rotation
#   VizPCA(so, 1:4)
# })
# dev.off()

# Plot genes with highest PC loadings
pgL <- lapply(names(lso), function(cluster) {
  so <- lso[[cluster]]
  # Genes with highest PC loadings
  ggL <- lapply(1:4, function(pc) {
    df <- rbind(data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc])[1:10])
    , data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc], decreasing = TRUE)[1:10])
      )
    df$GENE <- factor(row.names(df), levels = row.names(df))
    gg <- ggplot(df, aes(x = PC, y = GENE)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      ggtitle(paste0("Cluster: ", cluster, "\nPC ", pc))
    return(gg)
  })
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'r')
  return(pg)
})
pg <- plot_grid(plotlist = pgL, ncol = 1, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nGenes with highest PC loadings"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "PCAplots.pdf"), width = 13
  , height = 5*length(pgL), limitsize = FALSE)
################################################################################

### Determine statistically significant principal components

print("### Determine statistically significant principal components")

# PC variance plot
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  gg <- DimElbowPlot(so, reduction.type = "pca", dims.plot = 50)
  gg <- gg + ggtitle(paste0("Cluster: ", cl)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(0, max(so@dr$pca@sdev)))
  return(gg)
})
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 4)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPC variance for each cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PCElbowPlot.pdf")
  , width = 12, height = length(ggL)+4)

save(lso, file = paste0(outData, "seuratO.Robj"))
################################################################################

### Run Non-linear dimensional reduction (tSNE)

print("### Run Non-linear dimensional reduction (tSNE)")

# PC 1-20 test
lso <- lapply(lso, function(so) {
  # # Scale resolution by number of cells
  # res <- 2
  # print(str(res))
  # Set perplexity to 30 (default)
  perp <- 30
  # If < 1000 cells in cluster, set perplexity to 15
  if (ncol(so@scale.data) < 1000) {perp <- 15}
  # If < 100 cells in cluster, set perplexity to 5
  if (ncol(so@scale.data) < 100) {perp <- 5}
  so <- RunTSNE(so, dims.use = 1:20, do.fast = TRUE, perplexity = perp)
  so <- FindClusters(so, dims.use = 1:20, resolution = 0.6
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})
## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.0.6)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.0.6 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-20"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-20_Res06_tSNE.png")
  , width = 14, height = length(ggL))

## Perplexity 15 test PC 1-20
lso <- lapply(lso, function(so) {
  # Perplexity
  perp <- 15
  # If < 1000 cells in cluster, set perplexity to 15
  if (ncol(so@scale.data) < 1000) {perp <- 15}
  # If < 100 cells in cluster, set perplexity to 5
  if (ncol(so@scale.data) < 100) {perp <- 5}
  so <- RunTSNE(so, dims.use = 1:20, do.fast = TRUE, perplexity = perp)
  so <- FindClusters(so, dims.use = 1:20, resolution = 0.6
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})
## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.0.6)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.0.6 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-20 perplexity 15"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-20_Perp15_Res06_tSNE.png")
  , width = 14, height = length(ggL))

## Perplexity 60 test PC 1-20
lso <- lapply(lso, function(so) {
  # Perplexity
  perp <- 60
  # If < 1000 cells in cluster, set perplexity to 15
  if (ncol(so@scale.data) < 1000) {perp <- 15}
  # If < 100 cells in cluster, set perplexity to 5
  if (ncol(so@scale.data) < 100) {perp <- 5}
  so <- RunTSNE(so, dims.use = 1:20, do.fast = TRUE, perplexity = perp)
  so <- FindClusters(so, dims.use = 1:20, resolution = 0.6
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})
## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.0.6)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.0.6 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-20 perplexity 60"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-20_Perp60_Res06_tSNE.png")
  , width = 14, height = length(ggL))

## Clustering parameters to save PC 1-30
lso <- lapply(lso, function(so) {
  # # Scale resolution by number of cells
  # res <- 2
  # print(str(res))
  # Set perplexity to 30 (default)
  perp <- 30
  # If < 1000 cells in cluster, set perplexity to 15
  if (ncol(so@scale.data) < 1000) {perp <- 15}
  # If < 100 cells in cluster, set perplexity to 5
  if (ncol(so@scale.data) < 100) {perp <- 5}
  so <- RunTSNE(so, dims.use = 1:30, do.fast = TRUE, perplexity = perp)
  so <- FindClusters(so, dims.use = 1:30, resolution = 2
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})

## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.2)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.2 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-30 Resolution 2"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1), align = 'v', axis = 'l')
ggsave(paste0(outGraph, "PC1-40_PC1-30_Res2_tSNE.png")
  , width = 14, height = length(ggL))

## Recluster with resolution 0.6
# Clustering parameters to save
lso <- lapply(lso, function(so) {
  # # Scale resolution by number of cells
  # res <- 2
  # print(str(res))
  so <- FindClusters(so, dims.use = 1:30, resolution = 0.6
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})
## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.0.6)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.0.6 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-30 Resolution 0.6"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-30_Res06_tSNE.png")
  , width = 14, height = length(ggL))

## Recluster with resolution 0.4
# Clustering parameters to save
lso <- lapply(lso, function(so) {
  # # Scale resolution by number of cells
  # res <- 2
  # print(str(res))
  so <- FindClusters(so, dims.use = 1:30, resolution = 0.4
    , print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = so)
  return(so)
})
## Plot tSNE graphs
# Collect tSNE values
ggL <- lapply(names(lso), function(cl) {
  so <- lso[[cl]]
  # Specify cluster identity
  so@ident <- as.factor(so@meta.data$res.0.4)
  names(so@ident) <- so@meta.data$CELL
  # ggplot clustering
  ggC <- TSNEPlot(so, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE, alpha = 0.5)
  ggC <- ggC + guides(colour = guide_legend(override.aes = list(size = 3)))
  ggC <- ggC + ggtitle(paste0(
    "Cluster: ", cl
    , "\nNumber of cells: ", length(so@meta.data$res.0.4 == cl)
  ))
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-30 Resolution 0.4"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-30_Res04_tSNE.png")
  , width = 14, height = length(ggL))

# previous clustering
# Collect tSNE values
ggL <- lapply(names(lso), function(cluster) {
  so <- lso[[cluster]]
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  idx <- match(row.names(ggDF), row.names(centSO@meta.data))
  ggDF$CLUSTER_ROUND1 <- centSO@meta.data$res.0.6[idx]
  ggP <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = CLUSTER_ROUND1)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggtitle(cluster)
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nColor cluster round 2 by cluster round 1 cluster IDs"
  , "\nRemove cells < 200 genes detected in cluster"
  , "\nRemove genes detected in < 3 cells in cluster"
  , "\nOnly recluster clusters >= 100 cells"
  , "\nSeurat variable genes used for clustering"
  , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-30"
  , "\n"))
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "PC1-40_PC1-30_ClusterRound1_tSNE.png")
  , width = 14, height = length(ggL))

save(lso, file = paste0(outData, "seuratO.Robj"))
################################################################################


