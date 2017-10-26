# Damon Polioudakis
# 2017-05-28
# Clustering of Drop-seq cells by digital gene expression

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
# require(xlsx)

args <- commandArgs(trailingOnly = TRUE)
print(args)

## Input data
# Digital gene expression

# DS-002-003-004-005-006-007-008: exDF and metDF
load("../analysis/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7-8.Rdata")
# Subsetting for testing
# idx <- sample(1:ncol(exDF), 500)
# exDF <- exDF[idx]
# metDF <- metDF[idx, ]

# Seurat clustering of DS-002-003 to compare cluster identities
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

# Cell cycle markers from Macosko 2015 Table S2 to remove from variable gene
# list used for clustering
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_Cluster_DS-2-3-4-5-6-7-8.R"
outGraphPfx <- "../analysis/graphs/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMtLibBrain_"
outTablePfx <- "../analysis/tables/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMtLibBrain_"
outRdatPfx <- "../analysis/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMtLibBrain_"
# outGraphPfx <- "../analysis/graphs/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMt_"
# outTablePfx <- "../analysis/tables/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMt_"
# outRdatPfx <- "../analysis/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMt_"

## Output Directories
outDir <- dirname(outGraphPfx)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTablePfx)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdatPfx)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 16)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Basic QC and selecting cells for further analysis

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(exDF))
dense.size
sparse.size <- object.size(exDF)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where log-normalized values were passed in directly.
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
centSO <- new("seurat", raw.data = exDF)

# Number genes detected
print("Mean number of genes detected with no cells filtered:")
print(mean(colSums(centSO@raw.data >= 1)))

print("Number of cells input:")
print(ncol(exDF))
filtersDF <- data.frame(FILTER = "None", REMAINING = ncol(exDF)
  , MEAN_GENES_DETECTED = mean(colSums(centSO@raw.data >= 1)))

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
centSO <- Setup(centSO, min.cells = 3, min.genes = 200
  , do.logNormalize = TRUE, total.expr = 1e4, project = "DS-2-3-4-5-6-7-8"
  , do.scale = FALSE, do.center = FALSE)
print("Number of cells remaining after filtering cells <200 genes detected:")
# Number genes detected
print("Mean number of genes detected after filtering cells <200 genes detected:")
print(mean(colSums(centSO@data >= 1)))
print(ncol(centSO@data))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), "< 200 genes detected")
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED, mean(colSums(centSO@data >= 1))))

# While the setup function imposes a basic minimum gene-cutoff, you may want to
# filter out cells at this stage based on technical or biological parameters. 
# Seurat allows you to easily explore QC metrics and filter cells based on any 
# user-defined criteria. In the example below,  we visualize gene and molecule 
# counts, plot their relationship, and exclude cells with a clear outlier number 
# of genes detected as potential multiplets. Of course this is not a guarenteed 
# method to exclude cell doublets, but we include this as an example of 
# filtering user-defined outlier cells. We also filter cells based on the 
# percentage of mitochondrial genes present.

# nGene and nUMI are automatically calculated for every object by Seurat. For
# non-UMI data, nUMI represents the sum of the non-normalized values within a
# cell

# We calculate the percentage of mitochondrial genes here and store it in
# percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a
# common scRNA-seq QC metric.
mito.genes <- grep("^MT-", rownames(centSO@data), value = TRUE)
percent.mito <- apply(expm1(centSO@data[mito.genes, ]), 2, sum)/apply(expm1(centSO@data), 2, sum)

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
centSO <- AddMetaData(centSO, percent.mito, "percent.mito")
pdf(paste0(outGraphPfx, "QC_violinPlot.pdf"))
VlnPlot(centSO, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# GenePlot is typically used to visualize gene-gene relationships, but can be
# used for anything calculated by the object, i.e. columns in object@data.info,
# PC scores etc.

# Since there is a rare subset of cells with an outlier level of high
# mitochondrial percentage, and also low UMI content, we filter these as well
pdf(paste0(outGraphPfx, "QC_genePlot.pdf"))
par(mfrow = c(1, 2))
GenePlot(centSO, "nUMI", "percent.mito")
GenePlot(centSO, "nUMI", "nGene")
dev.off()

# We filter out cells that have unique gene counts over 2,500

# Note that accept.high and accept.low can be used to define a 'gate', and can
# filter cells not only based on nGene but on anything in the object (as in
# GenePlot above)
# centSO <- SubsetData(centSO, subset.name = "nGene", accept.high = 2500)
# Filter cells with gene counts > 3 sd above the mean
accept.high <- round(mean(centSO@data.info$nGene) + 3*sd(centSO@data.info$nGene), 0)
centSO <- SubsetData(centSO, subset.name = "nGene", accept.high = accept.high)
print(paste0("Number of cells remaining after filtering cells with >", accept.high, " (3 sd above mean) genes detected:"))
print(ncol(centSO@data))
print(paste0("Mean number of genes detected after after filtering cells with >", accept.high, " (3 sd above mean) genes detected:"))
print(mean(colSums(centSO@data >= 1)))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), paste0("> ", accept.high, " genes detected"))
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED, mean(colSums(centSO@data >= 1))))
centSO <- SubsetData(centSO, subset.name = "percent.mito", accept.high = 0.05)
print("Number of cells remaining after filtering cells with >5% of reads mapping to exons mapping to mitochondrial exons:")
print(ncol(centSO@data))
print("Number of genes detected after filtering cells with >5% of reads mapping to exons mapping to mitochondrial exons:")
print(mean(colSums(centSO@data >= 1)))
filtersDF <- data.frame(
  FILTERS = c(as.character(filtersDF$FILTER), ">5% Mt")
  , REMAINING = c(filtersDF$REMAINING, ncol(centSO@data))
  , MEAN_GENES_DETECTED = c(filtersDF$MEAN_GENES_DETECTED, mean(colSums(centSO@data >= 1))))

write.csv(filtersDF, file = paste0(outTablePfx, "Filters.csv"), quote = FALSE)
################################################################################

### Regress out unwanted sources of variation

# Your single cell dataset likely contains 'uninteresting' sources of variation.
# This could include not only technical noise, but batch effects, or even 
# biological sources of variation (cell cycle stage). As suggested in Buettner 
# et al, NBT, 2015, regressing these signals out of the analysis can improve 
# downstream dimensionality reduction and clustering. Seurat implements a basic 
# version of this by constructing linear models to predict gene expression based
# on user-defined variables. Seurat stores the z-scored residuals of these 
# models in the scale.data slot, and they are used for dimensionality reduction 
# and clustering.

# We typically regress out cell-cell variation in gene expression driven by
# batch (if applicable), cell alignment rate (as provided by Drop-seq tools for
# Drop-seq data),  the number of detected molecules, and mitochondrial gene
# expression. For cycling cells, we can also learn a 'cell-cycle' score (as in
# Macosko et al) and Regress this out as well. In this simple example here for
# post-mitotic blood cells, we simply regress on the number of detected
# molecules per cell as well as the percentage mitochondrial gene content an
# example.

# note that this overwrites centSO@scale.data. Therefore, if you intend to use
# RegressOut, you can set do.scale=F and do.center=F in the original object to
# save some time.

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

# No center or scale
noCentSO <- centSO
noCentSO <- RegressOut(centSO
  , latent.vars = c("nUMI", "percent.mito", "librarylab", "individual")
  # , latent.vars = c("nUMI", "percent.mito")
  , do.scale = FALSE, do.center = FALSE)
# Center and scale
centSO <- RegressOut(centSO
  # , latent.vars = c("nUMI", "percent.mito"))
  , latent.vars = c("nUMI", "percent.mito", "librarylab", "individual"))

### Detection of variable genes across the single cells

# Seurat calculates highly variable genes and focuses on these for downstream
# analysis. **MeanVarPlot()**, which works by calculating the average expression
# and dispersion for each gene, placing these genes into bins, and then
# calculating a z-score for dispersion within each bin. This helps control for
# the relationship between variability and average expression. This function is
# unchanged from (Macosko et al.), but new methods for variable gene expression
# identification are coming soon. We suggest that users set these parameters to
# mark visual outliers on the dispersion plot, but the exact parameter settings
# may vary based on the data type, heterogeneity in the sample, and normalization
# strategy. The parameters here identify ~2,000 variable genes, and represent
# typical parameter settings for UMI data that is normalized to a total of 1e4
# molecules.

pdf(paste0(outGraphPfx, "MeanVarPlot.pdf"))
centSO <- MeanVarPlot(centSO ,fxn.x = expMean, fxn.y = logVarDivMean
  , x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
dev.off()

print("Number of variable genes used for clustering:")
length(centSO@var.genes)

## Remove CC genes from variable genes

# Cleanup CC marker data frame
ccDF <- data.frame(lapply(ccDF, as.character), stringsAsFactors=FALSE)
cc <- c(unlist(ccDF))
cc <- gsub(" *", "", cc)
idx <- centSO@var.genes %in% cc
centSO@var.genes <- centSO@var.genes[! idx]

print("Number of variable genes used for clustering after removing CC genes:")
length(centSO@var.genes)

save(centSO, noCentSO, metDF, file = paste0(outRdatPfx, "seuratO.Robj"))
################################################################################

### Perform linear dimensional reduction

# Perform PCA on the scaled data. By default, the genes in object\@var.genes are
# used as input, but can be defined using pc.genes. We have typically found that
# running dimensionality reduction on genes with high-dispersion can improve
# performance. However, with UMI data - particularly after using RegressOut, we
# often see that PCA returns similar (albeit slower) results when run on much
# larger subsets of genes, including the whole transcriptome.
centSO <- PCA(centSO, pc.genes = centSO@var.genes, do.print = TRUE
  , pcs.print = 5, genes.print = 5)

# ProjectPCA scores each gene in the dataset (including genes not included in
# the PCA) based on their correlation with the calculated components. Though we
# don't use this further here, it can be used to identify markers that are
# strongly correlated with cellular heterogeneity, but may not have passed
# through variable gene selection. The results of the projected PCA can be
# explored by setting use.full=T in the functions below.
# Saves in slot @pca.x.full
centSO <- ProjectPCA(centSO)

# Seurat provides several useful ways of visualizing both cells and genes that
# define the PCA, including **PrintPCA()**, **VizPCA()**, **PCAPlot()**, and
# **PCHeatmap()** Examine  and visualize PCA results a few different ways
pdf(paste0(outGraphPfx, "PCAplots.pdf"))

PrintPCA(centSO, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
# Returns the top genes ranked by the score's absolute values
# From @pca.obj[[1]]$rotation
VizPCA(centSO, 1:4)
PCAPlot(centSO, 1, 2, pt.size = 0.5)

# In particular **PCHeatmap()** allows for easy exploration of the primary
# sources of heterogeneity in a dataset, and can be useful when trying to decide
# which PCs to include for further downstream analyses. Both cells and genes are
# ordered according to their PCA scores. Setting cells.use to a number plots the
# 'extreme' cells on both ends of the spectrum, which dramatically speeds
# plotting for large datasets. Though clearly a supervised analysis, we find
# this to be a valuable tool for exploring correlated gene sets.
# PCHeatmap(centSO, pc.use = 1, do.balanced = TRUE)
# PCHeatmap(centSO, pc.use = 1:6, do.balanced = TRUE, label.columns = FALSE
#   , use.full = FALSE)

dev.off()

save(centSO, noCentSO, metDF, file = paste0(outRdatPfx, "seuratO.Robj"))
################################################################################

### Determine statistically significant principal components

# To overcome the extensive technical noise in any single gene for scRNA-seq 
# data, Seurat clusters cells based on their PCA scores, with each PC
# essentially representing a 'metagene' that combines information across a
# correlated gene set. Determining how many PCs to include downstream is
# therefore an important step.

# In Macosko et al, we implemented a resampling test inspired by the jackStraw
# procedure. We randomly permute a subset of the data (1% by default) and rerun
# PCA, constructing a 'null distribution' of gene scores, and repeat this
# procedure. We identify 'significant' PCs as those who have a strong enrichment
# of low p-value genes.

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
# centSO <- JackStraw(centSO, num.replicate = 100, do.print = FALSE)

# The **JackStrawPlot()** function provides a visualization tool for comparing
# the distribution of p-values for each PC with a uniform distribution (dashed
# line). 'Significant' PCs will show a strong enrichment of genes with low
# p-values (solid curve above the dashed line). In this case it appears that PCs
# 1-10 are significant.

# JackStrawPlot(centSO, PCs = 1:12)

# A more ad hoc method for determining which PCs to use is to look at a plot of
# the standard deviations of the principle components and draw your cutoff where
# there is a clear elbow in the graph. This can be done with **PCElbowPlot()**.
# In this example, it looks like the elbow would fall around PC 9.

pdf(paste0(outGraphPfx, "PCElbowPlot.pdf"))
PCElbowPlot(centSO, num.pc = 50)
dev.off()

# # Calculate variance explained by each PC
# varExp <- (pcNmVar$sdev)^2 / sum(pcNmVar$sdev^2)
# topVar <- varExp[1:5]
# colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
#   , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

#PC selection -- identifying the true dimensionality of a dataset -- is an
#important step for Seurat, but can be challenging/uncertain for the user. We
#therefore suggest these three approaches to consider. The first is more
#supervised, exploring PCs to determine relevant sources of heterogeneity, and
#could be used in conjunction with GSEA for example. The second implements a
#statistical test based on a random null model, but is time-consuming for large
#datasets, and may not return a clear PC cutoff. The third is a heuristic that
#is commonly used, and can be calculated instantly. In this  example all three
#approaches yielded similar results, but we might have been justified in
#choosing anything between PC 7-10 as a cutoff. We followed the jackStraw  here,
#admittedly buoyed by seeing the PCHeatmap returning interpretable signals
#(including canonical dendritic cell markers) throughout these PCs. Though the
#results are only subtly affected by small shifts in this cutoff (you can test
#below), we strongly suggest always explore the PCs they choose to include
#downstream. 
################################################################################

### Cluster the cells

# Seurat now includes an graph-based clustering approach compared to (Macosko et
# al.). Importantly, the *distance metric* which drives the clustering analysis
# remains the same (based on previously identified PCs) remains the same.
# However, our approach to partioning the cellular distance matrix into clusters
# has dramatically improved. Our approach was heavily inspired by recent
# manuscripts which applied graph-based clustering approaches to scRNA-seq data
# [[SNN-Cliq, Xu and Su, Bioinformatics,
# 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract)
# and CyTOF data [[PhenoGraph, Levine et al., Cell,
# 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). Briefly, these methods
# embed cells in a graph structure - for example a K-nearest neighbor (KNN)
# graph, with edges drawn between cells with similar gene expression patterns,
# and then attempt to partition this graph into highly interconnected
# 'quasi-cliques' or 'communities'. As in PhenoGraph, we first construct a KNN
# graph based on the euclidean distance in PCA space, and refine the edge
# weights between any two cells based on the shared overlap in their local
# neighborhoods (Jaccard distance). To cluster the cells, we apply the smart
# local moving algorithm [[SLM, Blondel et al., Journal of Statistical
# Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), to
# iteratively group cell groupings together with the goal of optimizing the
# standard modularity function.

# The **FindClusters()** function implements the procedure, and contains a
# resolution parameter that sets the 'granularity' of the downstream clustering,
# with increased values leading to a greater number of clusters. We find that
# setting this parameter between 0.6-1.2 typically returns good results for
# single cell datasets of around 3K cells. Optimal resolution often increases
# for larger datasets. The clusters are saved in the object\@ident slot.

# save.SNN=T saves the SNN so that the  SLM algorithm can be rerun using the
# same graph, but with a different resolution value (see docs for full details)
# centSO <- FindClusters(centSO, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T)
################################################################################

### Run Non-linear dimensional reduction (tSNE)

# Seurat continues to use tSNE as a powerful tool to visualize and explore these
# datasets. While we no longer advise clustering directly on tSNE components,
# cells within the graph-based clusters determined above should  co-localize on
# the tSNE plot. This is because the tSNE aims to place cells with similar local
# neighborhoods in high-dimensional space together in low-dimensional space. As
# input to the tSNE, we suggest using the same PCs as input to the clustering
# analysis, although computing the tSNE based on scaled gene expression is also
# supported using the genes.use argument.

## Cluster and tSNE with different parameters

pdf(paste0(outGraphPfx, "tSNE_Res0.6_200gd_2500dg_Mt5_PCs.pdf"))

centSO <- RunTSNE(centSO, dims.use = 1:40, do.fast = TRUE)
centSO <- FindClusters(centSO, pc.use = 1:40, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
ggO <- TSNEPlot(centSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
ggO + ggtitle("PCs: 1-40")

centSO <- RunTSNE(centSO, dims.use = 1:6, do.fast = TRUE)
centSO <- FindClusters(centSO, pc.use = 1:6, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
ggO <- TSNEPlot(centSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
ggO + ggtitle("PCs: 1-6")

centSO <- RunTSNE(centSO, dims.use = 1:10, do.fast = TRUE)
centSO <- FindClusters(centSO, pc.use = 1:10, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
ggO <- TSNEPlot(centSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
ggO + ggtitle("PCs: 1-10")

centSO <- RunTSNE(centSO, dims.use = 1:15, do.fast = TRUE)
centSO <- FindClusters(centSO, pc.use = 1:10, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
ggO <- TSNEPlot(centSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
ggO + ggtitle("PCs: 1-15")

centSO <- RunTSNE(centSO, dims.use = 1:20, do.fast = TRUE)
centSO <- FindClusters(centSO, pc.use = 1:20, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
ggO <- TSNEPlot(centSO, pt.size = 0.75, do.return = TRUE, do.label = TRUE)
ggO + ggtitle("PCs: 1-20")

dev.off()

## Plot tSNE graph colored by lanes, samples, or GZ CP
# Collect tSNE values
ggDF <- centSO@tsne.rot

# metDF$LIBRARY <- "Plath"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N702"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N701"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-007" & metDF$NEXTERA == "N710"] <- "Geschwind"
# metDF$LIBRARY[metDF$SEQ_RUN == "DS-007" & metDF$NEXTERA == "N711"] <- "Geschwind"

# Add cluster identity
ggDF$CLUSTER <- centSO@ident
# Add metadata
ggDF <- data.frame(ggDF, metDF[match(row.names(ggDF), metDF$CELL), ])
ggDF$BRAIN <- as.factor(ggDF$BRAIN)

# Add DS-002-003 cluster identity
# Assign annotated cluster names to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
new.cluster.ids <- c(
  "Excitatory Upper Layer Neuron 1"
  , "Excitatory Neuron"
  , "Excitatory Upper Layer Neuron 2"
  , "Excitatory Deep Layer Neuron"
  , "Intermediate Progenitors"
  , "Interneuron"
  , "Mitotic Progenitors"
  , "oRG"
  , "Oligodendrocyte Precursor"
  , "Endothelial")
fetb@ident <- plyr::mapvalues(fetb@ident, from = current.cluster.ids
  , to = new.cluster.ids)
# Stash numerical cluster identities if want to use later
fetb <- StashIdent(fetb, save.name = "Cluster_Numbers")
preCl <- fetb@ident
names(preCl) <- gsub("_.*", "", names(preCl))
idx <- match(ggDF$CELL, names(preCl))
ggDF$PREVIOUS_CLUSTER <- preCl[idx]

# ggplot clustering
ggC <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = CLUSTER)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# ggplot clustering
ggP <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PREVIOUS_CLUSTER)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# ggsave(paste0(outGraphPfx, "tSNE_PreviousClustering_200gd_2500dg_Mt5_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Sample
ggB <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = BRAIN)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# ggsave(paste0(outGraphPfx, "tSNE_Sample_200gd_2500dg_Mt5_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Seq Run
ggS <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = SEQ_RUN)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# ggsave(paste0(outGraphPfx, "tSNE_Lane_200gd_2500dg_Mt5_PCA1to15.pdf")
#   , width = 8, height = 6)
# ggplot Region
ggR <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = REGION)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# ggplot Library
ggL <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = LIBRARY)) +
  geom_point(size = 0.3, alpha = 0.5) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-20"
    , "\n"))
# Save as separate pages
pdf(paste0(outGraphPfx, "tSNE_Metadata_200gd_2500dg_Mt5_PCA1to20.pdf")
  , width = 8, height = 6)
ggC; ggP; ggB; ggR; ggS; ggL
dev.off()
# Save as grid
pdf(paste0(outGraphPfx, "tSNE_MetadataGrid_200gd_2500dg_Mt5_PCA1to20.pdf")
  , width = 12, height = 16)
MultiPlotList(list(ggC, ggP, ggB, ggS, ggL, ggR), cols = 2)
dev.off()

# You can save the object at this point so that it can easily be loaded back in
# without having to rerun the computationally intensive steps performed above,
# or easily shared with collaborators.
save(centSO, noCentSO, metDF, file = paste0(outRdatPfx, "seuratO.Robj"))

# Heatmap of known markers genes
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
# kmDF <- kmDF[! is.na(kmDF$Gene.Symbol), ]
# # # Filter for gene symbols in data
# # kmDF <- kmDF[kmDF$Gene.Symbol %in% centSO@data@Dimnames[1][[1]], ]
# # pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"), width = 10, height = 4 + 0.2*nrow(kmDF))
# # kmDFL <- split(kmDF, kmDF$Population.Marked)
# # sapply(names(kmDFL), function(kmLName) {
# #   # print(kmLName)
# #   df <- kmDFL[[kmLName]]
# #   # print(df)
# #   print(as.character(df$Gene.Symbol))
# #   DoHeatmap(centSO, genes.use = as.character(df$Gene.Symbol)
# #     , order.by.ident = TRUE
# #     , slim.col.label = TRUE
# #     , remove.key = FALSE
# #     , main = (kmLName)
# #   )
# # })
# # dev.off()
# pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"), width = 16, height = 4 + 0.2*nrow(kmDF))
# DoHeatmap(centSO, genes.use = as.character(kmDF$Gene.Symbol)
#   , order.by.ident = TRUE
#   , slim.col.label = TRUE
#   , remove.key = FALSE
#   , main = "Known Markers")
# dev.off()
################################################################################

### Heatmap of expression of genes with highest PC loadings sorted by cluster

# Genes with highest PC loadings
ldf <- lapply(1:20, function(pc) {
  data.frame(GENES = names(sort(abs(centSO@pca.obj[[1]]$rotation[ ,pc])
    , decreasing = TRUE)[1:10]), PC = paste0("PC", pc))
})
genesDF <- do.call("rbind", ldf)

# Centered and scaled expression
# Subset to genes, merge to keep duplicate genes (from more than 1 PC)
ggDF <- merge(genesDF, centSO@scale.data, by.x = "GENES", by.y = "row.names"
  , all.x = TRUE)
# For gene ordering by cluster
# idx <- match(genesDF$GENES, ggDF$GENES)
# ggDF <- ggDF[idx, ]
row.names(ggDF) <- paste0(ggDF$GENES, " ", ggDF$PC)
ggDF <- ggDF[ ,-c(1:2)]
# Order by clustering
idx <- match(colnames(ggDF), names(sort(centSO@ident)))
ggDF <- ggDF[ ,idx]
# Format for ggplot2
ggDF <- as.matrix(ggDF)
ggDF <- melt(ggDF)
# Add clusters
idx <- match(ggDF$Var2, names(centSO@ident))
ggDF$CLUSTERS <- centSO@ident[idx]
# Add PCs
ggDF$PC <- gsub(".* ", "", ggDF$Var1)
# Set PC order
ggDF$PC <- factor(ggDF$PC, levels = unique(genesDF$PC))
# Set sample order by clustering
ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# Set expression limits
ggDF$value[ggDF$value > 2] <- 2
ggDF$value[ggDF$value < 0] <- 0
# ggplot
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_grid(PC~CLUSTERS, space = "free", scales = "free") +
  scale_fill_gradient2(name = "Normalized\nexpression"
   , high = "#d7191c", low = "white") +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
  #   , palette = 5, direction = -1) +
  # scale_fill_distiller(name = "Normalized\nexpression", palette = "Spectral") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  ylab("Genes") +
  xlab("Cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nExpression of genes with highest PC loadings"
    , "\nMean centered, variance scaled, normalized expression"
    , "\nLimits set to 0 and 2"
    , "\n"))
ggsave(paste0(outGraphPfx, "DoHeatmap_HighestLoading_CenterScale.png")
  , width = 10, height = 4 + 0.15*nrow(genesDF), dpi = 300)
# ggsave(paste0(outGraphPfx, "DoHeatmap_HighestLoading_CenterScale.pdf")
#   , width = 10, height = 4 + 0.11*nrow(genesDF))

# No centered and scaled expression
# Subset to genes, merge to keep duplicate genes (from more than 1 PC)
ggDF <- merge(genesDF, noCentSO@scale.data, by.x = "GENES", by.y = "row.names"
  , all.x = TRUE)
# For gene ordering by cluster
# idx <- match(genesDF$GENES, ggDF$GENES)
# ggDF <- ggDF[idx, ]
row.names(ggDF) <- paste0(ggDF$GENES, " ", ggDF$PC)
ggDF <- ggDF[ ,-c(1:2)]
# Order by clustering
idx <- match(colnames(ggDF), names(sort(centSO@ident)))
ggDF <- ggDF[ ,idx]
# Format for ggplot2
ggDF <- as.matrix(ggDF)
ggDF <- melt(ggDF)
# Add clusters
idx <- match(ggDF$Var2, names(centSO@ident))
ggDF$CLUSTERS <- centSO@ident[idx]
# Add PCs
ggDF$PC <- gsub(".* ", "", ggDF$Var1)
# Set PC order
ggDF$PC <- factor(ggDF$PC, levels = unique(genesDF$PC))
# Set sample order by clustering
ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# Set expression limits
ggDF$value[ggDF$value > 2] <- 2
ggDF$value[ggDF$value < 0] <- 0
# ggplot
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_grid(PC~CLUSTERS, space = "free", scales = "free") +
  scale_fill_gradient2(name = "Normalized\nexpression"
    , high = "#d7191c", low = "white") +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
  #   , palette = 5, direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  ylab("Genes") +
  xlab("Cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nExpression of genes with highest PC loadings"
    , "\nNormalized expression"
    , "\nLimits set to 0 and 2"
    , "\n"))
ggsave(paste0(outGraphPfx, "DoHeatmap_HighestLoading_NoCenterScale.png")
  , width = 10, height = 4 + 0.15*nrow(genesDF), dpi = 300)
# ggsave(paste0(outGraphPfx, "DoHeatmap_HighestLoading_NoCenterScale.pdf")
#   , width = 10, height = 4 + 0.11*nrow(genesDF))
################################################################################

### Finding differentially expressed genes (cluster biomarkers)

# Seurat can help you find markers that define clusters via differential
# expression. By default, it identifes positive and negative markers of a single
# cluster (specified in ident.1), compared to all other cells. 
# **FindAllMarkers()** automates this process for all clusters, but you can also
# test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a gene to be detected at a minimum percentage in
# either of the two groups of cells, and the thresh.test argument requires a
# gene to be differentially expressed (on average) by some amount between the
# two groups. You can set both of these to 0, but with a dramatic increase in
# time - since this will test a large number of genes that are unlikely to be
# highly discriminatory. As another option to speed up these computations,
# max.cells.per.ident can be set. This will downsample each identity class to
# have no more cells than whatever this is set to. While there is generally
# going to be a loss in power, the speed increases can be signficiant and the
# most highly differentially expressed genes will likely still rise to the top.

# # find all markers of cluster 1
# cluster1.markers <- FindMarkers(centSO, ident.1 = 1, min.pct = 0.25)
# print(head(cluster1.markers, 5))
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(centSO, 5, c(0,3), min.pct = 0.25)
# print(head(cluster1.markers, 5))
# find markers for every cluster compared to all remaining cells, report only the positive ones
centSO.markers <- FindAllMarkers(centSO, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# centSO.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
markersDF <- data.frame(centSO.markers %>% group_by(cluster) %>% top_n(20, avg_diff))

write.table(x = centSO.markers
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = markersDF
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All_Top20.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)

save(centSO, noCentSO, metDF, centSO.markers, file = paste0(outRdatPfx, "seuratO.Robj"))

# Seurat has four tests for differential expression which can be set with the
# test.use parameter: ROC test ("roc"), t-test ("t"), LRT test based on
# zero-inflated data ("bimod", default), LRT test based on tobit-censoring
# models ("tobit") The ROC test returns the 'classification power' for any
# individual marker (ranging from 0 - random, to 1 - perfect).

# cluster1.markers <- FindMarkers(centSO, ident.1 = 0, thresh.use = 0.25
# , test.use = "roc", only.pos = T)

# There are several tools for visualizing marker expression. **VlnPlot()**
# generates a violin plot which shows the probability density at different
# expression levels of the gene for each cluster. As seen below, good marker
# genes will show strongly in a single cluster. It can also be useful to look at
# gene/gene correlation with **GenePlot()** which returns a plot similar to a
# 'FACS' plot with cells colored by cluster identity. Also, the
# **FeaturePlot()** function is useful for viewing the expression of the gene in
# the context of all the cells and helps validate the specificity of the marker
# or the quality of the clustering.

# VlnPlot(centSO, c("STMN2","LIMCH1"))
# 
# #you can plot raw UMI counts as well
# VlnPlot(centSO, c("STMN2", "LIMCH1"),use.raw = T,y.log = T)
# 
# FeaturePlot(centSO, c("STMN2", "LIMCH1","SATB2","GPM6A","CSRP2","PLXNA4")
# , cols.use = c("grey","blue"))

# Heatmaps can also be a good way to examine heterogeneity within/between
# clusters. The **DoHeatmap()** function will generate a heatmap for given cells
# and genes. In this case, we are plotting the top 20 markers (or all markers if
# less than 20) for each cluster.

# Top 10 markers for each cluster
centSO.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10

# Heatmap of top 10 markers
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
# DoHeatmap(centSO, genes.use = top10$gene, order.by.ident = TRUE
#   , slim.col.label = TRUE, remove.key = FALSE)

# Centered and scaled
# Subset to genes, merge to keep duplicate genes (from more than 1 PC)
ggDF <- merge(top10, centSO@scale.data, by.x = "gene", by.y = "row.names"
  , all.x = TRUE)
# For gene ordering by cluster
idx <- match(rev(top10$gene), ggDF$gene)
ggDF <- ggDF[idx, ]
row.names(ggDF) <- paste0(length(ggDF$gene):1, "_", ggDF$gene)
ggDF <- ggDF[ ,-c(1:6)]
# Order by clustering
idx <- match(colnames(ggDF), names(sort(centSO@ident)))
ggDF <- ggDF[ ,idx]
# Format for ggplot2
ggDF <- as.matrix(ggDF)
ggDF <- melt(ggDF)
# Add clusters
idx <- match(ggDF$Var2, names(centSO@ident))
ggDF$CLUSTERS <- centSO@ident[idx]
# Set sample order by clustering
ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# Set expression limits
ggDF$value[ggDF$value > 2] <- 2
ggDF$value[ggDF$value < 0] <- 0
# ggplot
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_grid(~CLUSTERS, space = "free", scales = "free") +
  scale_fill_gradient2(name = "Normalized\nexpression"
    , high = "#d7191c", low = "white") +
  scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
  # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
  #   , palette = 5, direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  ylab("Genes") +
  xlab("Cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nTop 10 marker genes per cluster"
    , "\nMean centered, variance scaled, normalized expression"
    , "\n"))
ggsave(paste0(outGraphPfx, "DoHeatmap_Top10Markers_CenterScale.png")
  , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)

# No centered and scaled expression
# Subset to genes, merge to keep duplicate genes (from more than 1 PC)
ggDF <- merge(top10, noCentSO@scale.data, by.x = "gene", by.y = "row.names"
  , all.x = TRUE)
# For gene ordering by cluster
idx <- match(rev(top10$gene), ggDF$gene)
ggDF <- ggDF[idx, ]
row.names(ggDF) <- paste0(length(ggDF$gene):1, "_", ggDF$gene)
ggDF <- ggDF[ ,-c(1:6)]
# Order by clustering
idx <- match(colnames(ggDF), names(sort(centSO@ident)))
ggDF <- ggDF[ ,idx]
# Format for ggplot2
ggDF <- as.matrix(ggDF)
ggDF <- melt(ggDF)
# Add clusters
idx <- match(ggDF$Var2, names(centSO@ident))
ggDF$CLUSTERS <- centSO@ident[idx]
# Set sample order by clustering
ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
# Set expression limits
# Set expression limits
ggDF$value[ggDF$value > 2] <- 2
ggDF$value[ggDF$value < 0] <- 0
# ggplot
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_grid(~CLUSTERS, space = "free", scales = "free") +
  scale_fill_gradient2(name = "Normalized\nexpression"
    , high = "#d7191c", low = "white") +
  scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
  # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
  #   , palette = 5, direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  ylab("Genes") +
  xlab("Cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nTop 10 marker genes"
    , "\nNormalized expression"
    , "\n"))
ggsave(paste0(outGraphPfx, "DoHeatmap_Top10Markers.png")
  , width = 10, height = 4 + 0.15*nrow(top10), dpi = 300)

# Violin plots of top 10 markers
pdf(paste0(outGraphPfx, "ViolinPlot_Top10Markers.pdf"), width = 10)
top10L <- split(top10, top10$cluster)
lapply(top10L, function(top10cluster) {
  VlnPlot(centSO, top10cluster$gene, size.use = 0.5)
})
dev.off()

# Feature plot of top 10 markers
pdf(paste0(outGraphPfx, "FeaturePlot_Top10Markers.pdf"), width = 10)
top10L <- split(top10, top10$cluster)
lapply(top10L, function(top10cluster) {
  FeaturePlot(centSO, top10cluster$gene, cols.use = c("grey","blue")
    , pt.size = 0.7)
})
dev.off()
################################################################################
# 
# ###Assigning cell type identity to clusters
# # Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:
# 
# #   Cluster ID | Markers       | Cell Type
# # -----------|---------------|----------
# #   0          | IL7R          | CD4 T cells
# # 1          | CD14, LYZ     | CD14+ Monocytes
# # 2          | MS4A1         | B cells
# # 3          | CD8A          | CD8 T cells
# # 4          | FCGR3A, MS4A7 | FCGR3A+ Monocytes
# # 5          | GNLY, NKG7    | NK cells
# # 6          | FCER1A, CST3  | Dendritic Cells
# # 7          | PPBP          | Megakaryocytes
# 
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
# new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
# centSO@ident <- plyr::mapvalues(centSO@ident, from = current.cluster.ids, to = new.cluster.ids)
# TSNEPlot(centSO, do.label = T, pt.size = 0.5)
# 
# ###Further subdivisions within cell types If you perturb some of our parameter
# #choices above (for example, setting  `resolution=0.8` or changing the number of
# #PCs), you might see the CD4 T cells subdivide into two groups. You can explore
# #this subdivision to find markers separating the two T cell subsets. However,
# #before reclustering (which will overwrite `object@ident`), we can stash our
# #renamed identities to be easily recovered later.
# 
# #First lets stash our identities for later
# centSO <- StashIdent(centSO, save.name = "ClusterNames_0.6")
# 
# # Note that if you set save.snn=T above, you don't need to recalculate the SNN, and can simply put : centSO=FindClusters(centSO,resolution = 0.8)
# centSO <- FindClusters(centSO, pc.use = 1:10, resolution = 0.8, print.output = F)
# 
# # Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
# plot1 <- TSNEPlot(centSO, do.return = T, no.legend = TRUE, do.label = T)
# plot2 <- TSNEPlot(centSO, do.return = T, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = T)
# MultiPlotList(list(plot1, plot2), cols = 2)
# 
# #Find discriminating markers
# tcell.markers <- FindMarkers(centSO, 0, 1)
# 
# # Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can
# # see that CCR7 is upregulated in C0, strongly indicating that we can
# # differentiate memory from naive CD4 cells.
# 
# # cols.use demarcates the color palette from low to high expression
# FeaturePlot(centSO, c("LIMCH1", "MEF2C"), cols.use = c("green", "blue"))
# 
# # The memory/naive split is bit weak, and we would probably benefit from looking
# # at more cells to see if this becomes more convincing (stay tuned!). In the
# # meantime, we can restore our old cluster identities for downstream processing.
# 
# centSO <- SetAllIdent(centSO, id = "ClusterNames_0.6")





