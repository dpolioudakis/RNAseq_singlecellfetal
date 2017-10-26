# Damon Polioudakis
# 2016-12-06
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
# require(xlsx)

args <- commandArgs(trailingOnly = TRUE)
print(args)

## Input data
# Digital gene expression

# DS-002
cs1ExDF <- read.table(args[1], header = TRUE)
vs1ExDF <- read.table(args[2], header = TRUE)
vh1ExDF <- read.table(args[3], header = TRUE)
ch1ExDF <- read.table(args[4], header = TRUE)
# DS-003
cs2ExDF <- read.table(args[5], header = TRUE)
vs2ExDF <- read.table(args[6], header = TRUE)
vh2ExDF <- read.table(args[7], header = TRUE)
ch2ExDF <- read.table(args[8], header = TRUE)

# # DS-002
# cs1ExDF <- read.table("../DS-002/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsXbp083L1/N701/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# vs1ExDF <- read.table("../DS-002/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsXbp083L1/N702/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# vh1ExDF <- read.table("../DS-002/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsXbp083L1/N703/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# ch1ExDF <- read.table("../DS-002/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsXbp083L1/N704/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# # DS-003
# cs2ExDF <- read.table("../DS-003/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsVAP048L8/N706/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# vs2ExDF <- read.table("../DS-003/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsVAP048L8/N707/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# vh2ExDF <- read.table("../DS-003/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsVAP048L8/N708/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# ch2ExDF <- read.table("../DS-003/data/digital_gene_expression/GRCh37_75_assembly_NoERCC/SxaQSEQsVAP048L8/N709/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)

# Known cell type markers from Luis
# kmDF <- read.xlsx("../source/MarkersforSingleCell.xlsx", 1, header = TRUE)

## Variables
outGraphPfx <- args[9]
outTablePfx <- args[10]
outRdatPfx <- args[11]
graphCodeTitle <- "Cluster_Seurat.R"
# outGraphPfx <- "../analysis/graphs/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_"
# outTablePfx <- "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_"
# outRdatPfx <- "../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_"

## Output Directories
outDir <- dirname(args[9])
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(args[10])
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdatPfx)
dir.create(outRdatDir, recursive = TRUE)
# outDir <- dirname(outGraphPfx)
# outTableDir <- dirname(outTablePfx)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Compile sample data frames for Seurat

## Add sample name to cell ID
names(vh1ExDF)[-1] <- paste0(names(vh1ExDF)[-1], "_VZH1")
names(ch1ExDF)[-1] <- paste0(names(ch1ExDF)[-1], "_CPH1")
names(vs1ExDF)[-1] <- paste0(names(vs1ExDF)[-1], "_VZS1")
names(cs1ExDF)[-1] <- paste0(names(cs1ExDF)[-1], "_CPS1")
names(vh2ExDF)[-1] <- paste0(names(vh2ExDF)[-1], "_VZH2")
names(ch2ExDF)[-1] <- paste0(names(ch2ExDF)[-1], "_CPH2")
names(vs2ExDF)[-1] <- paste0(names(vs2ExDF)[-1], "_VZS2")
names(cs2ExDF)[-1] <- paste0(names(cs2ExDF)[-1], "_CPS2")

## Combine samples into 1 dataframe
exLDF <- list(cs1ExDF, vs1ExDF, vh1ExDF, ch1ExDF, cs2ExDF, vs2ExDF, vh2ExDF, ch2ExDF)
# exLDF <- list(chExDF, vhExDF)
for (i in 1:length(exLDF)) {
  df <- exLDF[[i]]
  if (i == 1) {
    exDF <- df}
  else {
    exDF <- merge(exDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
  }
}
str(exDF)
exDF[is.na(exDF)] <- 0
row.names(exDF) <- exDF$GENE
exDF <- exDF[ ,-1]
print("Number of cells input:")
print(ncol(exDF))
################################################################################

fetb.data <- exDF

#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(fetb.data))
dense.size
sparse.size <- object.size(fetb.data)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where
# log-normalized values were passed in directly. You can continue to pass in
# log-normalized values, just set do.logNormalize=F in the next step.
fetb <- new("seurat", raw.data = fetb.data)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes 
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules
# (as in Macosko et al. Cell 2015)
fetb <- Setup(fetb, min.cells = 3, min.genes = 200, do.logNormalize = T
  , total.expr = 1e4, project = "10X_fetb")
print("Number of cells remaining after filtering cells <200 genes detected:")
print(ncol(fetb@data))
################################################################################

### Basic QC and selecting cells for further analysis
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
mito.genes <- grep("^MT-", rownames(fetb@data), value = T)
percent.mito <- apply(expm1(fetb@data[mito.genes, ]), 2, sum)/apply(expm1(fetb@data), 2, sum)

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
fetb <- AddMetaData(fetb, percent.mito, "percent.mito")
pdf(paste0(outGraphPfx, "QC_violinPlot.pdf"))
VlnPlot(fetb, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# GenePlot is typically used to visualize gene-gene relationships, but can be
# used for anything calculated by the object, i.e. columns in object@data.info,
# PC scores etc.

# Since there is a rare subset of cells with an outlier level of high
# mitochondrial percentage, and also low UMI content, we filter these as well
pdf(paste0(outGraphPfx, "QC_genePlot.pdf"))
par(mfrow = c(1, 2))
GenePlot(fetb, "nUMI", "percent.mito")
GenePlot(fetb, "nUMI", "nGene")
dev.off()

# We filter out cells that have unique gene counts over 2,500

# Note that accept.high and accept.low can be used to define a 'gate', and can
# filter cells not only based on nGene but on anything in the object (as in
# GenePlot above)
fetb <- SubsetData(fetb, subset.name = "nGene", accept.high = 2500)
print("Number of cells remaining after filtering cells with >2500 genes detected:")
print(ncol(fetb@data))
fetb <- SubsetData(fetb, subset.name = "percent.mito", accept.high = 0.05)
print("Number of cells remaining after filtering cells with >5% of reads mapping to exons mapping to mitochondrial exons:")
print(ncol(fetb@data))
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

# note that this overwrites fetb@scale.data. Therefore, if you intend to use
# RegressOut, you can set do.scale=F and do.center=F in the original object to
# save some time.
fetb <- RegressOut(fetb, latent.vars = c("nUMI", "percent.mito"))

### Detection of variable genes across the single cells Seurat calculates highly
#variable genes and focuses on these for downstream analysis. **MeanVarPlot()**,
#which works by calculating the average expression and dispersion for each gene,
#placing these genes into bins, and then calculating a z-score for dispersion
#within each bin. This helps control for the relationship between variability
#and average expression. This function is unchanged from (Macosko et al.), but
#new methods for variable gene expression identification are coming soon. We
#suggest that users set these parameters to mark visual outliers on the
#dispersion plot, but the exact parameter settings may vary based on the data
#type, heterogeneity in the sample, and normalization strategy. The parameters
#here identify ~2,000 variable genes, and represent typical parameter settings
#for UMI data that is normalized to a total of 1e4 molecules.

pdf(paste0(outGraphPfx, "MeanVarPlot.pdf"))
fetb <- MeanVarPlot(fetb ,fxn.x = expMean, fxn.y = logVarDivMean
  , x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
dev.off()

length(fetb@var.genes)

# Write variable genes used by Seurat for clustering to table
write.table(x = fetbO@var.genes
  , file = paste0(outTablePfx, "VariableGenes.txt")
  , quote = FALSE, col.names = FALSE, row.names = FALSE)
################################################################################

### Perform linear dimensional reduction

# Perform PCA on the scaled data. By default, the genes in object\@var.genes are
# used as input, but can be defined using pc.genes. We have typically found that
# running dimensionality reduction on genes with high-dispersion can improve
# performance. However, with UMI data - particularly after using RegressOut, we
# often see that PCA returns similar (albeit slower) results when run on much
# larger subsets of genes, including the whole transcriptome.
fetb <- PCA(fetb, pc.genes = fetb@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

# ProjectPCA scores each gene in the dataset (including genes not included in
# the PCA) based on their correlation with the calculated components. Though we
# don't use this further here, it can be used to identify markers that are
# strongly correlated with cellular heterogeneity, but may not have passed
# through variable gene selection. The results of the projected PCA can be
# explored by setting use.full=T in the functions below.
fetb <- ProjectPCA(fetb)

# Seurat provides several useful ways of visualizing both cells and genes that
# define the PCA, including **PrintPCA()**, **VizPCA()**, **PCAPlot()**, and
# **PCHeatmap()** Examine  and visualize PCA results a few different ways
pdf(paste0(outGraphPfx, "PCAplots.pdf"))
PrintPCA(fetb, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
VizPCA(fetb, 1:2)
PCAPlot(fetb, 1, 2)

# In particular **PCHeatmap()** allows for easy exploration of the primary
# sources of heterogeneity in a dataset, and can be useful when trying to decide
# which PCs to include for further downstream analyses. Both cells and genes are
# ordered according to their PCA scores. Setting cells.use to a number plots the
# 'extreme' cells on both ends of the spectrum, which dramatically speeds
# plotting for large datasets. Though clearly a supervised analysis, we find
# this to be a valuable tool for exploring correlated gene sets.
# PCHeatmap(fetb, pc.use = 1, do.balanced = TRUE)
PCHeatmap(fetb, pc.use = 1:6, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()
################################################################################

###Determine statistically significant principal components To overcome the
#extensive technical noise in any single gene for scRNA-seq data, Seurat
#clusters cells based on their PCA scores, with each PC essentially representing
#a 'metagene' that combines information across a correlated gene set.
#Determining how many PCs to include downstream is therefore an important step.

# In Macosko et al, we implemented a resampling test inspired by the jackStraw
# procedure. We randomly permute a subset of the data (1% by default) and rerun
# PCA, constructing a 'null distribution' of gene scores, and repeat this
# procedure. We identify 'significant' PCs as those who have a strong enrichment
# of low p-value genes.

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
# fetb <- JackStraw(fetb, num.replicate = 100, do.print = FALSE)

# The **JackStrawPlot()** function provides a visualization tool for comparing
# the distribution of p-values for each PC with a uniform distribution (dashed
# line). 'Significant' PCs will show a strong enrichment of genes with low
# p-values (solid curve above the dashed line). In this case it appears that PCs
# 1-10 are significant.

# JackStrawPlot(fetb, PCs = 1:12)

# A more ad hoc method for determining which PCs to use is to look at a plot of
# the standard deviations of the principle components and draw your cutoff where
# there is a clear elbow in the graph. This can be done with **PCElbowPlot()**.
# In this example, it looks like the elbow would fall around PC 9.

pdf(paste0(outGraphPfx, "PCElbowPlot.pdf"))
PCElbowPlot(fetb)
dev.off()

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

###Cluster the cells

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
fetb <- FindClusters(fetb, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T)
################################################################################

###Run Non-linear dimensional reduction (tSNE)

# Seurat continues to use tSNE as a powerful tool to visualize and explore these
# datasets. While we no longer advise clustering directly on tSNE components,
# cells within the graph-based clusters determined above should  co-localize on
# the tSNE plot. This is because the tSNE aims to place cells with similar local
# neighborhoods in high-dimensional space together in low-dimensional space. As
# input to the tSNE, we suggest using the same PCs as input to the clustering
# analysis, although computing the tSNE based on scaled gene expression is also
# supported using the genes.use argument.

fetb <- FindClusters(fetb, pc.use = 1:40, resolution = 0.6, print.output = 0, save.SNN = T)
fetb <- RunTSNE(fetb, dims.use = 1:40, do.fast = T)
pdf(paste0(outGraphPfx, "tSNE_Res0.6_200gd_2500dg_Mt5_PCA1to40.pdf"))
TSNEPlot(fetb, pt.size = 0.75)
dev.off()

fetb <- FindClusters(fetb, pc.use = 1:6, resolution = 0.6, print.output = 0, save.SNN = T)
fetb <- RunTSNE(fetb, dims.use = 1:6, do.fast = T)
pdf(paste0(outGraphPfx, "tSNE_Res0.6_200gd_2500dg_Mt5_PCA1to6.pdf"))
TSNEPlot(fetb, pt.size = 0.75)
dev.off()

fetb <- FindClusters(fetb, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T)
fetb <- RunTSNE(fetb, dims.use = 1:10, do.fast = T)
# note that you can set do.label=T to help label individual clusters
pdf(paste0(outGraphPfx, "tSNE_Res0.6_200gd_2500dg_Mt5_PCA1to10.pdf"))
gg <- TSNEPlot(fetb, pt.size = 0.75, do.return = TRUE)
gg + scale_colour_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c'
  ,'#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))
gg
dev.off()

## Plot tSNE graph colored by cell batch
# Collect tSNE values
ggDF <- fetb@tsne.rot
# Format sample names for ggplot2
ggDF$SAMPLE <- gsub(".*_", "", rownames(ggDF))
ggDF$SAMPLE[ggDF$SAMPLE == "CPS1" | ggDF$SAMPLE == "CPH1" | ggDF$SAMPLE == "VZS1" | ggDF$SAMPLE == "VZH1"] <- "Batch 1"
ggDF$SAMPLE[ggDF$SAMPLE == "CPS2" | ggDF$SAMPLE == "CPH2" | ggDF$SAMPLE == "VZS2" | ggDF$SAMPLE == "VZH2"] <- "Batch 2"

ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = SAMPLE)) +
  geom_point(size = 0.75) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-10"
    , "\n"))
ggsave(paste0(outGraphPfx, "tSNE_CellBatch_200gd_2500dg_Mt5_PCA1to10.pdf"), width = 8, height = 6)

## Plot tSNE graph colored by VZ or CP samples hard and soft dissociation
# Collect tSNE values
ggDF <- fetb@tsne.rot
# Format sample names for ggplot2
ggDF$SAMPLE <- gsub(".*_", "", rownames(ggDF))
ggDF$SAMPLE[ggDF$SAMPLE == "CPH1" | ggDF$SAMPLE == "CPH2"] <- "Cortical plate - harsh dissociation"
ggDF$SAMPLE[ggDF$SAMPLE == "CPS1" | ggDF$SAMPLE == "CPS2"] <- "Cortical plate - soft dissociation"
ggDF$SAMPLE[ggDF$SAMPLE == "VZH1" | ggDF$SAMPLE == "VZH2"] <- "Ventricular zone - harsh dissociation"
ggDF$SAMPLE[ggDF$SAMPLE == "VZS1" | ggDF$SAMPLE == "VZS2"] <- "Ventricular zone - soft dissociation"

ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = SAMPLE)) +
  geom_point(size = 0.75) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-10"
    , "\n"))
ggsave(paste0(outGraphPfx, "tSNE_VZCPHardSoft_200gd_2500dg_Mt5_PCA1to10.pdf"), width = 8, height = 6)

## Plot tSNE graph colored by VZ or CP samples
# Collect tSNE values
ggDF <- fetb@tsne.rot
# Format sample names for ggplot2
ggDF$SAMPLE <- gsub(".*_", "", rownames(ggDF))
ggDF$SAMPLE[ggDF$SAMPLE == "CPH1" | ggDF$SAMPLE == "CPH2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "CPS1" | ggDF$SAMPLE == "CPS2"] <- "Cortical plate"
ggDF$SAMPLE[ggDF$SAMPLE == "VZH1" | ggDF$SAMPLE == "VZH2"] <- "Germinal zones"
ggDF$SAMPLE[ggDF$SAMPLE == "VZS1" | ggDF$SAMPLE == "VZS2"] <- "Germinal zones"

ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = SAMPLE)) +
  geom_point(size = 0.75) +
  scale_color_manual(values = c("#fb8072", "#80b1d3")) +
  ggtitle(paste0(graphCodeTitle
    , "\nSubset cells >200 genes detected"
    , "\nRemove cells with unique gene counts > 2,500"
    , "\nRemove cell with MT transcripts > 5%"
    , "\ntSNE PCA 1-10"
    , "\n"))
ggsave(paste0(outGraphPfx, "tSNE_VZCP_200gd_2500dg_Mt5_PCA1to10.pdf"), width = 8, height = 6)

# You can save the object at this point so that it can easily be loaded back in
# without having to rerun the computationally intensive steps performed above,
# or easily shared with collaborators.
print("Saving Seurat.Robj too:")
print(paste0(outRdatPfx, "fetb_seurat.Robj"))
save(fetb, file = paste0(outRdatPfx, "fetb_seurat.Robj"))

# Heatmap of known markers genes
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
# kmDF <- kmDF[! is.na(kmDF$Gene.Symbol), ]
# # Filter for gene symbols in data
# kmDF <- kmDF[kmDF$Gene.Symbol %in% fetb@data@Dimnames[1][[1]], ]
# pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"), width = 10, height = 4 + 0.2*nrow(kmDF))
# kmDFL <- split(kmDF, kmDF$Population.Marked)
# sapply(names(kmDFL), function(kmLName) {
#   # print(kmLName)
#   df <- kmDFL[[kmLName]]
#   # print(df)
#   print(as.character(df$Gene.Symbol))
#   DoHeatmap(fetb, genes.use = as.character(df$Gene.Symbol)
#     , order.by.ident = TRUE
#     , slim.col.label = TRUE
#     , remove.key = FALSE
#     , main = (kmLName)
#   )
# })
# dev.off()
# pdf(paste0(outGraphPfx, "DoHeatmap_KnownMarkers.pdf"), width = 8, height = 4 + 0.2*nrow(kmDF))
# DoHeatmap(fetb, genes.use = as.character(kmDF$Gene.Symbol)
#   , order.by.ident = TRUE
#   , slim.col.label = TRUE
#   # Turn off mean centering
#   , do.scale = FALSE
#   # Turn off histogram trace on key
#   # , density.info = "none"
#   , remove.key = FALSE
#   , col.use = my_palette
#   , main = "Known Markers")
# dev.off()
################################################################################

###Finding differentially expressed genes (cluster biomarkers)

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
# cluster1.markers <- FindMarkers(fetb, ident.1 = 1, min.pct = 0.25)
# print(head(cluster1.markers, 5))
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(fetb, 5, c(0,3), min.pct = 0.25)
# print(head(cluster1.markers, 5))
# find markers for every cluster compared to all remaining cells, report only the positive ones
fetb.markers <- FindAllMarkers(fetb, only.pos = TRUE, min.pct = 0.25
  , thresh.use = 0.25, return.thresh = 1)
# fetb.markers <- FindAllMarkers(fetb, only.pos = TRUE, min.pct = 0
#   , thresh.use = 0, return.thresh = 1)
# fetb.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
markersDF <- data.frame(fetb.markers %>% group_by(cluster) %>% top_n(20, avg_diff))

write.table(x = fetb.markers
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = markersDF
  , file = paste0(outTablePfx, "Marker_Genes_Clusters_Vs_All_Top20.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)


# Seurat has four tests for differential expression which can be set with the
# test.use parameter: ROC test ("roc"), t-test ("t"), LRT test based on
# zero-inflated data ("bimod", default), LRT test based on tobit-censoring
# models ("tobit") The ROC test returns the 'classification power' for any
# individual marker (ranging from 0 - random, to 1 - perfect).

# cluster1.markers <- FindMarkers(fetb, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = T)

# There are several tools for visualizing marker expression. **VlnPlot()**
# generates a violin plot which shows the probability density at different
# expression levels of the gene for each cluster. As seen below, good marker
# genes will show strongly in a single cluster. It can also be useful to look at
# gene/gene correlation with **GenePlot()** which returns a plot similar to a
# 'FACS' plot with cells colored by cluster identity. Also, the
# **FeaturePlot()** function is useful for viewing the expression of the gene in
# the context of all the cells and helps validate the specificity of the marker
# or the quality of the clustering.

# VlnPlot(fetb, c("STMN2","LIMCH1"))
# 
# #you can plot raw UMI counts as well
# VlnPlot(fetb, c("STMN2", "LIMCH1"),use.raw = T,y.log = T)
# 
# FeaturePlot(fetb, c("STMN2", "LIMCH1","SATB2","GPM6A","CSRP2","PLXNA4"), cols.use = c("grey","blue"))

# Heatmaps can also be a good way to examine heterogeneity within/between
# clusters. The **DoHeatmap()** function will generate a heatmap for given cells
# and genes. In this case, we are plotting the top 20 markers (or all markers if
# less than 20) for each cluster.

# Top 10 markers for each cluster
fetb.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10

# Heatmap of markers
# Top 10
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
pdf(paste0(outGraphPfx, "DoHeatmap_Top10Markers.pdf")
  , width = 8, height = 4 + 0.1*nrow(top10))
my_palette <- colorRampPalette(c("black", "yellow"))(n = 100)
DoHeatmap(fetb, genes.use = top10$gene
  , order.by.ident = TRUE
  , slim.col.label = TRUE
  # Turn off mean centering
  , do.scale = FALSE
  # Turn off histogram trace on key
  # , density.info = "none"
  , remove.key = FALSE
  , col.use = my_palette
  , main = "Top 10 Markers")
dev.off()

# Violin plots of top 10 markers
pdf(paste0(outGraphPfx, "ViolinPlot_Top10Markers.pdf"), width = 10)
top10L <- split(top10, top10$cluster)
lapply(top10L, function(top10cluster) {
  VlnPlot(fetb, top10cluster$gene, size.use = 0.5)
})
dev.off()

# Feature plot of top 10 markers
pdf(paste0(outGraphPfx, "FeaturePlot_Top10Markers.pdf"), width = 10)
top10L <- split(top10, top10$cluster)
lapply(top10L, function(top10cluster) {
  FeaturePlot(fetb, top10cluster$gene, cols.use = c("grey","blue"), pt.size = 0.7)
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
# fetb@ident <- plyr::mapvalues(fetb@ident, from = current.cluster.ids, to = new.cluster.ids)
# TSNEPlot(fetb, do.label = T, pt.size = 0.5)
# 
# ###Further subdivisions within cell types If you perturb some of our parameter
# #choices above (for example, setting  `resolution=0.8` or changing the number of
# #PCs), you might see the CD4 T cells subdivide into two groups. You can explore
# #this subdivision to find markers separating the two T cell subsets. However,
# #before reclustering (which will overwrite `object@ident`), we can stash our
# #renamed identities to be easily recovered later.
# 
# #First lets stash our identities for later
# fetb <- StashIdent(fetb, save.name = "ClusterNames_0.6")
# 
# # Note that if you set save.snn=T above, you don't need to recalculate the SNN, and can simply put : fetb=FindClusters(fetb,resolution = 0.8)
# fetb <- FindClusters(fetb, pc.use = 1:10, resolution = 0.8, print.output = F)
# 
# # Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
# plot1 <- TSNEPlot(fetb, do.return = T, no.legend = TRUE, do.label = T)
# plot2 <- TSNEPlot(fetb, do.return = T, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = T)
# MultiPlotList(list(plot1, plot2), cols = 2)
# 
# #Find discriminating markers
# tcell.markers <- FindMarkers(fetb, 0, 1)
# 
# # Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can
# # see that CCR7 is upregulated in C0, strongly indicating that we can
# # differentiate memory from naive CD4 cells.
# 
# # cols.use demarcates the color palette from low to high expression
# FeaturePlot(fetb, c("LIMCH1", "MEF2C"), cols.use = c("green", "blue"))
# 
# # The memory/naive split is bit weak, and we would probably benefit from looking
# # at more cells to see if this becomes more convincing (stay tuned!). In the
# # meantime, we can restore our old cluster identities for downstream processing.
# 
# fetb <- SetAllIdent(fetb, id = "ClusterNames_0.6")





