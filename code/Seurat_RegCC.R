

# Required hoffman2 modules: module load gcc/4.9.3; module load R/3.3.0

rm(list=ls())

require(Seurat)

# Log normalized, regressed nUMI and percent mito
load("../../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# Variables
outGraphPfx <- "Seurat_RegCC_"
graphCodeTitle <- "Seurat_RegCC.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 10))
################################################################################

### Functions

# Run PCA on expression data frame
PCA_Prcomp <- function(exDF) {
  meanCenteredM <- t(scale(t(exDF), scale = TRUE))
  # Remove NA rows
  print("Dim of mean centered matrix:")
  print(dim(meanCenteredM))
  meanCenteredM <- na.omit(meanCenteredM)
  print("Dim of mean centered matrix after removing NAs:")
  print(dim(meanCenteredM))
  ## Run PCA
  PCA <- prcomp(t(meanCenteredM), center = FALSE, scale = FALSE)
  return(PCA)
}

# Calculate mean expression of group of genes for each cell using seurat scaled
# expression values
Mean_Expression <- function(genes) {
  genesExDF <- RegCcNoScCtSO@scale.data[
    which(row.names(RegCcNoScCtSO@scale.data) %in% genes), ]
  # genesExDF <- as.matrix(regCcSO@data)[which(row.names(as.matrix(regCcSO@data)) %in% genes), ]
  # Only calculate column means if there are multiple genes
  print("Length genes:")
  print(length(genes))
  if (is.matrix(genesExDF)) {
    mnExDF <- colMeans(genesExDF)
  } else {
    mnExDF <- genesExDF
  }
  # Add to ggplot data frame
  ggDF$EXPRESSION <- mnExDF[match(row.names(ggDF), names(mnExDF))]
  return(ggDF)
}

# Transform data to desired limits for ggplot2
Set_Limits <- function(ggDF) {
  ggDF$EXPRESSION[ggDF$EXPRESSION < 0] <- 0
  ggDF$EXPRESSION[ggDF$EXPRESSION > 2] <- 2
  return(ggDF)  
}

# Color tSNE plot by expression from Mean_Expression()
Feature_Plot <- function(ggDF, grouping, genes) {
  print("Feature_Plot...")
  ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = 0.5) +
    scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
      , high = "#e31a1c", limits = c(0, 2)) +
    ggtitle(paste0(graphCodeTitle
      , "\n\n", grouping
      , "\n(", paste(as.character(genes), collapse = ", "), ")"
      , "\ntSNE plot, each point is a cell"
      , "\nColor indicates genes of interest mean expression"
      , "\n"))
  return(ggFp)
}

# tSNE plot colored by Seurat clustering
TSNE_Plot <- function(regCcSO, grouping) {
  # tSNE graph colored by cluster
  ggTsne <- TSNEPlot(regCcSO, do.label = TRUE, pt.size = 0.5, do.return = TRUE)
  ggTsne <- ggTsne + ggtitle(paste0(graphCodeTitle
    , "\n\n", grouping
    , "\ntSNE plot, each point is a cell"
    , "\nColor indicates cluster assignment"
    , "\nClusters annotated manually by expression of known marker genes"
    , "\n"))
  return(ggTsne)
}
################################################################################

### Format data

# Log normalized, regressed nUMI and percent mito expression matrix
exM <- seuratO@scale.data

# Cleanup CC marker data frame
ccDF <- data.frame(lapply(ccDF, as.character), stringsAsFactors=FALSE)
cc <- c(unlist(ccDF))
cc <- gsub(" *", "", cc)

# Subset to CC genes
# Log and read depth normalized
ccExM <- seuratO@scale.data[row.names(seuratO@scale.data) %in% cc, ]
################################################################################

### PCA

## PCA CC genes

# Log2 and read depth normalized
pcaCC <- PCA_Prcomp(ccExM)
str(pcaCC)
# topPCs <- ccPCA$rotation[ ,1:10];
# # Calculate variance explained by each PC
# varExp <- (pcaCC$sdev)^2 / sum(pcaCC$sdev^2)
# Plot variance of each PC
pdf(paste0(outGraphPfx, "CCpcaVar.pdf"))
plot(pcaCC, type = "line")
dev.off()
summary(pcaCC)$importance[ ,1:5]

## PCA all genes

# Log2 and read depth normalized
pcaAll <- PCA_Prcomp(exM)
pdf(paste0(outGraphPfx, "pcaVar.pdf"))
plot(pcaAll, type = "line")
dev.off()
summary(pcaAll)$importance[ ,1:5]

save(pcaAll, pcaCC, file = "Seurat_RegCC_PCA.rdat")

# Correlation PCs all genes vs CC genes
corDF <- data.frame(cor(pcaAll$x[ ,1:5], pcaCC$x[ ,1:5]))
rownames(corDF) <- paste0("AllGenes_", rownames(corDF))
colnames(corDF) <- paste0("CCgenes_", colnames(corDF))
corDF <- round(corDF, 2)
write.csv(corDF, file = paste0(outGraphPfx, "PCA_Correlation.csv")
  , quote = FALSE)
################################################################################

### Regress out cell cycle genes PC1

pc1CC <- pcaCC$x[ ,1]

# Add PC1 of CC genes to seurat object
seuratO <- AddMetaData(seuratO, pc1CC, "PC1CellCycleGenes")

# Regress out nUMI, percet Mt, and PC1 CC genes, scales and centers rows
regCcSO <- RegressOut(seuratO, latent.vars = c("nUMI", "percent.mito"
  , "PC1CellCycleGenes"))

# Regress out nUMI, percet Mt, and PC1 CC genes, no row or scaling or centering
RegCcNoScCtSO <- RegressOut(seuratO, latent.vars = c("nUMI", "percent.mito"
  , "PC1CellCycleGenes"), do.scale = FALSE, do.center = FALSE)
################################################################################

### Seurat PCA, clustering,  and tSNE

# Seurat variable genes calculation
pdf(paste0(outGraphPfx, "MeanVarPlot.pdf"))
regCcSO <- MeanVarPlot(regCcSO ,fxn.x = expMean, fxn.y = logVarDivMean
  , x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = TRUE
  , plot.both = TRUE, do.text = FALSE, contour.lwd = 0.1, contour.col = "red")
dev.off()

length(regCcSO@var.genes)

# Seurat PCA
regCcSO <- PCA(regCcSO, pc.genes = regCcSO@var.genes, do.print = TRUE
  , pcs.print = 5, genes.print = 5)

pdf(paste0(outGraphPfx, "PCAgeneScores.pdf"), width = 8, height = 10)
VizPCA(regCcSO, 1:6, font.size = 0.75)
PCAPlot(regCcSO, 1, 2)
# PCHeatmap(regCcSO, pc.use = 1, do.balanced = TRUE)
PCHeatmap(regCcSO, pc.use = 1:6, do.balanced = TRUE, remove.key = FALSE
  , label.columns = FALSE, use.full = FALSE)
dev.off()

# Seurat PC elbow plot
pdf(paste0(outGraphPfx, "ElbowPlot.pdf"))
PCElbowPlot(regCcSO)
dev.off()

# Stash no cell cycle regression cluster identities for later
# Have to run these lines by themselves or gets skipped
regCcSO <- StashIdent(regCcSO, save.name = "NoCcReg")
str(regCcSO@data.info)

# Seurat tSNE and clustering
# pdf(paste0(outGraphPfx, "tSNE.pdf"))
regCcSO <- FindClusters(regCcSO, pc.use = 1:10, resolution = 0.6
  , print.output = 0, save.SNN = TRUE)
regCcSO <- RunTSNE(regCcSO, dims.use = 1:10, do.fast = TRUE)
# TSNEPlot(regCcSO)
# dev.off()

# PCA log normalized, regressed nUMI, percent mito, and PC1 CC genes, expression matrix
pcaRegCCAll <- PCA_Prcomp(seuratO@scale.data)
pdf(paste0(outGraphPfx, "pcaVar.pdf"))
plot(pcaRegCCAll, type = "line")
dev.off()
summary(pcaRegCCAll)$importance[ ,1:5]

# Correlation PCs from CC reg expression matrix and CC genes PCs
# Should this be CC genes PCs after reg CC?
cor(pcaRegCCAll$x[ ,1:5], pcaCC$x[ ,1:5])

save(regCcSO, RegCcNoScCtSO, file = "Seurat_RegCC_SeuratObj.Robj")
save(pcaAll, pcaCC, pcaRegCCAll, file = "Seurat_RegCC_PCA.rdat")
################################################################################

### Plot tSNE and Seurat clustering

# tSNE graph colored by cluster
ggTsne1 <- TSNEPlot(regCcSO, do.label = TRUE, pt.size = 0.5, do.return = TRUE)
ggTsne1 <- ggTsne1 + ggtitle(paste0("Cluster + tSNE after regressing out cell cycle"
  , "\n"))
ggTsne2 <- TSNEPlot(regCcSO, do.label = TRUE, pt.size = 0.5, do.return = TRUE
  , group.by = "NoCcReg", no.legend = TRUE)
ggTsne2 <- ggTsne2 + ggtitle(paste0("tSNE after regressing out cell cycle"
  , "\nColored by previous clustering (no cell cycle correction)"
  , "\n"))
pdf(paste0(outGraphPfx, "tSNE.pdf"), width = 10, height = 5)
MultiPlotList(list(ggTsne2, ggTsne1), cols = 2)
dev.off()
################################################################################

### Feature plot of markers

# Collect tSNE values for ggplot
ggDF <- regCcSO@tsne.rot

## Subset to marker genes of interest for Luis' excel file
# Cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
kmDFL <- split(kmDF, kmDF$Grouping)

# Loop through and plot each group of genes
pdf(paste0(outGraphPfx, "FeaturePlot_Markers.pdf"), width = 13, height = 6)
lapply(names(kmDFL), function(grouping) {
  print(grouping)
  genes <- kmDF$Gene.Symbol[kmDF$Grouping == grouping]
  ggDF <- Mean_Expression(genes)
  ggDF <- Set_Limits(ggDF)
  ggFp <- Feature_Plot(ggDF, grouping, genes)
  ggTsne <- TSNE_Plot(regCcSO, grouping)
  MultiPlotList(list(ggTsne, ggFp), cols = 2)
})
dev.off()
################################################################################

### Feature plot of Macosko cell cycle genes

# Format cell cycle genes for Feature_Plot
ccLV <- sapply(ccDF, list)
ccLV <- lapply(ccLV, function(genes) {gsub(" *", "", genes)})

# Loop through and plot each group of genes
pdf(paste0(outGraphPfx, "FeaturePlot_CCgenes.pdf"), width = 13, height = 6)
lapply(names(ccLV), function(grouping) {
  # print(grouping)
  genes <- ccLV[[grouping]]
  # Remove blanks
  genes <- genes[! genes == ""]
  print(genes)
  ggDF <- Mean_Expression(genes)
  ggDF <- Set_Limits(ggDF)
  ggFp <- Feature_Plot(ggDF, grouping, genes)
  ggTsne <- TSNE_Plot(seuratO, grouping)
  MultiPlotList(list(ggTsne, ggFp), cols = 2)
})
dev.off()
################################################################################