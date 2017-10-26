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

# Set variable to gene of interest
genes <- "CHD2"
genes <- "KDM4B"

## Input
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb

## Variables
graphCodeTitle <- "Seurat_Plot_GOI_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Seurat_Plot_GOI_DS002_003_exon_FtMm250_"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size = 12)))
theme_update(plot.title = element_text(size = 14))
################################################################################

## Assign annotated cluster names to clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
new.cluster.ids <- c("Excitatory Upper Layer Neuron 1"
  , "Excitatory Neuron"
  , "Excitatory Upper Layer Neuron 2"
  , "Excitatory Deep Layer Neuron"
  , "Intermediate Progenitors"
  , "Interneuron"
  , "Mitotic Progenitors"
  , "oRG"
  , "Oligodendrocyte Precursor"
  , "Endothelial")
seuratO@ident <- plyr::mapvalues(seuratO@ident, from = current.cluster.ids
  , to = new.cluster.ids)

plyr::mapvalues(unique(seuratO@ident), from = current.cluster.ids
  , to = new.cluster.ids)
# Stash numerical cluster identities if want to use later
seuratO <- StashIdent(seuratO, save.name = "Cluster_Numbers")

## Violin plot
ggP <- VlnPlot(seuratO, genes, do.ret = TRUE)
ggP <- ggP[[1]] + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")
ggsave(paste0(outGraphPfx, "tSNE_ViolinPlot_", genes, ".pdf"), plot = ggP)

## tSNE graph colored by cluster and feature plot

# Collect tSNE values for ggplot
ggDF <- seuratO@tsne.rot

# Subset to genes of interest
genesExL <- seuratO@scale.data[which(row.names(seuratO@scale.data) == genes), ]
# Add to ggplot data frame
ggDF$EXPRESSION <- genesExL[match(row.names(ggDF), names(genesExL))]

# ggplot
ggFp <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
  geom_point(size = 0.75) +
  scale_colour_gradient(name = "Normalized\nExpression", low = "grey"
    , high = "blue") +
  ggtitle(paste0(graphCodeTitle
    , "\n\n", genes
    , "\ntSNE plot, each point is a cell"
    , "\nColor indicates gene of interest expression"
    , "\n"))

# tSNE graph colored by cluster
ggTsne <- TSNEPlot(seuratO, do.label = T, pt.size = 0.75, do.return = TRUE)
ggTsne <- ggTsne + ggtitle(paste0(graphCodeTitle
  , "\n\n", genes
  , "\ntSNE plot, each point is a cell"
  , "\nColor indicates cluster assignment"
  , "\nClusters annotated manually by expression of known marker genes"
  , "\n"))

# Arrange tSNE and feature plot next to each other and plot
pdf(paste0(outGraphPfx, "tSNE_FeaturePlot_", genes, ".pdf"), width = 13, height = 6)
MultiPlotList(list(ggTsne, ggFp), cols = 2)
dev.off()
################################################################################