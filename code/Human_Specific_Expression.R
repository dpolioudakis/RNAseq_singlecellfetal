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
require(reshape2)
require(gridExtra)
require(ggplot2)
require(cowplot)
require(fdrtool)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# # Log normalized, regressed nUMI and percent mito
# # seuratO
# load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")
# # Log normalized, regressed nUMI and percent mito, mean centered and scaled
# # fetb
# load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

# Seurat
# PC 1-40
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Cluster DE table
deDF <- read.table(
  "../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Luis metaMat results
mmapDF <- read.csv("../source/metaMat/Overlapped-Genes.csv", header = TRUE)

# Allen Developmental Macaque human specific genes
hsDF <- read.csv("../source/Bakken_2016_AllenDevMacaque_ST10_HumanSpecific.csv"
  , header = TRUE, skip = 1)

# Known marker Luis table
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv")

## Variables
graphCodeTitle <- "Human_Specific_Expression.R"
outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_Expression_"
outTable <- "../analysis/tables/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Human_Specific_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_Expression_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

Subset_To_Specific_Clusters <- function(deDF, cluster, okayClusters, fcHigh, fcLow) {
  # Clusters gene cannot be DE in
  clsNo <- c(0:17)[! c(0:17) %in% okayClusters]
  # Gene is > X FC in cluster
  genes1 <- deDF$GENE[deDF$LOG_FC > fcHigh & deDF$CLUSTER == cluster]
  # Genes in clusters genes cannot be DE in > 0.3
  genes2 <- deDF$GENE[deDF$LOG_FC > fcLow & deDF$CLUSTER %in% clsNo]
  # Check
  print(table(genes1 %in% genes2))
  # Remove genes in clusters genes cannot be DE in
  genes1 <- genes1[! genes1 %in% genes2]
  # Filter DE DF
  utdeDF <- deDF[deDF$GENE %in% genes1, ]
  utdeDF <- utdeDF[utdeDF$CLUSTER == cluster, ]
  return(utdeDF)
}

Combine_DE_and_Expression <- function(deDF, exDF) {
  # Column for setting order of genes
  utdeDF$ORDER <- seq(1, nrow(utdeDF))
  # Merge with expression data frame
  ggDF <- merge(utdeDF[c("GENE", "CLUSTER", "ORDER")], exDF
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
  # Set order
  ggDF <- ggDF[order(-ggDF$ORDER), ]
  # Remove order variable now set
  ggDF <- ggDF[ ,! colnames(ggDF) == "ORDER"]
}
################################################################################

### metaMat cluster DE oRG human specific genes and expression ranking in RG

genes <- mmapDF[mmapDF$X == "Human-specific", "X7"]
genes <- unlist(strsplit(as.character(genes), split = "\\|"))

## Expression ranking in RG
v1 <- rowMeans(noCentExM[ ,centSO@ident %in% c(7)])
df <- data.frame(v1, rank(-v1))
df <- df[row.names(df) %in% genes, ]
df <- df[order(-df[ ,1]), ]
# PTN     3.7742088         3
# SERF2   1.0234892       205
# GNG5    0.7509851       312
# ID2     0.5910675       453
# ARL6IP5 0.5367576       507
# PTPRA   0.5233324       522
# EEF1D   0.4840341       591
# PTTG1IP 0.4240280       703
# ITGA6   0.4094136       734
# PHGDH   0.3610526       875
# LGALS1  0.3597090       881
# FAM63B  0.3568268       889
# ALDH7A1 0.3537301       901
# LYN     0.3138616      1042
# NPC2    0.2648934      1272
# VEGFA   0.2635715      1279
# UTRN    0.2339614      1491
# JPH1    0.2312890      1512

genes <- row.names(df)
################################################################################

### Human specific genes DE in oRG cluster

geneGroupDF <- data.frame(GENE = genes, GROUP = "")

ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = noCentExM
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = 0
  , upperLimit = 3
  , geneOrder = genes)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific genes DE in oRG (cluster 7)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
    , "\n")
  )
ggsave(paste0(outGraph, "HumanSpecific_Heatmap_Normalized.png")
  , width = 16, height = 10)

# Heatmap
# Normalized, mean centered and scaled
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = genes)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific genes DE in oRG (cluster 7)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n"))
ggsave(paste0(outGraph, "HumanSpecific_Heatmap_NormalizedCenteredScaled.png")
  , width = 16, height = 10)
################################################################################

### Uniquely expressed in oRG

# 18
length(genes)
# Number of DE genes 7768
nrow(deDF)
# Intersect DE gene lists and TFs, co-factors, chromatin remodelers list: 408
df <- deDF[deDF$GENE %in% genes, ]

## expressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia

# Subset to expressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia
ldf <- lapply(c(7), function(cluster) {
  Subset_To_Specific_Clusters(deDF = df, cluster = cluster
    , okayClusters = c(7, 8, 10, 11, 13, 15, 16), fcHigh = 0.2, fcLow = 0.1)
})
ssDeDF <- do.call("rbind", ldf)

# Heatmap
# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = ssDeDF$GENE, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = genes)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters: > 0.2 log fold change in cluster; < 0.1 for other clusters"
    , "\n")
  )
ggsave(paste0(
  outGraph, "DeUniqueRG_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 10, limitsize = FALSE)

# Heatmap
# Normalized
geneGroupDF <- data.frame(GENE = ssDeDF$GENE, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1
  , upperLimit = 3
  , geneOrder = genes)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters: > 0.2 log fold change in cluster; < 0.1 for other clusters"
    , "\n")
)
ggsave(paste0(
  outGraph, "DeUniqueRG_Heatmap_Normalized.png")
  , width = 12, height = 10, limitsize = FALSE)


## Feature plots

# Collect tSNE values for ggplot
tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

# Normalized
ggL <- FeaturePlot(
  genes = ssDeDF$GENE
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1
  , limHigh = 2
  , geneGrouping = NULL
  , centScale = FALSE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nNormalized expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeUniqueRG_FeaturePlot_Normalized.png")
  , width = 20, height = 35, limitsize = FALSE)

# Normalized centered scaled
ggL <- FeaturePlot(
  genes = ssDeDF$GENE
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nNormalized centered scaled expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeUniqueRG_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 20, height = 35, limitsize = FALSE)
################################################################################


DE_oRG_vs_vRG_Boxplot <- function(genes, title) {
  df <- deOvDF[row.names(deOvDF) %in% genes, ]
  df[ ,2] <- round(df[ ,2], 2)
  ggplot(df, aes(x = Gene, y = Log_Fold_Change_oRGvsvRG)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = signif(Pvalue_oRGvsvRG, 2)
      , y = Log_Fold_Change_oRGvsvRG, x = Gene), angle = 90) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Genes") +
    ylab("Log fold change") +
    ggtitle(title)
}

## oRG vs vRG

clusterIDs1 <- 7
clusterIDs2 <- 9

exDF <- centSO@data
ids <- names(centSO@ident)[centSO@ident %in% c(clusterIDs1, clusterIDs2)]
exDF <- exDF[ ,colnames(exDF) %in% ids]

# DE Linear model
termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual", "res.0.6")]
# Subset to clusters of interest
termsDF <- termsDF[termsDF$res.0.6 %in% c(clusterIDs1, clusterIDs2), ]
# Add term TRUE/FALSE cell is in cluster
termsDF$cluster <- FALSE
termsDF$cluster[termsDF$res.0.6 %in% clusterIDs1] <- TRUE
# Linear model
deLM1 <- DE_Linear_Model(
  exDatDF = exDF
  , termsDF = termsDF
  , mod = "y ~ cluster+nUMI+librarylab+individual")

# Format DE
deOvDF <- data.frame(
  Gene = row.names(deLM1$coefmat)
  , Log_Fold_Change_oRGvsvRG = deLM1$coefmat[ ,"clusterTRUE"]
  , Pvalue_oRGvsvRG = deLM1$pvalmat[ ,"clusterTRUE"]
)

# Plot fold changes of oRG marker gene lists

DE_oRG_vs_vRG_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG"]
  , title = paste0(graphCodeTitle
    , "\n\noRG known markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue")
  )
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_oRGknownMarkers.pdf")
  , width = 7, height = 6)

DE_oRG_vs_vRG_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG-PollenS3"]
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue")
)
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_oRGPollenS3.pdf")
  , width = 11, height = 7)

DE_oRG_vs_vRG_Boxplot(
  genes = c("ITGA6", "LGALS1", "LYN", "NPC2")
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue")
)
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_ITGA6_LGALS1_LYN_NPC2.pdf")
  , width = 4, height = 7)


## oRG vs Neurons, IPCs, vRGs

clusterIDs1 <- 7
clusterIDs2 <- c(0,1,2,3,4,5,6,9,12,14)

exDF <- DE_Filters_ClustersAvsB_ExpMatrix(
  so = centSO
  , minPercent = 0.1
  , foldChange = 0.2
  , clusterIDs1 = clusterIDs1
  , clusterIDs2 = clusterIDs2
)

# DE Linear model
termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual", "res.0.6")]
# Subset to clusters of interest
termsDF <- termsDF[termsDF$res.0.6 %in% c(clusterIDs1, clusterIDs2), ]
# Add term TRUE/FALSE cell is in cluster
termsDF$cluster <- FALSE
termsDF$cluster[termsDF$res.0.6 %in% clusterIDs1] <- TRUE
deLM2 <- DE_Linear_Model(
  exDatDF = exDF
  , termsDF = termsDF
  , mod = "y ~ cluster+nUMI+librarylab+individual")

df2 <- deLM1$coefmat
df2[row.names(df) %in% c("GNG5"), ]



deLM1$coefmat[ ,"clusterTRUE"]
deLM1$pvalmat[ ,"clusterTRUE"] < 0.05
corrected <- fdrtool(deLM1$pvalmat[ ,"clusterTRUE"], statistic = "pvalue", plot = FALSE)


genes <- intersect(row.names(deLM1$coefmat)[deLM1$coefmat[ ,"clusterTRUE"] > 0.3]
  , row.names(deLM2$coefmat)[deLM2$coefmat[ ,"clusterTRUE"] > 0.3]
  )

# Subset to human specific
intersect(row.names(deLM1$coefmat)[deLM1$coefmat[ ,"clusterTRUE"] > 0.2]
  , hsDF$Gene[hsDF$Set == "Human-specific"])
intersect(row.names(deLM2$coefmat)[deLM1$coefmat[ ,"clusterTRUE"] > 0.2]
  , hsDF$Gene[hsDF$Set == "Human-specific"])
intersect(genes, hsDF$Gene[hsDF$Set == "Human-specific"])
genes

# Heatmap
# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = genes, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(
  geneGroupDF = geneGroupDF
  , exprM = as.matrix(centSO@scale.data)
  , seuratO = centSO
  , clusters1 = c(0:1)
  , clusters2 = c(2:10)
  , clusters3 = c(11:17)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = genes)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters: > 0.2 log fold change in cluster; < 0.1 for other clusters"
    , "\n")
)
ggsave(paste0(
  outGraph, "DeRGspecific_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 10, limitsize = FALSE)
################################################################################
