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

DE_oRG_vs_vRG_Boxplot <- function(genes, title) {
  df <- deOvDF[row.names(deOvDF) %in% genes, ]
  df[ ,2] <- round(df[ ,2], 2)
  ggplot(df, aes(x = Gene, y = Log_Fold_Change_oRGvsvRG, fill = Human_specific)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = signif(Pvalue_oRGvsvRG, 2)
      , y = Log_Fold_Change_oRGvsvRG, x = Gene), angle = 90) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Genes") +
    ylab("Log fold change") +
    ggtitle(title)
}

DE_oRG_vs_vRG_Neuron_IPC_Boxplot <- function(genes, title) {
  df <- deOnivDF[row.names(deOnivDF) %in% genes, ]
  df[ ,2] <- round(df[ ,2], 2)
  ggplot(df, aes(x = Gene, y = Log_Fold_Change_oRG_vs_vRG_Neuron_IPC
    , fill = Human_specific)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = signif(Pvalue_oRG_vs_vRG_Neuron_IPC, 2)
      , y = Log_Fold_Change_oRG_vs_vRG_Neuron_IPC, x = Gene), angle = 90) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Genes") +
    ylab("Log fold change") +
    ggtitle(title)
}

Intersection_Of_Cells_Expressing_Both <- function(gene1, gene2){
  v1 <- noCentExM[row.names(noCentExM) %in% c(gene1), ] > 0.5
  v2 <- noCentExM[row.names(noCentExM) %in% c(gene2), ] > 0.5
  df1 <- rbind(v1, v2)
  int <- apply(df1, 2, all)
  df2 <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  df2$Intersection <- int[match(row.names(df2), names(int))]
  return(df2)
}

Intersection_tSNE_Plots <- function(genes) {
  genesM <- t(combn(genes, 2))
  genesM <- rbind(genesM, cbind(Gene1 = genes, Gene2 = genes))
  ggL <- apply(genesM, 1, function(row){
    # print(row)})
    gene1 <- row[["Gene1"]]
    gene2 <- row[["Gene2"]]
    df1 <- Intersection_Of_Cells_Expressing_Both(gene1, gene2)
    gg <- ggplot(df1, aes(x = tSNE_1, y = tSNE_2, color = Intersection)) +
      geom_point(size = 0.05) +
      scale_color_manual(name = "Intersection", values = c("grey", "red")) +
      guides(colour = guide_legend(override.aes = list(size = 7))) +
      ggtitle(paste(gene1, gene2))
    return(gg)
  })
  return(ggL)
}


Number_Of_Cells_Intersection_Heatmap <- function(genes, title){
  m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
  l1 <- apply(m1, 2, function(col) row.names(m1)[col])
  
  ldf <- lapply(l1, function(genes) {
    v2 <- genes[genes %in% genes]
    v3 <- genes[genes %in% genes]
    expand.grid(v2, v3)
  })
  df1 <- do.call("rbind", ldf)
  
  # Format for matrix
  df3 <- dcast(df1, Var1 ~ Var2)
  row.names(df3) <- df3$Var1
  df3 <- df3[ ,-1]
  
  # Format for ggplot
  df3$Gene <- row.names(df3)
  df3 <- melt(df3)
  df3$Gene <- factor(df3$Gene, levels = genes)
  df3$variable <- factor(df3$variable, levels = genes)
  
  # Plot counts of intersections
  gg <- ggplot(df3, aes(x = Gene, y = variable, fill = value)) +
    geom_tile() +
    geom_text(data = df3, aes(x = Gene, y = variable, label = value), color = "black") +
    scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "Number") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Genes") +
    ylab("Genes") +
    ggtitle(title)
  
  return(gg)
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

### Human specific genes that are DE oRG vs vRG and DE oRG vs neuron, IPC, vRG

## oRG vs vRG

# Clusters of interest
clusterIDs1 <- 7
clusterIDs2 <- 9

# Subset exDF to clusters of interest
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

# Add human specific
deOvDF$Human_specific <- FALSE
deOvDF$Human_specific[
  deOvDF$Gene %in% hsDF$Gene[hsDF$Set == "Human-specific"]] <- TRUE

# Plot fold changes of oRG marker gene lists

DE_oRG_vs_vRG_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG"]
  , title = paste0(graphCodeTitle
    , "\n\noRG known markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
  )
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_oRGknownMarkers.pdf")
  , width = 7, height = 6)

DE_oRG_vs_vRG_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG-PollenS3"]
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_oRGPollenS3.pdf")
  , width = 11, height = 7)

DE_oRG_vs_vRG_Boxplot(
  genes = c("ITGA6", "LGALS1", "LYN", "NPC2")
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_ITGA6_LGALS1_LYN_NPC2.pdf")
  , width = 4, height = 7)

DE_oRG_vs_vRG_Boxplot(
  genes = deOvDF$Gene[deOvDF$Human_specific == TRUE]
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG (cluster 9)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRG_intersect_HumanSpecific.pdf")
  , width = 28, height = 7)


## oRG vs Neurons, IPCs, vRGs

# Clusters of interest
clusterIDs1 <- 7
clusterIDs2 <- c(0,1,2,3,4,5,6,9,12,14)

# Subset exDF to clusters of interest
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
deLM2 <- DE_Linear_Model(
  exDatDF = exDF
  , termsDF = termsDF
  , mod = "y ~ cluster+nUMI+librarylab+individual")


# Format DE
deOnivDF <- data.frame(
  Gene = row.names(deLM2$coefmat)
  , Log_Fold_Change_oRG_vs_vRG_Neuron_IPC = deLM2$coefmat[ ,"clusterTRUE"]
  , Pvalue_oRG_vs_vRG_Neuron_IPC = deLM2$pvalmat[ ,"clusterTRUE"]
)
# Save as csv
write.csv(deOnivDF, file = paste0(outTable, "DE_oRG_vs_Neuron_IPC_vRG.csv")
  , quote = FALSE)
deOnivDF <- read.csv(paste0(outTable, "DE_oRG_vs_Neuron_IPC_vRG.csv")
  , header = TRUE, row.names = 1)

# Add human specific
deOnivDF$Human_specific <- FALSE
deOnivDF$Human_specific[
  deOnivDF$Gene %in% hsDF$Gene[hsDF$Set == "Human-specific"]] <- TRUE

# Plot fold changes of oRG marker gene lists

DE_oRG_vs_vRG_Neuron_IPC_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG"]
  , title = paste0(graphCodeTitle
    , "\n\noRG known markers"
    , "\n\nDE of oRG (cluster 7) vs vRG, neuron, IPC clusters (clusters 0,1,2,3,4,5,6,9,12,14)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRGneuronIPC_intersect_oRGknownMarkers.pdf")
  , width = 7, height = 6)

DE_oRG_vs_vRG_Neuron_IPC_Boxplot(
  genes = kmDF$Gene.Symbol[kmDF$Grouping == "oRG-PollenS3"]
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG, neuron, IPC clusters (clusters 0,1,2,3,4,5,6,9,12,14)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRGneuronIPC_intersect_oRGPollenS3.pdf")
  , width = 11, height = 7)

DE_oRG_vs_vRG_Neuron_IPC_Boxplot(
  genes = c("ITGA6", "LGALS1", "LYN", "NPC2")
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG, neuron, IPC clusters (clusters 0,1,2,3,4,5,6,9,12,14)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRGneuronIPC_intersect_ITGA6_LGALS1_LYN_NPC2.pdf")
  , width = 4, height = 7)

DE_oRG_vs_vRG_Neuron_IPC_Boxplot(
  genes = deOnivDF$Gene[deOnivDF$Human_specific == TRUE]
  , title = paste0(graphCodeTitle
    , "\n\noRG Pollen Table S3 markers"
    , "\n\nDE of oRG (cluster 7) vs vRG, neuron, IPC clusters (clusters 0,1,2,3,4,5,6,9,12,14)"
    , "\nLog fold change (oRG vs vRG)"
    , "\nText indicates pvalue"
    , "\nHuman specific: Allen list")
)
ggsave(paste0(outGraph, "DE_oRGvsvRGneuronIPC_intersect_HumanSpecific.pdf")
  , width = 28, height = 7)

# Intersect DE oRG vs vRG, DE oRG vs neuron IPC vRG, and human specific genes
genes <- intersect(row.names(deOvDF)[deOvDF$Log_Fold_Change_oRGvsvRG > 0.15]
  , row.names(deOnivDF)[deOnivDF$Log_Fold_Change_oRG_vs_vRG_Neuron_IPC > 0.3]
)
genes <- intersect(genes, hsDF$Gene[hsDF$Set == "Human-specific"])

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
    , "\n\nExpression of Allen human specific differentially expressed in oRG cluster"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.15 log fold change oRG vs vRG cluster"
    , "\n> 0.3 log fold change oRG vs vRG, neuron, IPC clusters"
    , "\n")
)
ggsave(paste0(
  outGraph, "DeRGspecific015_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 15, limitsize = FALSE)

# Intersect DE oRG vs vRG, DE oRG vs neuron IPC vRG, and human specific genes
genes <- intersect(row.names(deOvDF)[deOvDF$Log_Fold_Change_oRGvsvRG > 0.3]
  , row.names(deOnivDF)[deOnivDF$Log_Fold_Change_oRG_vs_vRG_Neuron_IPC > 0.3]
  )
genes <- intersect(genes, hsDF$Gene[hsDF$Set == "Human-specific"])

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
    , "\n\nExpression of Allen human specific differentially expressed in oRG cluster"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.3 log fold change oRG vs vRG cluster"
    , "\n> 0.3 log fold change oRG vs vRG, neuron, IPC clusters"
    , "\n")
)
ggsave(paste0(
  outGraph, "DeRGspecific_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 15, limitsize = FALSE)


# Feature plot

# Collect tSNE values for ggplot
tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

# Normalized
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1
  , limHigh = 2
  , geneGrouping = NULL
  , centScale = FALSE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nDE filters:"
    , "\n> 0.3 log fold change oRG vs vRG cluster"
    , "\n> 0.3 log fold change oRG vs vRG, neuron, IPC clusters"
    , "\nNormalized expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeRGspecific_FeaturePlot_Normalized.png")
  , width = 16, height = 16, limitsize = FALSE)

# Normalized centered scaled
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(
  ggPlotsL = ggL, ncol = 2, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific differentially expressed in RG clusters"
    , "\nExpressed in oRG, and may be in cycling cells, endo, OPCs, pericytes, microglia"
    , "\nDE filters:"
    , "\n> 0.3 log fold change oRG vs vRG cluster"
    , "\n> 0.3 log fold change oRG vs vRG, neuron, IPC clusters"
    , "\nNormalized centered scaled expression"
    , "\n")
)
ggsave(paste0(outGraph, "DeRGspecific_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 16, height = 16, limitsize = FALSE)


## tSNE colored by intersection and heatmap of numbers of intersections

# RG
genes <- c(genes, "ITGA6", "SOX2", "PAX6", "HOPX", "CRYAB")

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"
    , "\nRG"))
)
ggsave(paste0(outGraph, "DeRGspecific_And_Marker_Intersection_tSNE.png")
  , width = 20, height = 75, limitsize = FALSE)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\nRG")
)
ggsave(paste0(outGraph, "DeRGspecific_And_Markers_NumberIntersect_Heatmap.png")
  , width = 7, height = 7)

# Violin plot of genes of interest
Gene_Expression_By_Cluster_ViolinPlot(genes = genes
  , exprM = noCentExM, clusterIDs = centSO@ident
  , geneOrder = genes, grouping = NULL)
ggsave(paste0(outGraph, "DeRGspecific_ExprViolinPlot.png")
  , width = 13, height = 10)
################################################################################