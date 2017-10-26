# Damon Polioudakis
# 2017-05-28
# Clustering of Drop-seq cells by digital gene expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(reshape2)
require(ggplot2)
require(cowplot)
require(WGCNA)
require(pheatmap)

options(stringsAsFactors = FALSE)

# Regress CC score
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
regCcSO <- centSO
regCcExM <- noCentExM
rm(noCentExM)
rm(centSO)
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# regCcSO <- ssCentSO
# regCcExM <- ssNoCentExM

# Remove CC genes from variable gene list used for clustering
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")
rmCcSO <- centSO
rmCcExM <- noCentExM
rm(noCentExM)
rm(centSO)
# load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_TEST_seuratO.Robj")
# rmCcSO <- ssCentSO
# rmCcExM <- ssNoCentExM

# Keep CC genes from variable gene list used for clustering
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
keepCcSO <- centSO
keepCcExM <- noCentExM
rm(noCentExM)
rm(centSO)
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# keepCcSO <- ssCentSO
# keepCcExM <- ssNoCentExM

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccGenes <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
sGenes <- ccGenes[1:43]
g2mGenes <- ccGenes[44:98]

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_Cluster_CellCycleCorrection_DS2-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11_CellCycleCorrection/Seurat_Cluster_CellCycleCorrection_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11_CellCycleCorrection/Seurat_Cluster_CellCycleCorrection_DS2-11_"
outData <- "../analysis/Seurat_Cluster_DS2-11_CellCycleCorrection/Seurat_Cluster_CellCycleCorrection_DS2-11_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outData)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

Seurat_CCexpression_by_Cluster_Boxplot <- function (exM, seuratO, title) {
  l1 <- colMeans(exM[row.names(exM) %in% sGenes, ])
  df <- data.frame(S = l1)
  l2 <- colMeans(exM[row.names(exM) %in% g2mGenes, ])
  df$G2M <- l2
  df$CLUSTER <- seuratO@ident[match(row.names(df), names(seuratO@ident))]
  df <- melt(df)
  gg <- ggplot(df, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_boxplot() +
    scale_fill_discrete(name = "Phase") +
    ylab("Mean log2 normalized expression") +
    xlab("Clusters") +
    ggtitle(title)
  return(gg)
}

Seurat_CC_by_Cluster_Barplot <- function (seuratO, title) {
  df <- seuratO@meta.data[ ,c("Phase", "res.0.6")]
  df <- do.call("rbind", tapply(df$Phase, df$res.0.6, table))
  df <- melt(df)
  gg <- ggplot(df, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity") +
    ylab("Number of cells") +
    xlab("Clusters") +
    ggtitle(title)
  return(gg)
}

Jaccard_Index <- function(v1, v2) {
  length(intersect(v1, v2)) / (length(union(v1, v2)))
}

Seurat_Cluster_Compare_Jaccard_Index <- function (
  seuratO1, seuratO2, label1, label2) {
  clids1 <- as.numeric(as.character(sort(unique(seuratO1@ident))))
  clids2 <- as.numeric(as.character(sort(unique(seuratO2@ident))))
  
  # Empty matrix
  jiM <- matrix(NA
    , length(clids1)
    , length(clids2))
  
  # Column and row names are cluster IDs
  row.names(jiM) <- paste(clids1, label1)
  colnames(jiM) <- paste(clids2, label2)
  
  # Fill with Jaccard index
  for (i in 1:length(clids1)){
    for (j in 1:length(clids2)){
      clid1 <- clids1[i]
      clid2 <- clids2[j]
      v1 <- names(seuratO1@ident)[seuratO1@ident == clid1]
      v2 <- names(seuratO2@ident)[seuratO2@ident == clid2]
      jiM[i,j] <- Jaccard_Index(v1, v2)
    }
  }
  
  return(jiM)
}

Seurat_Cluster_Compare_Jaccard_Index_Heatmap <- function (
  jiM, xlab, ylab, title = NULL) {
  
  # Format
  df <- melt(jiM)
  df$value <- round(df$value, 3)
  
  # Plot
  ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(axis.text.y = element_text(size = 8)) +
    scale_fill_distiller(
      name = "Jaccard index", type = "seq", palette = 1, direction = 1) +
    geom_text(data = df, aes(x = Var2, y = Var1, label = value)
      , color = "black") +
    scale_x_discrete(breaks = unique(df$Var2)) +
    #   , limits = c(0, length(unique(df$Var2)))) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)
}

Seurat_Cluster_Compare_Percent_Overlap <- function (
  seuratO1, seuratO2, label1, label2) {
  clids1 <- as.numeric(as.character(sort(unique(seuratO1@ident))))
  clids2 <- as.numeric(as.character(sort(unique(seuratO2@ident))))
  
  # Empty matrix
  jiM <- matrix(NA
    , length(clids1)
    , length(clids2))
  
  # Column and row names are cluster IDs
  row.names(jiM) <- paste(clids1, label1)
  colnames(jiM) <- paste(clids2, label2)
  
  # Fill with Jaccard index
  for (i in 1:length(clids1)){
    for (j in 1:length(clids2)){
      clid1 <- clids1[i]
      clid2 <- clids2[j]
      v1 <- names(seuratO1@ident)[seuratO1@ident == clid1]
      v2 <- names(seuratO2@ident)[seuratO2@ident == clid2]
      jiM[i,j] <- (sum(v1 %in% v2) / length(v1)) * 100
    }
  }
  
  return(jiM)
}

Seurat_Cluster_Compare_Percent_Overlap_Heatmap <- function (
  olM, xlab, ylab, title) {
  
  # Format
  df <- melt(olM)
  df$value <- round(df$value, 1)
  
  # Plot
  ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(axis.text.y = element_text(size = 8)) +
    scale_fill_distiller(
      name = "Percent overlap", type = "seq", palette = 1, direction = 1) +
    geom_text(data = df, aes(x = Var2, y = Var1, label = value)
      , color = "black") +
    scale_x_discrete(breaks = unique(df$Var2)) +
    #   , limits = c(0, length(unique(df$Var2)))) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle("Percent of each column cluster in each row cluster")
}

Mean_Gene_Group_Expression <- function(exM, grouping) {
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% grouping]
  colMeans(
    exM[row.names(exM) %in% genes, ]
  )
}

Expression_Subset <- function (mnExM, threshold) {
  df <- mnExM
  df[ ,1:4] <- df[ ,1:4] > threshold
  df <- melt(df, measure.vars = colnames(df[1:4]))
  # Number of genes above threshold
  df <- aggregate(df$value, list(df$variable, df$CLUSTER), sum)
  # Factor levels for ggplot
  df$Group.2 <- factor(df$Group.2
    , levels = sort(unique(as.numeric(df$Group.2))))
  return(df)
}

Format_Number_Cell_Types_Cluster_Dataframe <- function (exM, seuratO) {
  
  mnExDF <- data.frame(
    vRG = Mean_Gene_Group_Expression(exM = exM, grouping = "vRG")
    , oRG = Mean_Gene_Group_Expression(exM = exM, grouping = "oRG")
    , RG = Mean_Gene_Group_Expression(exM = exM, grouping = "RG")
    , IP = Mean_Gene_Group_Expression(exM = exM, grouping = "IP")
    , Endothelial = Mean_Gene_Group_Expression(exM = exM, grouping = "Endothelial Cell")
    , Neuron = Mean_Gene_Group_Expression(exM = exM, grouping = "Neuron")
  )
  
  idx <- match(row.names(mnExDF), row.names(seuratO@meta.data))
  mnExDF$PHASE <- seuratO@meta.data$Phase[idx]
  mnExDF$CLUSTER <- seuratO@meta.data$res.0.6[idx]
  mnExDF$nUMI <- seuratO@meta.data$nUMI[idx]
  
  return(mnExDF)
}

Number_Cell_Types_Cluster_Barplot <- function(
  exM, seuratO, threshold, title) {
  
  df <- Format_Number_Cell_Types_Cluster_Dataframe(exM, seuratO)
  df <- Expression_Subset(df, threshold = threshold)
  
  gg <- ggplot(df, aes(x = Group.2, y = x, fill = Group.1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab("Cluster") +
    ylab("Number of cells") +
    ggtitle(title)
  
  return(gg)
}

Positive_Negative_Expression_Flag <- function(
  exDF, highThreshold, lowThreshold) {
  
  # Expected input data frame format
  # (from Format_Number_Cell_Types_Cluster_Dataframe)
  # vRG         oRG          RG          IP Endothelial
  # ACCAGCTAGCCT -0.12961098 -0.04274303 -0.16390677 -0.05747777  0.08723860
  # GGAAGGACTGCA  0.17994922 -0.08862593  0.47211654 -0.02872334  1.65468242
  # Neuron PHASE CLUSTER nUMI
  # ACCAGCTAGCCT  0.8625643    G1       3 6165
  # GGAAGGACTGCA -0.9443603   G2M      11 4166
  
  df <- exDF
  df$TYPE <- NA
  df$TYPE[df[ ,c("IP")] > highThreshold] <- "IP+"
  df$TYPE[df[ ,c("RG")] > highThreshold &
      apply((df[ ,c("vRG", "oRG", "IP")] < lowThreshold), 1, all)] <- "RG+ vRG- oRG- IP-"
  df$TYPE[df[ ,c("vRG")] > highThreshold &
      df[ ,c("IP")] < lowThreshold] <- "vRG+ IP-"
  df$TYPE[df[ ,c("oRG")] > highThreshold &
      df[ ,c("IP")] < lowThreshold] <- "oRG+ IP-"
  df$TYPE[df[ ,c("IP")] > highThreshold &
      apply((df[ ,c("vRG", "oRG", "RG")] < lowThreshold), 1, all)] <- "IP+ vRG- oRG- RG-"
  df$TYPE[apply((df[ ,c("IP", "RG")] > highThreshold), 1, all)] <- "IP+ RG+"
  # Controls
  df$TYPE[df[ ,c("Endothelial")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Endothelial+ IP+"
  df$TYPE[df[ ,c("Neuron")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Neuron+ IP+"
  
  df$TYPE <- factor(df$TYPE, levels = c("IP+", "RG+ vRG- oRG- IP-", "vRG+ IP-"
    , "oRG+ IP-", "IP+ vRG- oRG- RG-", "IP+ RG+", "Endothelial+ IP+", 
    "Neuron+ IP+"))
  df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))
  
  return(df)
}

Number_RgIp_Cluster_Barplot <- function(
  exM, seuratO, highThreshold, lowThreshold, title) {
  
  df <- Format_Number_Cell_Types_Cluster_Dataframe(exM, seuratO)
  df <- Positive_Negative_Expression_Flag(
    exDF = df, highThreshold = 0.5, lowThreshold = 0.5)
  
  gg <- ggplot(df, aes(x = CLUSTER, fill = TYPE)) +
    geom_bar(stat = "count") +
    xlab("Cluster") + 
    ylab("Number of cells") +
    ggtitle(title)
  
  return(gg)
}

Percent_Of_Table <- function (tableArg) {
  tb <- table(tableArg, exclude = NULL) / sum(table(tableArg, exclude = NULL)) * 100
  v1 <- data.frame(tb)$Freq
  return(tb)
}

Percent_RgIpCluster_Barplot <- function(
  exM, seuratO, highThreshold, lowThreshold, title) {
  
  df <- Format_Number_Cell_Types_Cluster_Dataframe(exM, seuratO)
  df <- Positive_Negative_Expression_Flag(
    exDF = df, highThreshold = highThreshold, lowThreshold = lowThreshold)
  df <- aggregate(df$TYPE, list(df$CLUSTER), Percent_Of_Table)
  df2 <- as.data.frame(df[ ,2])
  df2$CLUSTER <- df[ ,1]
  print(df2)
  df <- melt(df2)
  
  gg <- ggplot(df, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    xlab("Cluster") + 
    ylab("Percent of cells") +
    ggtitle(title)
  
  return(gg)
}


Positive_Negative_Expression_Subset <- function(
  exDF, highThreshold, lowThreshold) {
  
  # Expected input data frame format
  # vRG         oRG          RG          IP Endothelial
  # ACCAGCTAGCCT -0.12961098 -0.04274303 -0.16390677 -0.05747777  0.08723860
  # GGAAGGACTGCA  0.17994922 -0.08862593  0.47211654 -0.02872334  1.65468242
  # Neuron PHASE CLUSTER nUMI
  # ACCAGCTAGCCT  0.8625643    G1       3 6165
  # GGAAGGACTGCA -0.9443603   G2M      11 4166
  
  df <- exDF
  
  Threshold_High <- function(groupings, highThreshold, anyOrAll) {
    passFilter <- apply((data.frame(df[ ,groupings] > highThreshold))
      , 1, anyOrAll)
    return(passFilter)
  }
  
  Threshold_Low <- function(groupings, lowThreshold, anyOrAll) {
    passFilter <- apply((data.frame(df[ ,groupings] < lowThreshold))
      , 1, anyOrAll)
    return(passFilter)
  }
  
  ldf <- list(
    data.frame(
      nUMI = df$nUMI[df[ ,c("IP")] > highThreshold]
      , SUBSET = "IP+")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("RG"), highThreshold, all) &
          Threshold_Low(groupings = c("vRG", "oRG", "IP"), lowThreshold, all)]
      , SUBSET = "RG+ vRG- oRG- IP-")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("vRG"), highThreshold, all) &
          Threshold_Low(groupings = c("IP"), lowThreshold, all)]
      , SUBSET = "vRG+ IP-")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("oRG"), highThreshold, all) &
          Threshold_Low(groupings = c("IP"), lowThreshold, all)]
      , SUBSET = "oRG+ IP-")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("IP"), highThreshold, all) &
          Threshold_Low(groupings = c("vRG", "oRG", "RG"), lowThreshold, all)]
      , SUBSET = "IP+ vRG- oRG- RG-")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("IP", "RG"), highThreshold, all)]
      , SUBSET = "IP+ RG+")
    
    , data.frame(
      nUMI = df$nUMI[
        Threshold_High(groupings = c("IP", "Neuron"), highThreshold, all)]
      , SUBSET = "Neuron+ IP+")
    
    , data.frame(
      nUMI = df$nUMI
      , SUBSET = "All cells")
    
    # , data.frame(
    #   nUMI = df$nUMI[
    #     Threshold_High(groupings = c("IP", "Endothelial"), highThreshold, all)]
    #   , SUBSET = "Endothelial+ IP+")
  )
  
  # df$TYPE <- factor(df$TYPE, levels = c("IP+", "RG+ vRG- oRG- IP-", "vRG+ IP-"
  #   , "oRG+ IP-", "IP+ vRG- oRG- RG-", "IP+ RG+", "Endothelial+ IP+",
  #   "Neuron+ IP+"))
  
  df <- do.call("rbind", ldf)
  return(df)
}

nUMI_RgIP_Boxplot <- function(exM, seuratO, title) {
  df <- Format_Number_Cell_Types_Cluster_Dataframe(exM, seuratO)
  df <- Positive_Negative_Expression_Subset(
    exDF = df, highThreshold = 0.5, lowThreshold = 0.5)
  gg <- ggplot(df, aes(x = SUBSET, y = nUMI)) +
    geom_boxplot(aes(fill = SUBSET)) +
    geom_jitter(aes(x = SUBSET, y = nUMI)
      , size = 0.02, height = 0, alpha = 0.2) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
  return(gg)
}
################################################################################

### Cell Cycle phase plots

# Add phase info to remove CC genes Seurat object
# (have not rerun calculating CC phase for this clustering yet)
rmCcSO@meta.data$Phase <- keepCcSO@meta.data$Phase


## PCA using CC genes

# PCA
tmp <- RunPCA(object = rmCcSO, pc.genes = c(sGenes, g2mGenes), do.print = FALSE)
gg1 <- PCAPlot(object = tmp, do.return = TRUE, group.by = "Phase", pt.size = 0.1) +
  ggtitle("Remove CC genes")
rm(tmp)
tmp <- RunPCA(object = keepCcSO, pc.genes = c(sGenes, g2mGenes), do.print = FALSE)
gg2 <- PCAPlot(object = tmp, do.return = TRUE, group.by = "Phase", pt.size = 0.1) +
  ggtitle("Keep CC genes")
rm(tmp)
tmp <- RunPCA(object = regCcSO, pc.genes = c(sGenes, g2mGenes), do.print = FALSE)
gg3 <- PCAPlot(object = tmp, do.return = TRUE, group.by = "Phase", pt.size = 0.1) +
  ggtitle("Regress CC")
rm(tmp)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPCA using CC genes"
  , "\n\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(outGraph, "PCA_CCgenes.png")
  , width = 14, height = 6)


## Expression of CC genes by cluster
gg1 <- Seurat_CCexpression_by_Cluster_Boxplot(keepCcExM, keepCcSO, "Keep CC genes")
gg2 <- Seurat_CCexpression_by_Cluster_Boxplot(rmCcExM, rmCcSO, "Remove CC genes")
gg3 <- Seurat_CCexpression_by_Cluster_Boxplot(regCcExM, regCcSO, "Regress out CC")
pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nExpression of cell cycle genes by cluster"
  , "\n\nMean expression of cell cycle genes calculated for each cell"
  , "\n\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(outGraph, "CCphase_Cluster_Expression_Barplot.pdf")
  , width = 15, height = 7)

## Expression of CC genes by cluster using remove CC genes clusters
gg1 <- Seurat_CCexpression_by_Cluster_Boxplot(keepCcExM, rmCcSO, "Keep CC genes")
gg2 <- Seurat_CCexpression_by_Cluster_Boxplot(rmCcExM, rmCcSO, "Remove CC genes")
gg3 <- Seurat_CCexpression_by_Cluster_Boxplot(regCcExM, rmCcSO, "Regress out CC")
pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nExpression of cell cycle genes by cluster (using remove CC clusters)"
  , "\n\nMean expression of cell cycle genes calculated for each cell"
  , "\n\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(outGraph, "CCphase_ClusterRemoveCC_Expression_Barplot.pdf")
  , width = 15, height = 7)


## Number of cells assigned to each CC phase
df <- data.frame(table(keepCcSO@meta.data$Phase))
# Check numbers are the same for other CC correction Seurat objects
table(regCcSO@meta.data$Phase)
# Plot
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = Var1, y = Freq+300)) +
  ylab("Number of cells in CC phase") +
  xlab("Cell cycle phase") +
  ggtitle(paste0(graphCodeTitle, 
    "\n\nNumber of cells assigned to each cell cycle phase"))
ggsave(paste0(outGraph, "CCphase_Barplot.pdf"), width = 5, height = 6)


## Number of cells assigned to each CC phase for each cluster
gg1 <- Seurat_CC_by_Cluster_Barplot(keepCcSO, "Keep CC genes")
gg2 <- Seurat_CC_by_Cluster_Barplot(rmCcSO, "Remove CC genes")
gg3 <- Seurat_CC_by_Cluster_Barplot(regCcSO, "Regress out CC")
pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nNumber of cells assigned to cell cycle phase by cluster"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(outGraph, "CCphase_Cluster_Barplot.pdf"), width = 13, height = 7)

## tSNE colored by CC phase
# tSNE values from Seurat object
ldf <- list(
  as.data.frame(keepCcSO@dr$tsne@cell.embeddings)
  , as.data.frame(rmCcSO@dr$tsne@cell.embeddings)
  , as.data.frame(regCcSO@dr$tsne@cell.embeddings)
)
# Add CC phase
ldf <- lapply(ldf, function(df) {
  df$PHASE <- keepCcSO@meta.data$Phase
  return(df)})
names(ldf) <- c("Keep CC genes", "Remove CC genes", "Regress out CC")
# Plot tSNE from each cell cycle correction tested colored by CC phase 
ggL <- lapply(names(ldf), function(name) {
  df <- ldf[[name]]
  # Plot
  gg <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, col = PHASE)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggtitle(name)
  return(gg)
})
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE colored by cell cycle phase"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.35, 1))
ggsave(paste0(outGraph, "tSNE_CCphase.png"), width = 16, height = 7)
################################################################################

### Compare clustering

## tSNE colored by clustering

# Remove CC genes tSNE
ggDF <- as.data.frame(rmCcSO@dr$tsne@cell.embeddings)
ggDF$RM_CC_CLUSTER <- rmCcSO@ident

idx <- match(names(keepCcSO@ident), row.names(ggDF))
ggDF$KEEP_CC_CLUSTER <- keepCcSO@ident[idx]

idx <- match(names(regCcSO@ident), row.names(ggDF))
ggDF$REG_CC_CLUSTER <- regCcSO@ident[idx]

# Keep CC genes tSNE
gg2DF <- as.data.frame(keepCcSO@dr$tsne@cell.embeddings)
gg2DF$KEEP_CC_CLUSTER <- keepCcSO@ident

idx <- match(names(rmCcSO@ident), row.names(gg2DF))
gg2DF$RM_CC_CLUSTER <- rmCcSO@ident[idx]

idx <- match(names(regCcSO@ident), row.names(gg2DF))
gg2DF$REG_CC_CLUSTER <- regCcSO@ident[idx]

# Regress CC genes tSNE
gg3DF <- as.data.frame(regCcSO@dr$tsne@cell.embeddings)
gg3DF$REG_CC_CLUSTER <- regCcSO@ident

idx <- match(names(rmCcSO@ident), row.names(gg3DF))
gg3DF$RM_CC_CLUSTER <- rmCcSO@ident[idx]

idx <- match(names(keepCcSO@ident), row.names(gg3DF))
gg3DF$KEEP_CC_CLUSTER <- keepCcSO@ident[idx]


gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = RM_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Remove CC genes tSNE, Remove CC genes clustering")

gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = KEEP_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Remove CC genes tSNE, Keep CC genes clustering")

gg3 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = REG_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Remove CC genes tSNE, Regress CC clustering")

gg4 <- ggplot(gg2DF, aes(x = tSNE_1, y = tSNE_2, col = RM_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Keep CC genes tSNE, Remove CC genes clustering")

gg5 <- ggplot(gg2DF, aes(x = tSNE_1, y = tSNE_2, col = KEEP_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Keep CC genes tSNE, Keep CC genes clustering")

gg6 <- ggplot(gg2DF, aes(x = tSNE_1, y = tSNE_2, col = REG_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Keep CC genes tSNE, Regress CC clustering")

gg7 <- ggplot(gg3DF, aes(x = tSNE_1, y = tSNE_2, col = RM_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Reg CC genes tSNE, Remove CC genes clustering")

gg8 <- ggplot(gg3DF, aes(x = tSNE_1, y = tSNE_2, col = KEEP_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Reg CC genes tSNE, Keep CC genes clustering")

gg9 <- ggplot(gg3DF, aes(x = tSNE_1, y = tSNE_2, col = REG_CC_CLUSTER)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 7))) +
  theme(legend.title = element_blank()) +
  ggtitle("Reg CC genes tSNE, Regress CC clustering")

pg <- plot_grid(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE colored by clustering after different cell cycle corrections"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "tSNE_RemCC_RegCC.png"), width = 21, height = 20)


## Jacard index


# Keep CC and Remove CC
m1 <- Seurat_Cluster_Compare_Jaccard_Index(
  rmCcSO, keepCcSO, "RemoveCC", "KeepCC")
gg1 <- Seurat_Cluster_Compare_Jaccard_Index_Heatmap(
  jiM = m1
  , xlab = "Keep CC"
  , ylab = "Remove CC"
)

# Remove CC and Regress CC
m2 <- Seurat_Cluster_Compare_Jaccard_Index(
  rmCcSO, regCcSO, "RemoveCC", "RegCC")
gg2 <- Seurat_Cluster_Compare_Jaccard_Index_Heatmap(
  jiM = m2
  , xlab = "RegCC"
  , ylab = "Remove CC"
)

# Keep CC and Regress CC
m3 <- Seurat_Cluster_Compare_Jaccard_Index(
  keepCcSO, regCcSO, "KeepCC", "RegCC")
gg3 <- Seurat_Cluster_Compare_Jaccard_Index_Heatmap(
  jiM = m3
  , xlab = "RegCC"
  , ylab = "KeepCC"
)

# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 1)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nJaccard index cluster overlap after different cell cycle corrections"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
# Save
ggsave(paste0(outGraph, "Jaccard_Heatmap_RemoveVsKeepCC.pdf")
  , height = 38, width = 12)


## Percent overlap of cells in cluster A versus cluster B

# Keep CC and Remove CC
m1 <- Seurat_Cluster_Compare_Percent_Overlap(
  rmCcSO, keepCcSO, "RemoveCC", "KeepCC")
gg1 <- Seurat_Cluster_Compare_Percent_Overlap_Heatmap(
  olM = m1
  , xlab = "Keep CC"
  , ylab = "Remove CC"
)

# Remove CC and Regress CC
m2 <- Seurat_Cluster_Compare_Percent_Overlap(
  rmCcSO, regCcSO, "RemoveCC", "RegCC")
gg2 <- Seurat_Cluster_Compare_Percent_Overlap_Heatmap(
  olM = m2
  , xlab = "RegCC"
  , ylab = "Remove CC"
)

# Keep CC and Regress CC
m3 <- Seurat_Cluster_Compare_Percent_Overlap(
  keepCcSO, regCcSO, "KeepCC", "RegCC")
gg3 <- Seurat_Cluster_Compare_Percent_Overlap_Heatmap(
  olM = m3
  , xlab = "RegCC"
  , ylab = "KeepCC"
)

# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 1)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent overlap of clusters after different cell cycle corrections"
  , "\nPercent of each column cluster in each row cluster"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
# Save
ggsave(paste0(outGraph, "PercentOverlap_Heatmap_RemoveVsKeepCC.pdf")
  , height = 32, width = 10)
################################################################################

### Number of cell types in CC phase by cluster

## Number cell types per cluster barplot

gg1 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 1
  , title = "Remove CC genes\n>1 normalized expression")

gg2 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 0.5
  , title = "Remove CC genes\n>0.5 normalized expression")

gg3 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 1
  , title = "Keep CC genes\n>1 normalized expression")

gg4 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 0.5
  , title = "Keep CC genes\n>0.5 normalized expression")

gg5 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 1
  , title = "Regress CC\n>1 normalized expression")

gg6 <- Number_Cell_Types_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, threshold = 0.5
  , title = "Regress CC\n>0.5 normalized expression")

Plot_Grid <- function(ggPlotsL, ncol, title, rel_height) {
  # Plot grid
  pg <- plot_grid(plotlist = ggPlotsL, ncol = 2)
  # now add the title
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(rel_height, 1))
}
# Plot grid (cowplot)
Plot_Grid(list(gg1, gg2, gg3, gg4, gg5, gg6), ncol = 2, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells above expression threshold of mean expression of marker genes"
    , "\nPercent of each column cluster in each row cluster"
    , "\nCell cycle corrections tested:"
    , "\nUse cell cycle genes for clustering"
    , "\nRemove cell cycle genes from genes used for clustering"
    , "\nRegress out cell cycle scores"
    , "\n"))
ggsave(paste0(outGraph, "NumberCellTypesCluster_Barplot.pdf")
  , width = 13, height = 13)


## Number of cells passing expression filters per cluster barplots

gg1 <- Number_RgIp_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, highThreshold = 0.5, lowThreshold = 0.5
  , title = "Remove CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Number_RgIp_Cluster_Barplot(
  exM = rmCcExM, seuratO = rmCcSO, highThreshold = 0.5, lowThreshold = 0.25
  , title = "Remove CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Number_RgIp_Cluster_Barplot(
  exM = keepCcExM, seuratO = keepCcSO, highThreshold = 0.5, lowThreshold = 0.5
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg4 <- Number_RgIp_Cluster_Barplot(
  exM = keepCcExM, seuratO = keepCcSO, highThreshold = 0.5, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg5 <- Number_RgIp_Cluster_Barplot(
  exM = regCcExM, seuratO = regCcSO, highThreshold = 0.5, lowThreshold = 0.5
  , title = "Regress CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg6 <- Number_RgIp_Cluster_Barplot(
  exM = regCcExM, seuratO = regCcSO, highThreshold = 0.5, lowThreshold = 0.25
  , title = "Regress CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid (cowplot)
Plot_Grid(list(gg1, gg2, gg3, gg4, gg5, gg6), ncol = 2, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells passing combinations of RG, oRG, vRG, IPC expression filters"
    , "\nCell cycle corrections tested:"
    , "\nUse cell cycle genes for clustering"
    , "\nRemove cell cycle genes from genes used for clustering"
    , "\nRegress out cell cycle scores"
    , "\n")
  )
# Save
ggsave(paste0(outGraph, "NumberRgIpCluster_Barplot.pdf")
  , width = 12, height = 18)


## Percent of cells passing expression filters per cluster

gg1 <- Percent_RgIpCluster_Barplot(exM = rmCcExM, seuratO = rmCcSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Remove CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Percent_RgIpCluster_Barplot(exM = rmCcExM, seuratO = rmCcSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Remove CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Percent_RgIpCluster_Barplot(exM = keepCcExM, seuratO = keepCcSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg4 <- Percent_RgIpCluster_Barplot(exM = keepCcExM, seuratO = keepCcSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg5 <- Percent_RgIpCluster_Barplot(exM = regCcExM, seuratO = regCcSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Regress CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg6 <- Percent_RgIpCluster_Barplot(exM = regCcExM, seuratO = regCcSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Regress CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)

# Plot grid
pg <- plot_grid(gg1, gg2, gg3, gg4, gg5, gg6, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of RG, oRG, vRG, IPC expression filters"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentRgIpCluster_Barplot.pdf")
  , width = 12, height = 18)


## nUMI per cluster box plots

df <- Format_Number_Cell_Types_Cluster_Dataframe(rmCcExM, rmCcSO)
df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))
gg1 <- ggplot(df, aes(x = CLUSTER, y = nUMI)) +
  geom_boxplot(aes(fill = CLUSTER)) +
  geom_jitter(aes(x = CLUSTER, y = nUMI)
    , size = 0.02, height = 0, alpha = 0.2) +
  theme(legend.position="none") +
  ggtitle("Remove CC genes")

df <- Format_Number_Cell_Types_Cluster_Dataframe(keepCcExM, keepCcSO)
df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))
gg2 <- ggplot(df, aes(x = CLUSTER, y = nUMI)) +
  geom_boxplot(aes(fill = CLUSTER)) +
  geom_jitter(aes(x = CLUSTER, y = nUMI)
    , size = 0.02, height = 0, alpha = 0.2) +
  theme(legend.position="none") +
  ggtitle("Keep CC genes")

df <- Format_Number_Cell_Types_Cluster_Dataframe(regCcExM, regCcSO)
df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))
gg3 <- ggplot(df, aes(x = CLUSTER, y = nUMI)) +
  geom_boxplot(aes(fill = CLUSTER)) +
  geom_jitter(aes(x = CLUSTER, y = nUMI)
    , size = 0.02, height = 0, alpha = 0.2) +
  theme(legend.position = "none") +
  ggtitle("Remove CC genes")

pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nnUMI per cluster"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
# Save
ggsave(paste0(outGraph, "nUMI_Cluster_Boxplot.png")
  , width = 15, height = 7)
################################################################################

## nUMI per expression filter box plots

# Remove CC genes
gg1 <- nUMI_RgIP_Boxplot(
  exM = rmCcExM, seuratO = rmCcSO, title = "Remove CC genes")
# Keep CC genes
gg2 <- nUMI_RgIP_Boxplot(
  exM = keepCcExM, seuratO = keepCcSO, title = "Keep CC genes")
# Regress CC
gg3 <- nUMI_RgIP_Boxplot(
  exM = regCcExM, seuratO = regCcSO, title = "Reg CC")
pg <- plot_grid(gg1, gg2, gg3, ncol = 3)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nnUMI per expression subset"
  , "\nCell cycle corrections tested:"
  , "\nUse cell cycle genes for clustering"
  , "\nRemove cell cycle genes from genes used for clustering"
  , "\nRegress out cell cycle scores"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
# Save
ggsave(paste0(outGraph, "nUMI_RgIP_Boxplot.png")
  , width = 13, height = 7)
################################################################################

### Number / Percent of cells in CC phase subset by markers

## Percent of cells in CC phase
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = keepCcExM, seuratO = keepCcSO)

df <- rbind(
  data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$IP > 0.5]), SUBSET = "vRG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$oRG > 0.5 & df1$IP > 0.5]), SUBSET = "oRG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "vRG+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$oRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "oRG+ Neuron+")
)
df$tableArg <- factor(df$tableArg, levels = c("G1", "S", "G2M"))
# Plot
ggplot(df, aes(x = SUBSET, y = Freq, fill = tableArg)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_discrete(name = "CC phase") + 
  xlab("Cell subset") +
  ylab("Percent of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nPercent of cell subsets in cell cycle phases"
    , "\nSubset by mean expression levels of groups of marker genes"
    , "\n(> 0.5 log normalized mean expression)")
  )
ggsave(paste0(outGraph, "PercentCCphase_RG_IP_Subset_Barplot.pdf")
  , width = 7, height = 7)


## Number of cells in CC phase

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = keepCcExM, seuratO = keepCcSO)

df <- rbind(
  data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$IP > 0.5]), SUBSET = "vRG+ IP+")
  , data.frame(table(df1$PHASE[df1$oRG > 0.5 & df1$IP > 0.5]), SUBSET = "oRG+ IP+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "vRG+ Neuron+")
  , data.frame(table(df1$PHASE[df1$oRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "oRG+ Neuron+")
)
df$Var1 <- factor(df$Var1, levels = c("G1", "S", "G2M"))
# Plot
ggplot(df, aes(x = SUBSET, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_discrete(name = "CC phase") + 
  xlab("Cell subset") +
  ylab("Number of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cell subsets in cell cycle phases"
    , "\nSubset by mean expression levels of groups of marker genes"
    , "\n(> 0.5 log normalized mean expression)")
  )
ggsave(paste0(outGraph, "NumberCCphase_RG_IP_Subset_Barplot.pdf")
  , width = 7, height = 7)
################################################################################

### Correlation of cluster EG and RG IP cells

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = keepCcExM, seuratO = keepCcSO)

ME_By_Markers_And_Phase <- function(markerGroup, phase)

cellType <- NA
cellType[df1$IP > 0.5 & df1$PHASE %in% "G2M"] <- "IPC"
eg1 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$RG > 0.5 & df1$PHASE == "G2M"] <- "RG"
eg2 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$IP > 0.5 & df1$PHASE %in% "S"] <- "IPC"
eg3 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$RG > 0.5 & df1$PHASE == "S"] <- "RG"
eg4 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$IP > 0.5 & df1$PHASE %in% "G1"] <- "IPC"
eg5 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$RG > 0.5 & df1$PHASE == "G1"] <- "RG"
eg6 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$IP > 0.5] <- "IPC"
eg7 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$RG > 0.5] <- "RG"
eg8 <- moduleEigengenes(keepCcExM, cellType)$eigengenes




cellType[df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE == "G2M"] <- "G2M"
cellType[df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE == "S"] <- "S"
cellType[df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE == "G1"] <- "G1"
eg10 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

cellType <- NA
cellType[df1$RG > 0.5 & df1$IP > 0.5] <- "All CC phase"
eg11 <- moduleEigengenes(keepCcExM, cellType)$eigengenes

df <- data.frame(eg1, eg2)




cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, egIPC$MEIPC)
cor(mnEx, egRG$MERG)

cellIDs <- row.names(df1)[df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, egIPC$MEIPC)
cor(mnEx, egRG$MERG)

cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, egIPC$MEIPC)
cor(mnEx, egRG$MERG)








cellType <- NA
cellType[keepCcSO@ident == 2] <- "IPC"
cellType[keepCcSO@ident %in% c(7,9)] <- "RG"
eg <- moduleEigengenes(keepCcExM, cellType)$eigengenes




cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)

cellIDs <- row.names(df1)[df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)

cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE == "G2M"]
exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)
cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)














df2 <- Positive_Negative_Expression_Flag(
  exDF = df1, highThreshold = 0.5, lowThreshold = 0.5)

cellIDs <- row.names(df2)[df2$TYPE == "IP+ RG+"]
cellIDs <- cellIDs[! is.na(cellIDs)]

exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)

cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)


cellIDs <- row.names(df2)[df2$TYPE == "IP+ vRG- oRG- RG-"]
cellIDs <- cellIDs[! is.na(cellIDs)]

exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)

cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)


cellIDs <- row.names(df2)[df2$TYPE == "oRG+ IP-"]
cellIDs <- cellIDs[! is.na(cellIDs)]

exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)

cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)

cellIDs <- row.names(df2)[df2$TYPE == "vRG+ IP-"]
cellIDs <- cellIDs[! is.na(cellIDs)]

exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)

cor(mnEx, eg$MEIPC)
cor(mnEx, eg$MERG)





cellIDs <- row.names(df2)[df2$TYPE == "oRG+ IP-"]
cellIDs <- cellIDs[! is.na(cellIDs)]

exM <- keepCcExM[ ,colnames(keepCcExM) %in% cellIDs]
mnEx <- rowMeans(exM)

cor(mnEx, eg$MEIPC)

mnIPC <- rowMeans(keepCcExM[ ,keepCcSO@ident == 2])
cor(mnEx, mnIPC)



# randomIDs <- sample(colnames(keepCcExM), size = length(cellIDs), replace = FALSE)
# exM <- keepCcExM[ ,colnames(keepCcExM) %in% randomIDs]
# mnEx <- rowMeans(exM)
# cor(mnEx, eg$MEIPC)


################################################################################


df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = keepCcExM, seuratO = keepCcSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 | df1$IP > 0.5 & df1$PHASE == "G2M"]
exM <- keepCcSO@scale.data[ ,colnames(keepCcSO@scale.data) %in% cellIDs]
mns <- rowMeans(exM)
genes <- names(sort(mns, decreasing = TRUE))[1:5000]
exM <- exM[row.names(exM) %in% genes, ]

png(paste0(outGraph, "pheatmap.png"))
pheatmap(exM, 
  cluster_row = TRUE,
  cluster_cols = TRUE
)
dev.off()
################################################################################