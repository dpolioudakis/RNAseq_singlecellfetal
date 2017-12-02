# Damon Polioudakis
# 2017-11-01
# Separate cycling vRG and oRG

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
require(RColorBrewer)
require(viridis)
require(fdrtool)
source("Function_Library.R")

options(stringsAsFactors = FALSE)

# Keep CC genes from variable gene list used for clustering
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

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
graphCodeTitle <- "Seurat_Cluster_Cycling_vRGoRG_DS2-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outData <- "../analysis/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"

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
  df$TYPE[df[ ,c("Neuron")] > highThreshold &
      df[ ,c("RG")] > highThreshold] <- "Neuron+ RG+"
  
  df$TYPE <- factor(df$TYPE, levels = c("IP+", "RG+ vRG- oRG- IP-", "vRG+ IP-"
    , "oRG+ IP-", "IP+ vRG- oRG- RG-", "IP+ RG+", "Endothelial+ IP+", 
    "Neuron+ IP+", "Neuron+ RG+"))
  df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))
  
  return(df)
}

Mean_Gene_Group_Expression <- function(exM, grouping) {
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% grouping]
  colMeans(
    exM[row.names(exM) %in% genes, ]
  )
}

Format_Number_Cell_Types_Cluster_Dataframe <- function (exM, seuratO) {
  
  mnExDF <- data.frame(
    vRG = Mean_Gene_Group_Expression(exM = exM, grouping = "vRG")
    , vRG_PollenS3 = Mean_Gene_Group_Expression(exM = exM, grouping = "vRG-PollenS3")
    , oRG = Mean_Gene_Group_Expression(exM = exM, grouping = "oRG")
    , oRG_PollenS3 = Mean_Gene_Group_Expression(exM = exM, grouping = "oRG-PollenS3")
    , RG = Mean_Gene_Group_Expression(exM = exM, grouping = "RG")
    , RG_PollenS3 = Mean_Gene_Group_Expression(exM = exM, grouping = "RG-PollenS3")
    , IP = Mean_Gene_Group_Expression(exM = exM, grouping = "IP")
    , Endothelial = Mean_Gene_Group_Expression(exM = exM, grouping = "Endothelial Cell")
    , Neuron = Mean_Gene_Group_Expression(exM = exM, grouping = "Neuron")
  )
  
  idx <- match(row.names(mnExDF), row.names(seuratO@meta.data))
  mnExDF$PHASE <- seuratO@meta.data$Phase[idx]
  mnExDF$CLUSTER <- seuratO@meta.data$res.0.6[idx]
  mnExDF$nUMI <- seuratO@meta.data$nUMI[idx]
  mnExDF$G2Mscore <- seuratO@meta.data$G2M.Score[idx]
  mnExDF$Sscore <- seuratO@meta.data$S.Score[idx]
  
  return(mnExDF)
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

MeanExprRank_Stdev_Variance_ScatterPlot <- function (exM, title) {
  mns <- rowMeans(exM)
  sdev <- apply(exM, 1, sd)
  variance <- apply(exM, 1, var)
  df <- data.frame(Mean = mns, StdDev = sdev)
  # df$CoefficentOfVariation = Coefficent_Of_Variation(df$Mean, df$StdDev)
  df$Variance <- variance
  df$MeanRank <- rank(df$Mean)
  gg1 <- ggplot(df, aes(x = MeanRank, y = Mean)) +
    geom_point(size = 0.1)
  gg2 <- ggplot(df, aes(x = MeanRank, y = StdDev)) +
    geom_point(size = 0.1)
  gg3 <- ggplot(df, aes(x = MeanRank, y = Variance)) +
    geom_point(size = 0.1)
  pg <- Plot_Grid(ggPlotsL = list(gg1, gg2, gg3), ncol = 3, rel_height = 0.2
    , title = title)
  return(pg)
}

Prcomp_Loadings_Plot <- function(pca, nGenes, nPCs, title) {
  
  # Example function call:
  # Prcomp_Loadings_Plot(pca = pca, nGenes = 1:20, nPCs = 1:8)
  
  # Plot highest loading genes
  ggL <- lapply(nPCs, function(pc) {
    df <- rbind(data.frame(PC = sort(pca$rotation[ ,pc])[nGenes])
      , data.frame(PC = sort(pca$rotation[ ,pc], decreasing = TRUE)[nGenes])
    )
    df$GENE <- factor(row.names(df), levels = row.names(df))
    gg <- ggplot(df, aes(x = PC, y = GENE)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      ggtitle(paste0("PC ", pc))
    return(gg)
  })
  Plot_Grid(ggPlotsL = ggL
    , ncol = 4
    , align = 'v'
    , axis = 'r'
    , rel_height = 0.1
    , title = title
  )
}

PCA_Format_For_GGplot <- function(pca) {
  
  df <- as.data.frame(pca$x)
  
  varExpL <- (pca$sdev)^2 / sum(pca$sdev^2)
  names(varExpL) <- paste0("PC", 1:length(varExpL))
  
  df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
    exM = noCentExM, seuratO = centSO)
  
  # Flag RG+, IP+, and RG+ IP+
  df$Cell_Subset_075 <- NA
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$Neuron > 0.75]] <- "Neuron"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$RG > 0.75]] <- "RG"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$IP > 0.75]] <- "IP"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$RG > 0.75 & df1$Neuron > 0.75]] <- "RG Neuron"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$RG > 0.75 & df1$IP > 0.75]] <- "RG IP"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$IP > 0.75 & df1$Neuron > 0.75]] <- "IP Neuron"
  df$Cell_Subset_075[rownames(df) %in% row.names(df1)[df1$RG > 0.75 & df1$IP > 0.75 & df1$Neuron > 0.75]] <- "RG IP Neuron"
  
  df$Cell_Subset <- NA
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$Neuron > 0.5]] <- "Neuron"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5]] <- "IP"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$Neuron > 0.5]] <- "RG Neuron"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5]] <- "RG IP"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5 & df1$Neuron > 0.5]] <- "IP Neuron"
  df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron > 0.5]] <- "RG IP Neuron"
  
  df$Cell_Subset_025 <- NA
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$Neuron > 0.25]] <- "Neuron"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$RG > 0.25]] <- "RG"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$IP > 0.25]] <- "IP"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$RG > 0.25 & df1$Neuron > 0.25]] <- "RG Neuron"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$RG > 0.25 & df1$IP > 0.25]] <- "RG IP"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$IP > 0.25 & df1$Neuron > 0.25]] <- "IP Neuron"
  df$Cell_Subset_025[rownames(df) %in% row.names(df1)[df1$RG > 0.25 & df1$IP > 0.25 & df1$Neuron > 0.25]] <- "RG IP Neuron"
  
  # Flag vRG+, oRG+, and vRG+ oRG+
  df$vRG_oRG_Subset <- NA
  df$vRG_oRG_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
  df$vRG_oRG_Subset[rownames(df) %in% row.names(df1)[df1$vRG > 0.5]] <- "vRG"
  df$vRG_oRG_Subset[rownames(df) %in% row.names(df1)[df1$oRG > 0.5]] <- "oRG"
  df$vRG_oRG_Subset[rownames(df) %in% row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]] <- "vRG oRG"
  
  # vRG and oRG marker expression
  idx <- match(rownames(df), row.names(df1))
  df$vRG <- df1$vRG[idx]
  df$oRG <- df1$oRG[idx]
  df$vRG_PollenS3 <- df1$vRG_PollenS3[idx]
  df$oRG_PollenS3 <- df1$oRG_PollenS3[idx]
  
  # RG IP Neuron marker expression
  idx <- match(rownames(df), row.names(df1))
  df$RG <- df1$RG[idx]
  df$RG_PollenS3 <- df1$RG_PollenS3[idx]
  df$IP <- df1$IP[idx]
  df$Neuron <- df1$Neuron[idx]
  
  # Seurat S phase and G2M scores
  idx <- match(rownames(df), row.names(centSO@meta.data))
  df$G2M_Score <- centSO@meta.data$G2M.Score[idx]
  df$S_Score <- centSO@meta.data$S.Score[idx]
  
  pcaL <- list(PCA_for_ggplot = df, Variance_Explained = varExpL)
  
  return(pcaL)
}

PCA_Plot <- function(ggDF, varExpL, PCx, PCy){
  gg <- ggplot(ggDF, aes(x = ggDF[[PCx]], y = ggDF[[PCy]], color = colorBy)) +
    geom_point() +
    xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
    ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
  return(gg)
}

PCA_Plot_PC1to8 <- function(pcaL, colorBy, limLow = NULL, limHigh = NULL) {
  
  # PCA data frame for ggplot
  ggDF <- pcaL[["PCA_for_ggplot"]]
  
  varExpL <- pcaL[["Variance_Explained"]]
  
  # Color by
  ggDF$colorBy = ggDF[ ,colnames(ggDF) %in% colorBy]
  # Set expression limits
  if (class(ggDF$colorBy) == "numeric" & ! is.null(limLow) & ! is.null(limHigh)) {
    ggDF$colorBy[ggDF$colorBy < limLow] <- limLow
    ggDF$colorBy[ggDF$colorBy > limHigh] <- limHigh
  }
  
  # Plot
  ggL <- list(
    PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC2")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC3")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC4")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC2", PCy = "PC3")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC2", PCy = "PC4")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC3", PCy = "PC4")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC5")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC6")
    , PCA_Plot(ggDF = ggDF, varExpL = varExpL, PCx = "PC1", PCy = "PC7")
    )
  # Use continuous color scale if value to color by is continuous
  if (class(ggDF$colorBy) == "numeric") {
    # ggL <- lapply(ggL, function(gg) {gg + scale_color_viridis()})
    ggL <- lapply(ggL, function(gg) {
      gg <- gg +
        scale_color_distiller(name = "Normalized\nexpression"
          , type = "div", palette = 5, direction = -1, limits = c(limLow, limHigh)) +
        geom_point(size = 0.05)
      return(gg)})
  }
  return(ggL)
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

### Percent of cells passing expression filters per cluster

gg1 <- Percent_RgIpCluster_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Percent_RgIpCluster_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)

# Plot grid
pg <- plot_grid(gg1, gg2, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of RG, oRG, vRG, IPC expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentRgIpCluster_Barplot.pdf")
  , width = 12, height = 7)
################################################################################

### Number / Percent of cells in CC phase subset by markers
print("### Number / Percent of cells in CC phase subset by markers")

## Percent of cells in CC phase
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)

df <- rbind(
  data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$oRG > 0.5]), SUBSET = "vRG+ oRG+")
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
  exM = noCentExM, seuratO = centSO)

df <- rbind(
  data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$oRG > 0.5]), SUBSET = "vRG+ oRG+")
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

### Correlation of cluster mean expression profile and RG IP cells
print("### Correlation of cluster EG and RG IP cells")

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)

cellSubsetsL <- list(
  IP_G2M = df1$IP > 0.5 & df1$PHASE %in% "G2M"
  , RG_G2M = df1$RG > 0.5 & df1$PHASE %in% "G2M"
  
  , IP_S = df1$IP > 0.5 & df1$PHASE %in% "S"
  , RG_S = df1$RG > 0.5 & df1$PHASE %in% "S"
  
  , IP_G1 = df1$IP > 0.5 & df1$PHASE %in% "G1"
  , RG_G1 = df1$RG > 0.5 & df1$PHASE %in% "G1"
  
  , IP = df1$IP > 0.5
  , RG = df1$RG > 0.5
  
  , RG_IPn_G2M = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "G2M"
  , RG_IP_G2M = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G2M"
  , RGn_IP_G2M = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G2M"
  
  , RG_IPn_S = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "S"
  , RG_IP_S = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "S"
  , RGn_IP_S = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "S"
  
  , RG_IPn_G1 = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "G1"
  , RG_IP_G1 = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G1"
  , RGn_IP_G1 = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G1"
  
  , RG_IP = df1$RG > 0.5 & df1$IP > 0.5
  , RG_IPn = df1$RG > 0.5 & df1$IP < 0.5
  , RGn_IP = df1$RG < 0.5 & df1$IP > 0.5
)
egL <- lapply(names(cellSubsetsL), function(name) {
  cellSubset <- cellSubsetsL[[name]]
  exM <- noCentExM
  mns <- rowMeans(exM)
  genes <- names(sort(mns, decreasing = TRUE))[1:2500]
  exM <- exM[row.names(exM) %in% genes, cellSubset]
  rowMeans(exM)
  # print(str(exM))
  # colors <- rep("TRUE", nrow(exM))
  # eg <- moduleEigengenes(t(exM), colors)$eigengenes
  # eg <- eg["METRUE"]
  # colnames(eg) <- name
  # return(eg)
})
names(egL) <- names(cellSubsetsL)
corM <- round(cor(data.frame(egL)), 2)

write.csv(corM
  , file = paste0(outTable, "RG_IP_subsets_Mean_correlation_matrix.csv")
  , quote = FALSE)

# randomIDs <- sample(colnames(noCentExM), size = length(cellIDs), replace = FALSE)
# exM <- noCentExM[ ,colnames(noCentExM) %in% randomIDs]
# mnEx <- rowMeans(exM)
# cor(mnEx, eg$MEIPC)
################################################################################

### vRG and oRG gene expression by cluster
print("### vRG and oRG gene expression by cluster")

ggL <- lapply(c("vRG", "oRG", "RG"), function(grouping) {
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% grouping]
  genes <- genes[! duplicated(genes)]
  gg <- Gene_Expression_Facet_By_Cluster_ViolinPlot(
    genes = genes
    , exprM = noCentExM[ ,colnames(noCentExM) %in% names(centSO@ident)[centSO@ident %in% c(7,8,9,10)]]
    , clusterIDs = centSO@ident[centSO@ident %in% c(7,8,9,10)]
    , geneOrder = genes
    , ncol = 4
    , ggtitle = grouping)
  gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(gg)
})
Plot_Grid(ggPlotsL = ggL
  , ncol = 1
  , title = paste0(graphCodeTitle
    , "\n\nvRG, oRG, and RG markers expression by cluster")
  , rel_height = 0.1)
ggsave(paste0(outGraph, "vRG_oRG_RG_Expression_ViolinPlot.png")
  , width = 13, height = 11)
################################################################################

### Heatmap and heirarchical clustering of RG and IP cells
print("### Heatmap and heirarchical clustering of RG and IP cells")

## Heatmap and heirarchical clustering of RG and IP G2M cells

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.5 | df1$IP > 0.5 & df1$Neuron < 0.5]
cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "G2M"])
exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
mns <- rowMeans(exM)
genes <- names(sort(mns, decreasing = TRUE))[1:2500]
exM <- exM[row.names(exM) %in% genes, ]

annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
row.names(annotation_col) <- row.names(df1)
annotation_col$Type[df1$RG > 0.5] <- "RG+"
annotation_col$Type[df1$IP > 0.5] <- "IP+"
annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

annotation_col$vRGoRG <- NA
ids <- row.names(df1)[df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
ids <- row.names(df1)[df1$vRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"

idx <- match(row.names(annotation_col), row.names(df1))
annotation_col$G2M <- df1$G2Mscore[idx]
annotation_col$Sscore <- df1$Sscore[idx]

exM[exM > 3] <- 3
exM[exM < -3] <- 3
breaks <- seq(-2, 2, by = 0.1)

png(paste0(outGraph, "RG_IP_G2M_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
pheatmap(exM, 
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
)
dev.off()


## Heatmap and heirarchical clustering of RG cluster 10 cells

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
df1 <- df1[row.names(df1) %in% cellIDs, ]
exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
mns <- rowMeans(exM)
genes <- names(sort(mns, decreasing = TRUE))[1:2500]
exM <- exM[row.names(exM) %in% genes, ]

annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
row.names(annotation_col) <- row.names(df1)
annotation_col$Type[df1$RG > 0.5] <- "RG+"
annotation_col$Type[df1$IP > 0.5] <- "IP+"
annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

annotation_col$vRGoRG <- NA
ids <- row.names(df1)[df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
ids <- row.names(df1)[df1$vRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"

idx <- match(row.names(annotation_col), row.names(df1))
annotation_col$G2M <- df1$G2Mscore[idx]
annotation_col$Sscore <- df1$Sscore[idx]

exM[exM > 3] <- 3
exM[exM < -3] <- 3
breaks <- seq(-3, 3, by = 0.1)

png(paste0(outGraph, "RG_Cluster10_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
pheatmap(exM, 
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
  , treeheight_col = 150
  , cutree_cols = 10
)
dev.off()


## Heatmap and heirarchical clustering of RG cluster 8 cells

df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
df1 <- df1[row.names(df1) %in% cellIDs, ]
exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
mns <- rowMeans(exM)
genes <- names(sort(mns, decreasing = TRUE))[1:2500]
exM <- exM[row.names(exM) %in% genes, ]

annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
row.names(annotation_col) <- row.names(df1)
annotation_col$Type[df1$RG > 0.5] <- "RG+"
annotation_col$Type[df1$IP > 0.5] <- "IP+"
annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

annotation_col$vRGoRG <- NA
ids <- row.names(df1)[df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
ids <- row.names(df1)[df1$vRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"

idx <- match(row.names(annotation_col), row.names(df1))
annotation_col$G2M <- df1$G2Mscore[idx]
annotation_col$Sscore <- df1$Sscore[idx]

exM[exM > 3] <- 3
exM[exM < -3] <- 3
breaks <- seq(-3, 3, by = 0.1)

png(paste0(outGraph, "RG_Cluster8_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
pheatmap(exM, 
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
  , treeheight_col = 150
  , cutree_cols = 10
)
dev.off()
################################################################################

### PCA of RG+ and / or IP+ cells
print("### PCA of RG+ and / or IP+ cells")

exLM <- list()

# Identify cluster 8 cells
cellIDs <- names(centSO@ident[centSO@ident == 8])
exLM[["Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify cluster 10 cells
cellIDs <- names(centSO@ident[centSO@ident == 10])
exLM[["Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or IP+ but Neuron- cluster 8 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IP_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or IP+ but Neuron- cluster 10 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IP_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 8 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuron_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 10 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuron_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ IP- Neuron- cluster 8 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+  IP- Neuron- cluster 10 cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# # Identify RG+ cluster 8 cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
# exLM[["RG_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# 
# # Identify RG+ cluster 10 cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
# exLM[["RG_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# 
# # Identify RG+ and / or IP+ G2M cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 | df1$IP > 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "G2M"])
# exLM[["RG_IP_G2M"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# 
# # Identify RG+ and / or IP+ S phase cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 | df1$IP > 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "S"])
# exLM[["RG_IP_Sphase"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# 
# # Identify RG+ G2M cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "G2M"])
# exLM[["RG_G2M"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# 
# # Identify RG+ S phase cells
# df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "S"])
# exLM[["RG_Sphase"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Plot mean expression rank, standard deviation, and variance
pgL <- lapply(names(exLM), function(name) {
  exM <- exLM[[name]]
  MeanExprRank_Stdev_Variance_ScatterPlot(
    exM = exM
    , title = name
  )
})
# Combine plots from each cell subset
Plot_Grid(ggPlotsL = pgL, ncol = 1, rel_height = 0.1
  , title = paste0(
    graphCodeTitle
    , "\n\nSubset cells by RG+ and / or IP+ expression and G2M or S phase"))
ggsave(paste0(outGraph, "RG_IP_MeanExpr_StdDev.png"), width = 13, height = 20)

# Subset to top 2500 expressed genes
exLM <- lapply(exLM, function(exM){
  mns <- rowMeans(exM)
  genes <- names(sort(mns, decreasing = TRUE))[1:2500]
  exM <- exM[row.names(exM) %in% genes, ]
  return(exM)
})

# PCA
pcaL <- lapply(exLM, function(exM){
  pca <- prcomp(t(exM))
  return(pca)
})

# Plot PCA

# Plot PCA loadings
lapply(names(pcaL), function(name){
  pca <- pcaL[[name]]
  Prcomp_Loadings_Plot(pca = pca, nGenes = c(1:20), nPCs = c(1:8)
    , title = paste0(graphCodeTitle, "\n\n", name, "\nGenes with highest PC loadings"))
  ggsave(paste0(outGraph, name, "_PCAloadings.pdf"), width = 13
    , height = 20, limitsize = FALSE)
})

# Format for ggplot
ldf <- lapply(pcaL, function(pca){
  PCA_Format_For_GGplot(pca)
})

# mns <- rowMeans(noCentExM[ ,centSO@ident == 7])
# corCluster7 <- apply(noCentExM, 2, function(col){cor(col, mns)})
# 
# mns <- rowMeans(noCentExM[ ,centSO@ident == 9])
# corCluster9 <- apply(noCentExM, 2, function(col){cor(col, mns)})
# 
# ldf <- lapply(ldf, function(df){
#   
#   idx <- match(row.names(df), names(corCluster7))
#   df$Cor_Cluster7 <- corCluster7[idx]
#   
#   idx <- match(row.names(df), names(corCluster9))
#   df$Cor_Cluster9 <- corCluster7[idx]
#   
#   return(df)
# })

# PCA plots
lapply(names(ldf), function(name){
  
  pcaL <- ldf[[name]]
  
  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_MarkerLabel05.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "Cell_Subset_025")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.25 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_MarkerLabel025.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "Cell_Subset_075")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.75 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_MarkerLabel075.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "vRG_oRG_Subset")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG+ and / or oRG+"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_vRGoRG.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "S_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by S phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_Sscore.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "G2M_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by G2M phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_G2Mscore.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "vRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_vRG.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "vRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_vRGPollenS3.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "oRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_oRG.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "oRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_oRGPollenS3.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "RG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_RG.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "RG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_RGPollenS3.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "IP", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_IP.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_Neuron.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "IP", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_IP025scale.png"), width = 14, height = 14)
  
  ggL <- PCA_Plot_PC1to8(pcaDF = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraph, name, "_PCA_Neuron025scale.png"), width = 14, height = 14)
  
  # # PCA plots
  # lapply(names(ldf), function(name){
  #   
  # df <- ldf[[name]]
  
  # ggL <- PCA_Plot_PC1to8(pcaDF = df, colorBy = "Cor_Cluster7", limLow = -1, limHigh = 1)
  # Plot_Grid(ggPlotsL = ggL
  #   , ncol = 3
  #   , rel_height = 0.1
  #   , title = paste0(graphCodeTitle
  #     , "\n\n", name, " PCA"
  #     , "\nColored by correlation to Cluster 7 (oRG) mean expression profile"
  #     , "\n+ = > 0.5 log normalized expression")
  # )
  # ggsave(paste0(outGraph, name, "_PCA_CorCluster7.png"), width = 14, height = 14)
  # 
  # ggL <- PCA_Plot_PC1to8(pcaDF = df, colorBy = "Cor_Cluster9", limLow = -1, limHigh = 1)
  # Plot_Grid(ggPlotsL = ggL
  #   , ncol = 3
  #   , rel_height = 0.1
  #   , title = paste0(graphCodeTitle
  #     , "\n\n", name, " PCA"
  #     , "\nColored by correlation to Cluster 9 (vRG) mean expression profile"
  #     , "\n+ = > 0.5 log normalized expression")
  # )
  # ggsave(paste0(outGraph, name, "_PCA_CorCluster9.png"), width = 14, height = 14)
  
})
################################################################################

### Expression of vRG, oRG, RG markers versus cell cycle score
print("### Expression of vRG, oRG, RG markers by cell cycle phase")

# Violin plots of Seurat cell cycle scores
df1 <- centSO@meta.data[c("S.Score", "G2M.Score", "res.0.6")]
df1$res.0.6 <- factor(df1$res.0.6, levels = sort(as.numeric(unique(df1$res.0.6))))
ggplot(df1, aes(x = res.0.6, y = S.Score)) +
  geom_violin() +
  geom_jitter(size = 0.01, alpha = 0.25) +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nS phase score by cluster"))
ggsave(paste0(outGraph, "CellCycleSphaseScore_By_Cluster_Violin.png")
  , width = 7, height = 6)
ggplot(df1, aes(x = res.0.6, y = G2M.Score)) +
  geom_violin() +
  geom_jitter(size = 0.01, alpha = 0.25) +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nG2M phase score by cluster"))
ggsave(paste0(outGraph, "CellCycleG2MScore_By_Cluster_Violin.png")
  , width = 7, height = 6)

df <- data.frame(Phase = rep("G0", ncol(noCentExM)))
df$Phase[centSO@meta.data["S.Score"] > 0.2] <- "S phase"
df$Phase[centSO@meta.data["G2M.Score"] > 0.5] <- "G2/M"

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "RG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$RG <- mnEx

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "vRG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$vRG <- mnEx

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "oRG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$oRG <- mnEx

# Subset to RG+ cells
df <- df[df$RG > 0.5, ]

# Format for ggplot
df <- melt(df)
df$Phase <- factor(df$Phase, levels = c("G0", "S phase", "G2/M"))

# Plot
ggplot(df, aes(x = variable, y = value, fill = Phase)) +
  geom_boxplot() +
  xlab("Marker gene group") +
  ylab("Mean normalized expression") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMean expression of vRG, oRG, or RG markers in RG+ cells by cell cycle"
    , "\nSubset to cells to > 0.5 mean normalized expression of RG markers")
  )
ggsave(paste0(outGraph, "RGexpr_By_CellCycle.png"), width = 7, height = 6)

# Plot
ggplot(df, aes(x = variable, y = value, fill = Phase)) +
  geom_violin() +
  geom_jitter(size = 0.25, alpha = 0.25) +
  xlab("Marker gene group") +
  ylab("Mean normalized expression") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMean expression of vRG, oRG, or RG markers in RG+ cells by cell cycle"
    , "\nSubset to cells to > 0.5 mean normalized expression of RG markers")
  )
ggsave(paste0(outGraph, "RGexpr_By_CellCycle_Violin.png"), width = 7, height = 6)
################################################################################

### DE RG vs IP vs IP+ RG+

## DE RG vs IP clusters

# Subset expression matrix to RG and IP clusters
clusters1 <- 2
clusters2 <- c(7,9)
ids <- names(centSO@ident)[centSO@ident %in% c(clusters1, clusters2)]
exM <- as.matrix(centSO@data)
exM <- exM[ ,colnames(exM) %in% ids]
# DE Linear model
termsDF <- centSO@meta.data[
  row.names(centSO@meta.data) %in% ids
  , c("nUMI", "librarylab", "individual", "res.0.6")]
# Add term TRUE/FALSE cell is in cluster
termsDF$clusters <- "clusters1"
termsDF$clusters[termsDF$res.0.6 %in% clusters2] <- "clusters2"
deLM <- DE_Linear_Model(
  exDatDF = exM
  , termsDF = termsDF
  , mod = "y ~ clusters+nUMI+librarylab+individual")
# Format LM DE
deDF <- data.frame(Log_FC = deLM$coefmat[ ,"clustersclusters2"]
  , Pvalue = deLM$pvalmat[ ,"clustersclusters2"])
deDF <- deDF[order(deDF$Log_FC), ]
deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1

# FDR correct
# NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
corrected <- fdrtool(deDF$Pvalue, statistic = "pvalue", plot = FALSE)
deDF$FDR <- corrected$lfdr
# Check
table(deDF$Pvalue < 0.05)
table(deDF$FDR < 0.05)
print(head(deDF))

# Save as csv
write.csv(deDF, file = paste0(outTable, "DE_RGvsIP.csv")
  , quote = FALSE)
deDF <- read.csv(paste0(outTable, "DE_RGvsIP.csv")
  , header = TRUE, row.names = 1)


## Heatmap of RG vs IP DE genes in IP+ and Neuron- or RG+ and Neuron- cells

# S phase cells

# Subset expression matrix to DE genes and IP+ and Neuron- or RG+ and Neuron- cells
# DE genes
genes <- row.names(deDF)[deDF$FDR < 0.05]
# IP+ and Neuron- or RG+ and Neuron- cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
ids <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$Neuron < 0.25
  | df1$CLUSTER == 8 & df1$IP > 0.5 & df1$Neuron < 0.25]
# Subset
exM <- as.matrix(centSO@scale.data)
exM <- exM[row.names(exM) %in% genes, colnames(exM) %in% ids]

# Heatmap column color bars
df2 <- df1[row.names(df1) %in% ids, ]
annotation_col <- data.frame(Type = rep(NA, nrow(df2)))
row.names(annotation_col) <- row.names(df2)
annotation_col$Type[df2$RG > 0.5] <- "RG+"
annotation_col$Type[df2$IP > 0.5] <- "IP+"
annotation_col$Type[df2$RG > 0.5 & df2$IP > 0.5] <- "RG+ IP+"
annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

# Heatmap row color bars
annotation_row <- data.frame(DE = rep(NA, nrow(exM)))
row.names(annotation_row) <- row.names(exM)
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC > 0]] <- "RG"
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC < -0]] <- "IP"

breaks <- seq(-2, 2, by = 0.1)

png(paste0(outGraph, "DE_RG_IP_S_pheatmap.png")
  , width = 9, height = 9, units = "in", res = 300)
pheatmap(exM, 
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , annotation_row = annotation_row
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
)
dev.off()

# G2M phase cells

# Subset expression matrix to DE genes and IP+ and Neuron- or RG+ and Neuron- cells
# DE genes
genes <- row.names(deDF)[deDF$FDR < 0.05]
# IP+ and Neuron- or RG+ and Neuron- cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)
ids <- row.names(df1)[df1$CLUSTER == 10 & df1$RG > 0.5 & df1$Neuron < 0.25
  | df1$CLUSTER == 10 & df1$IP > 0.5 & df1$Neuron < 0.25]
# Subset
exM <- as.matrix(centSO@scale.data)
exM <- exM[row.names(exM) %in% genes, colnames(exM) %in% ids]

# Heatmap column color bars
df2 <- df1[row.names(df1) %in% ids, ]
annotation_col <- data.frame(Type = rep(NA, nrow(df2)))
row.names(annotation_col) <- row.names(df2)
annotation_col$Type[df2$RG > 0.5] <- "RG+"
annotation_col$Type[df2$IP > 0.5] <- "IP+"
annotation_col$Type[df2$RG > 0.5 & df2$IP > 0.5] <- "RG+ IP+"
annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

# Heatmap row color bars
annotation_row <- data.frame(DE = rep(NA, nrow(exM)))
row.names(annotation_row) <- row.names(exM)
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC > 0]] <- "RG"
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC < -0]] <- "IP"

breaks <- seq(-2, 2, by = 0.1)

png(paste0(outGraph, "DE_RG_IP_G2M_pheatmap.png")
  , width = 9, height = 9, units = "in", res = 300)
pheatmap(exM, 
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , annotation_row = annotation_row
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
)
dev.off()


## DE RG+ vs RG+ IP+ S phase

DE_CellGroup1_Vs_CellGroup2 <- function(ids1, ids2) {
  
  # Subset expression matrix
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% c(ids1, ids2)]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% c(ids1, ids2)
    , c("nUMI", "librarylab", "individual", "res.0.6")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$groups <- "ids1"
  termsDF$groups[row.names(termsDF) %in% ids2] <- "ids2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ groups+nUMI+librarylab+individual")
  
  print(head(deLM$coefmat))
  
  # Format LM DE
  deDF <- data.frame(Log_FC = deLM$coefmat[ ,2]
    , Pvalue = deLM$pvalmat[ ,2])
  deDF <- deDF[order(deDF$Log_FC), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  
  # FDR correct
  # NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
  corrected <- fdrtool(deDF$Pvalue, statistic = "pvalue", plot = FALSE)
  deDF$FDR <- corrected$lfdr
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))
  
  return(deDF)
}

# Mean expression of marker genes to use to subset groups of cells
df1 <- Format_Number_Cell_Types_Cluster_Dataframe(
  exM = noCentExM, seuratO = centSO)

# S phase RG+ Neuron- vs RG+ IP+ Neuron-
ids1 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron < 0.25]
ids2 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
# DE
deRGvsRgIpDF <- DE_CellGroup1_Vs_CellGroup2(ids1, ids2)
# Save as csv
write.csv(deRGvsRgIpDF, file = paste0(outTable, "DE_RGvsRGIP_S.csv")
  , quote = FALSE)
deRGvsRgIpDF <- read.csv(paste0(outTable, "DE_RGvsRGIP_S.csv")
  , header = TRUE, row.names = 1)


# S phase RG+ Neuron- vs RG+ IP+ Neuron-
ids1 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron < 0.25]
ids2 <- row.names(df1)[df1$CLUSTER == 8 & df1$IP > 0.5 & df1$RG < 0.25 & df1$Neuron < 0.25]
# DE
deIPvsRgIpDF <- DE_CellGroup1_Vs_CellGroup2(ids1, ids2)
# Save as csv
write.csv(deIPvsRgIpDF, file = paste0(outTable, "DE_IPvsRGIP_S.csv")
  , quote = FALSE)
deIPvsRgIpDF <- read.csv(paste0(outTable, "DE_IPvsRGIP_S.csv")
  , header = TRUE, row.names = 1)


intersect(
  row.names(deDF)[deDF$Log_FC > 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC > 0 & deDF$FDR < 0.05]
)

intersect(
  row.names(deDF)[deDF$Log_FC < 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC < 0 & deDF$FDR < 0.05]
)

intersect(
  row.names(deDF)[deDF$Log_FC > 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC < 0 & deDF$FDR < 0.05]
)



intersect(names(deLM1$coefmat[ ,"clustersclusters2"])[deLM1$coefmat[ ,"clustersclusters2"] < 0.5]
, names(deLM$coefmat[,"groupsids2"])[deLM$coefmat[,"groupsids2"] > 0.5])
################################################################################

### S phase gene expression

# Sort genes by mean expression
mnEx <- rowMeans(noCentExM[row.names(noCentExM) %in% sGenes, ])
mnEx <- sort(mnEx, decreasing = TRUE)
genes <- names(mnEx)

# Number of cells each marker expressed in
m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
df1 <- melt(m1)
df1$Cluster <- centSO@ident[match(df1$Var2, names(centSO@ident))]
df1 <- aggregate(value~Var1+Cluster, df1, sum)
df1$Var1 <- factor(df1$Var1, levels = genes)

# Plot number of cells each marker expressed in
ggplot(df1, aes(x = Cluster, y = value)) +
  facet_wrap(~Var1, scales ="free_x", ncol = 3) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cells expressing genes"
    , "\n(> 0.5 normalized expression)"))
ggsave(paste0(outGraph, "Sgenes_Number_Barplot.png")
  , width = 12, height = 34)

# Feature plot
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1, limHigh = 3)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.05, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of S phase genes sorted by mean expression"
    , "\nNormalized expression"
    , "\nGene list from Tirosh et al. 2016"))
ggsave(paste0(outGraph, "Sgenes_FeaturePlot.png"), width = 20, height = 70
  , limitsize = FALSE)
################################################################################

### Intersection of S phase genes, RG markers, IP markers for FISH probes

## tSNE colored by intersection and heatmap of numbers of intersections

# Genes to intersect
genes <- c("PCNA", "EOMES", "SOX2", "PAX6", "HOPX", "CRYAB")

# Number of cells each marker expressed in
m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
df1 <- melt(m1)
df1$Cluster <- centSO@ident[match(df1$Var2, names(centSO@ident))]
df1 <- aggregate(value~Var1+Cluster, df1, sum)

# Plot number of cells each marker expressed in
ggplot(df1, aes(x = Cluster, y = value)) +
  facet_wrap(~Var1, scales ="free_x") +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cells expressing genes"
    , "\n(> 0.5 normalized expression)"))
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Number_Barplot.png")
  , width = 12, height = 10)

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.05, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"))
)
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Intersection_tSNE.png")
  , width = 20, height = 38)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\n(> 0.5 normalized expression)")
)
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Intersection_Heatmap.png")
  , width = 7, height = 7)
################################################################################


