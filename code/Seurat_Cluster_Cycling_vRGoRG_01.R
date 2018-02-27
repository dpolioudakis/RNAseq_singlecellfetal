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
# require(pheatmap)
require(RColorBrewer)
require(viridis)
# require(fdrtool)
require(monocle)
require(WGCNA)
require(gridExtra)
source("Function_Library.R")

options(stringsAsFactors = FALSE)

# Keep CC genes from variable gene list used for clustering
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TESTcluster0278910_seuratO.Robj")
# centSO <- ssCentSO; rm(ssCentSO)
# noCentExM <- ssNoCentExM; rm(ssNoCentExM)

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccGenes <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
sGenes <- ccGenes[1:43]
g2mGenes <- ccGenes[44:98]

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_Cluster_Cycling_vRGoRG_DS2-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

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
  # (from Average_MarkersExp_Per_Cell)
  # vRG         oRG          RG          IP Endothelial
  # ACCAGCTAGCCT -0.12961098 -0.04274303 -0.16390677 -0.05747777  0.08723860
  # GGAAGGACTGCA  0.17994922 -0.08862593  0.47211654 -0.02872334  1.65468242
  # Neuron PHASE CLUSTER nUMI
  # ACCAGCTAGCCT  0.8625643    G1       3 6165
  # GGAAGGACTGCA -0.9443603   G2M      11 4166

  df <- exDF
  df$TYPE <- NA
  df$TYPE[df[ ,c("Neuron")] > highThreshold] <- "Neuron+"
  df$TYPE[df[ ,c("IP")] > highThreshold] <- "IP+"
  df$TYPE[df[ ,c("RG")] > highThreshold] <- "RG+"
  df$TYPE[df[ ,c("RG")] > highThreshold &
      apply((df[ ,c("vRG", "oRG")] < lowThreshold), 1, all)] <- "RG+ vRG- oRG-"
  # df$TYPE[df[ ,c("vRG")] > highThreshold &
  #     df[ ,c("IP")] < lowThreshold] <- "vRG+ IP-"
  # df$TYPE[df[ ,c("oRG")] > highThreshold &
  #     df[ ,c("IP")] < lowThreshold] <- "oRG+ IP-"
  # df$TYPE[df[ ,c("IP")] > highThreshold &
  #     apply((df[ ,c("vRG", "oRG", "RG")] < lowThreshold), 1, all)] <- "IP+ vRG- oRG- RG-"
  df$TYPE[apply((df[ ,c("IP", "RG")] > highThreshold), 1, all)] <- "IP+ RG+"
  df$TYPE[apply((df[ ,c("Neuron", "IP")] > highThreshold), 1, all)] <- "Neuron+ IP+"
  df$TYPE[apply((df[ ,c("Neuron", "RG")] > highThreshold), 1, all)] <- "Neuron+ RG+"

  # Controls
  df$TYPE[df[ ,c("Endothelial")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Endothelial+ IP+"
  df$TYPE[df[ ,c("Interneuron")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Interneuron+ IP+"

  df$TYPE <- factor(df$TYPE, levels = c("Neuron+", "IP+", "RG+"
    , "RG+ vRG- oRG-", "IP+ RG+", "Neuron+ IP+", "Neuron+ RG+"
    , "Endothelial+ IP+", "Interneuron+ IP+"))
  df$CLUSTER <- factor(df$CLUSTER, levels = sort(unique(as.numeric(df$CLUSTER))))

  return(df)
}

Mean_Gene_Group_Expression <- function(exM, grouping) {
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% grouping]
  colMeans(
    exM[row.names(exM) %in% genes, ]
  )
}

Average_MarkersExp_Per_Cell <- function(exM, seuratO) {

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
    , Interneuron = Mean_Gene_Group_Expression(exM = exM
      , grouping = "GABAergic interneuron ")
  )

  idx <- match(row.names(mnExDF), row.names(seuratO@meta.data))
  mnExDF$PHASE <- seuratO@meta.data$Phase[idx]
  mnExDF$CLUSTER <- seuratO@meta.data$res.0.6[idx]
  mnExDF$nUMI <- seuratO@meta.data$nUMI[idx]
  mnExDF$G2Mscore <- seuratO@meta.data$G2M.Score[idx]
  mnExDF$Sscore <- seuratO@meta.data$S.Score[idx]

  return(mnExDF)
}

Percent_Of_Table <- function(tableArg) {
  tb <- table(tableArg, exclude = NULL) / sum(table(tableArg, exclude = NULL)) * 100
  v1 <- data.frame(tb)$Freq
  return(tb)
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

  df1 <- Average_MarkersExp_Per_Cell(
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
    # xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
    # ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
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
        geom_point(size = 0.01)
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

Marker_Expression_Flag <- function(df, df1){

  # Flag RG+, IP+, and RG+ IP+
  df$Cell_Subset_075 <- NA
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$Neuron > 0.75]] <- "Neuron"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$RG > 0.75]] <- "RG"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$IP > 0.75]] <- "IP"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$RG > 0.75 & df1$Neuron > 0.75]] <- "RG Neuron"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$RG > 0.75 & df1$IP > 0.75]] <- "RG IP"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$IP > 0.75 & df1$Neuron > 0.75]] <- "IP Neuron"
  df$Cell_Subset_075[df$CellID %in% row.names(df1)[df1$RG > 0.75 & df1$IP > 0.75 & df1$Neuron > 0.75]] <- "RG IP Neuron"

  df$Cell_Subset <- NA
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$Neuron > 0.5]] <- "Neuron"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$IP > 0.5]] <- "IP"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$RG > 0.5 & df1$Neuron > 0.5]] <- "RG Neuron"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5]] <- "RG IP"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$IP > 0.5 & df1$Neuron > 0.5]] <- "IP Neuron"
  df$Cell_Subset[df$CellID %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron > 0.5]] <- "RG IP Neuron"

  df$Cell_Subset_025 <- NA
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$Neuron > 0.25]] <- "Neuron"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$RG > 0.25]] <- "RG"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$IP > 0.25]] <- "IP"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$RG > 0.25 & df1$Neuron > 0.25]] <- "RG Neuron"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$RG > 0.25 & df1$IP > 0.25]] <- "RG IP"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$IP > 0.25 & df1$Neuron > 0.25]] <- "IP Neuron"
  df$Cell_Subset_025[df$CellID %in% row.names(df1)[df1$RG > 0.25 & df1$IP > 0.25 & df1$Neuron > 0.25]] <- "RG IP Neuron"

  # Flag vRG+, oRG+, and vRG+ oRG+
  df$vRG_oRG_Subset <- NA
  df$vRG_oRG_Subset[df$CellID %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
  df$vRG_oRG_Subset[df$CellID %in% row.names(df1)[df1$vRG > 0.5]] <- "vRG"
  df$vRG_oRG_Subset[df$CellID %in% row.names(df1)[df1$oRG > 0.5]] <- "oRG"
  df$vRG_oRG_Subset[df$CellID %in% row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]] <- "vRG oRG"

  # vRG and oRG marker expression
  idx <- match(df$CellID, row.names(df1))
  df$vRG <- df1$vRG[idx]
  df$oRG <- df1$oRG[idx]
  df$vRG_PollenS3 <- df1$vRG_PollenS3[idx]
  df$oRG_PollenS3 <- df1$oRG_PollenS3[idx]

  # RG IP Neuron marker expression
  idx <- match(df$CellID, row.names(df1))
  df$RG <- df1$RG[idx]
  df$RG_PollenS3 <- df1$RG_PollenS3[idx]
  df$IP <- df1$IP[idx]
  df$Neuron <- df1$Neuron[idx]

  # Seurat S phase and G2M scores
  idx <- match(df$CellID, row.names(centSO@meta.data))
  df$G2M_Score <- centSO@meta.data$G2M.Score[idx]
  df$S_Score <- centSO@meta.data$S.Score[idx]

  print(range(df$RG))
  print(range(df$IP))
  print(range(df$Neuron))

  return(df)
}
################################################################################

### Percent of cells passing expression filters per cluster

Mixed_Marker_By_Cluster_Percent_Barplot <- function(
  exM, seuratO, highThreshold, lowThreshold, title, cluster_order = NULL) {

  df <- Average_MarkersExp_Per_Cell(noCentExM, centSO)
  df <- Positive_Negative_Expression_Flag(
    exDF = df, highThreshold = 0.5, lowThreshold = 0.5)
  df$TYPE[df$TYPE %in% c("Neuron+", "IP+", "RG+", "RG+ vRG- oRG-")] <- "NA"
  df$TYPE <- droplevels(df$TYPE)
  df <- aggregate(df$TYPE, list(df$CLUSTER), Percent_Of_Table)

  df2 <- as.data.frame(df[ ,2])
  df2$CLUSTER <- df[ ,1]
  df <- melt(df2)

  # Set cluster order
  if (! is.null(cluster_order)){
    df$CLUSTER <- factor(df$CLUSTER, levels = cluster_order)
  }
  df <- df[! is.na(df$CLUSTER), ]

  gg <- ggplot(df, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    xlab("Cluster") +
    ylab("Percent of cells") +
    ggtitle(title)

  return(gg)
}

# Plot
gg1 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.75, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.75 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of RG, oRG, vRG, IPC expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot.pdf")
  , width = 12, height = 14)

# Plot with reordered clusters
gg1 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , cluster_order = c(9,7,8,10,2,0,1,12,4,3,14,5,6,11,13,15,16)
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , cluster_order = c(9,7,8,10,2,0,1,12,4,3,14,5,6,11,13,15,16)
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Mixed_Marker_By_Cluster_Percent_Barplot(exM = noCentExM, seuratO = centSO
  , highThreshold = 0.75, lowThreshold = 0.25
  , cluster_order = c(9,7,8,10,2,0,1,12,4,3,14,5,6,11,13,15,16)
  , title = "Keep CC\n+ = > 0.75 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of RG, oRG, vRG, IPC expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot_paper.pdf")
  , width = 11, height = 11)
################################################################################

### Mixed color tSNE

Mixed_tSNE_Format_Data_For_GGplot <- function(green_genes, red_genes){

  genes1 <- green_genes
  genes2 <- red_genes

  # Collect tSNE and expression values
  tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  exp1_DF <- Mean_Expression(tsneDF, genes = genes1, exM = noCentExM)
  tsneDF$Expression_1 <- exp1_DF$EXPRESSION
  exp2_DF <- Mean_Expression(tsneDF, genes = genes2, exM = noCentExM)
  tsneDF$Expression_2 <- exp2_DF$EXPRESSION

  # Expression limits
  tsneDF$Expression_1[tsneDF$Expression_1 > 1] <- 1
  tsneDF$Expression_1[tsneDF$Expression_1 < 0] <- 0
  tsneDF$Expression_2[tsneDF$Expression_2 > 1] <- 1
  tsneDF$Expression_2[tsneDF$Expression_2 < 0] <- 0

  # Transform expression to 1-255 range for rgb function
  tsneDF$Expression_1 <- tsneDF$Expression_1/1
  tsneDF$Expression_2 <- tsneDF$Expression_2/1
  tsneDF$Expression_1 <- round(tsneDF$Expression_1 * 255, 0)
  tsneDF$Expression_2 <- round(tsneDF$Expression_2 * 255, 0)

  # alpha <- rowMeans(tsneDF[c("Expression_1", "Expression_2")])
  # alpha <- alpha[alpha < 0.25] <- 0.25
  # alpha <- round(alpha * 255, 0)

  # Convert to rgb
  tsneDF <- within(tsneDF
    , mix <- rgb(green = Expression_1, red = Expression_2, blue = 0
      , maxColorValue = 255
      # , alpha = alpha
    )
  )
  return(tsneDF)
}

# RG IP
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "RG"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "IP"]
)
gg1 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = RG; Red = IP")
# RG Neuron
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "RG"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "Neuron"]
)
gg2 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = RG; Red = Neuron")
# IP Neuron
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "IP"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "Neuron"]
)
gg3 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = IP; Red = Neuron")
# Combine
Plot_Grid(list(gg1, gg2, gg3), ncol = 3, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nTwo channel mixed expression of marker genes"
  )
)
ggsave(paste0(outGraph, "tSNE_Mixed.png")
  , width = 12, height = 5)
################################################################################

### Number / Percent of cells in CC phase subset by markers
print("### Number / Percent of cells in CC phase subset by markers")

## Percent of cells in CC phase
df1 <- Average_MarkersExp_Per_Cell(
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

df1 <- Average_MarkersExp_Per_Cell(
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

### Number of cells positive for RG, IP, or Neuron markers

df <- data.frame(Cell_Subset = rep(NA, nrow(df1)))
row.names(df) <- row.names(df1)

df$Cell_Subset <- NA
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$Neuron > 0.5]] <- "Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5]] <- "IP"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$Neuron > 0.5]] <- "RG Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5]] <- "RG IP"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5 & df1$Neuron > 0.5]] <- "IP Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron > 0.5]] <- "RG IP Neuron"

table(df)
################################################################################

### DE RG(cluster 7,9) vs IP(2); RG(7,9) vs Neuron(0); IP(2) vs Neuron(0)

DE_Clusters_Vs_Clusters <- function(clusters1, clusters2) {
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
  deDF <- data.frame(Log_FC_C2vC1 = deLM$coefmat[ ,"clustersclusters2"]
    , Pvalue = deLM$pvalmat[ ,"clustersclusters2"])
  deDF <- deDF[order(deDF$Log_FC), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method="BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))
  return(deDF)
}

DE_By_CellIDs <- function(cell_IDs_1, cell_IDs_2) {
  ids <- c(cell_IDs_1, cell_IDs_2)
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% ids]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% ids
    , c("nUMI", "librarylab", "individual", "CELL")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$groups <- "cell_IDs_1"
  termsDF$groups[termsDF$CELL %in% cell_IDs_2] <- "cell_IDs_2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ groups+nUMI+librarylab+individual")
  # Format LM DE
  deDF <- data.frame(
    Log2_FC_Group1_vs_Group2 = deLM$coefmat[ ,"groupscell_IDs_2"]
    , Pvalue = deLM$pvalmat[ ,"groupscell_IDs_2"]
  )
  deDF$Gene = row.names(deDF)
  # Make DE cell_IDs_1 positive fold change
  deDF$Log2_FC_Group1_vs_Group2 <- deDF$Log2_FC_Group1_vs_Group2 * -1
  # Order by fold change
  deDF <- deDF[order(deDF$Log2_FC_Group1_vs_Group2), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method="BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))
  return(deDF)
}

# DE cell types
# RG vs IP clusters
de_RG_v_IP_DF <- DE_Clusters_Vs_Clusters(clusters1 = 2, clusters2 = c(7,9))
# RG vs Migrating Neuron clusters
de_RG_v_Ne_DF <- DE_Clusters_Vs_Clusters(clusters1 = 0, clusters2 = c(7,9))
# IP vs Migrating Neuron clusters
de_IP_v_Ne_DF <- DE_Clusters_Vs_Clusters(clusters1 = 0, clusters2 = c(2))

# Format and save
de_RG_v_IP_DF$Comparison <- "RG_vs_IP"
de_RG_v_Ne_DF$Comparison <- "RG_vs_Neuron"
de_IP_v_Ne_DF$Comparison <- "IP_vs_Neuron"
deDF <- rbind(de_RG_v_IP_DF, de_RG_v_Ne_DF, de_IP_v_Ne_DF)
# Save as csv
write.csv(deDF, file = paste0(outTable, "DE_CellTypes.csv")
  , quote = FALSE)


Format_Cell_IDs_1_2 <- function(cellIDs1, cellIDs2, comparison_label){
  de_to_run_DF <- data.frame(
    Cell_IDs = c(cellIDs1, cellIDs2)
    , Group = c(rep("cell_IDs_1", length(cellIDs1))
      , rep("cell_IDs_2", length(cellIDs2)))
    , Comparison = comparison_label
  )
  return(de_to_run_DF)
}

# RG+ vs RG+IP+
# S phase
Subset_Cell_IDs_RG_RGIP_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.25 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 8])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 8])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2, "RG_vs_RGIP_Sphase")
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_RG_RGIP_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.25 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 10])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 10])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2, "RG_vs_RGIP_G2Mphase")
  return(de_to_run_DF)
}

# RG+ vs RG+Neuron+
# S phase
Subset_Cell_IDs_RG_RGNeuron_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.25 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 8])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.25 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 8])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2
    , "RG_vs_RGNeuron_Sphase")
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_RG_RGNeuron_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.25 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 10])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 10])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2
    , "RG_vs_RGNeuron_G2Mphase")
  return(de_to_run_DF)
}

# IP+ vs IP+Neuron+
# S phase
Subset_Cell_IDs_IP_IPNeuron_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.25 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 8])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.25 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 8])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2
    , "IP_vs_IPNeuron_Sphase")
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_IP_IPNeuron_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.25 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.25
    ]
  cellIDs1 <- intersect(cellIDs1, names(centSO@ident)[centSO@ident == 10])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.25 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(cellIDs2, names(centSO@ident)[centSO@ident == 10])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(cellIDs1, cellIDs2
    , "IP_vs_IPNeuron_G2Mphase")
  return(de_to_run_DF)
}

# Run DE for transition states
de_to_run_DF <- rbind(
  Subset_Cell_IDs_RG_RGIP_Sphase()
  , Subset_Cell_IDs_RG_RGIP_G2Mphase()
  , Subset_Cell_IDs_RG_RGNeuron_Sphase()
  , Subset_Cell_IDs_RG_RGNeuron_G2Mphase()
  , Subset_Cell_IDs_IP_IPNeuron_Sphase()
  , Subset_Cell_IDs_IP_IPNeuron_G2Mphase()
)
# Format DE into table
transition_state_DE_DF <- do.call("rbind"
  , lapply(unique(de_to_run_DF$Comparison)
    , function(comparison){
    cellIDs1 <- de_to_run_DF$Cell_IDs[
      de_to_run_DF$Group == "cell_IDs_1"
        & de_to_run_DF$Comparison == comparison]
    cellIDs2 <- de_to_run_DF$Cell_IDs[
      de_to_run_DF$Group == "cell_IDs_2"
        & de_to_run_DF$Comparison == comparison]
    deDF <- DE_By_CellIDs(cell_IDs_1 = cellIDs1, cell_IDs_2 = cellIDs2)
    deDF$Comparison <- comparison
    return(deDF)
  })
)
# Save as csv
write.csv(transition_state_DE_DF, file = paste0(outTable, "DE_CellTransitionStates.csv")
  , quote = FALSE)



# DE RG intersection DE RG from RG+ vs IP+
# versus DE RG+G1 vs RG+Sphase
# DE RG+IP+ intersection DE IP from RG+ vs IP+

Percent_Signature_DE_Genes <- function(Comparison_1, Comparison_2){
  genes1 <- row.names(deDF)[
    deDF$Log_FC_C2vC1 > 0.25 &
      deDF$FDR < 0.05 &
      deDF$Comparison == Comparison_1]
  genes2 <- row.names(transition_state_DE_DF)[
    transition_state_DE_DF$Log2_FC_Group1_vs_Group2 > 0.25 &
      # transition_state_DE_DF$FDR < 0.05 &
      transition_state_DE_DF$Comparison == Comparison_2]
  percent_signature <- length(intersect(genes1, genes2))/length(genes1) * 100
  return(percent_signature)
}

Percent_Signature_DE_Genes(
  Comparison_1 = "RG_vs_IP"
  , Comparison_2 = "RG_vs_RGIP_Sphase"
)
Percent_Signature_DE_Genes(
  Comparison_1 = "RG_vs_IP"
  , Comparison_2 = "RG_vs_RGIP_G2Mphase"
)

Percent_Signature_DE_Genes(
  Comparison_1 = "RG_vs_Neuron"
  , Comparison_2 = "RG_vs_RGIP_Sphase"
)
Percent_Signature_DE_Genes(
  Comparison_1 = "RG_vs_Neuron"
  , Comparison_2 = "RG_vs_RGNeuron_G2Mphase"
)

Percent_Signature_DE_Genes(
  Comparison_1 = "IP_vs_Neuron"
  , Comparison_2 = "RG_vs_IPNeuron_Sphase"
)
Percent_Signature_DE_Genes(
  Comparison_1 = "IP_vs_Neuron"
  , Comparison_2 = "RG_vs_IPNeuron_G2Mphase"
)



genes1 <- row.names(deDF)[
  deDF$Log_FC_C2vC1 > 0.25 &
    deDF$FDR < 0.05 &
    deDF$Comparison == "RG_vs_IP"]
genes2 <- row.names(transition_state_DE_DF)[
  transition_state_DE_DF$Log2_FC_Group1_vs_Group2 > 0.25 &
    # transition_state_DE_DF$FDR < 0.05 &
    transition_state_DE_DF$Comparison == "RG_vs_RGIP_Sphase"]
length(intersect(genes1, genes2))/length(genes1)

genes1 <- row.names(deDF)[
  deDF$Log_FC_C2vC1 < -0.25 &
    deDF$FDR < 0.05 &
    deDF$Comparison == "RG_vs_IP"]
genes2 <- row.names(transition_state_DE_DF)[
  transition_state_DE_DF$Log2_FC_Group1_vs_Group2 < -0.25 &
    # transition_state_DE_DF$FDR < 0.05 &
    transition_state_DE_DF$Comparison == "RG_vs_RGIP_Sphase"]
length(intersect(genes1, genes2))/length(genes1)


# de_RG_v_IP_DF <- read.csv(paste0(outTable, "DE_RGvsIP.csv")
#   , header = TRUE, row.names = 1)
# de_RG_v_Ne_DF <- read.csv(paste0(outTable, "DE_RGvsNeuron.csv")
#   , header = TRUE, row.names = 1)
# de_IP_v_Ne_DF <- read.csv(paste0(outTable, "DE_IPvsNeuron.csv")
#   , header = TRUE, row.names = 1)
################################################################################

### ME of DE genes in RG+ IP+ Neuron+

ME_CellType_Genes <- function(
  cellIDs, deDF, mk_exp_DF, foldChange, Class1_label, Class2_label) {

  df2 <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

  names(foldChange) <- foldChange
  genesDFL <- lapply(foldChange, function(foldChange){

    genesDF <- data.frame(Genes = row.names(df2)
      , Class = rep("Neither", length(row.names(df2))))
    genesDF$Class[genesDF$Genes %in% row.names(deDF)[
      deDF$Log_FC_C2vC1 > foldChange & deDF$FDR < 0.05]] <- Class2_label
    genesDF$Class[genesDF$Genes %in% row.names(deDF)[
      deDF$Log_FC_C2vC1 < -(foldChange) & deDF$FDR < 0.05]] <- Class1_label

    return(genesDF)
  })

  ggL <- lapply(names(genesDFL), function(name){

    genesDF <- genesDFL[[name]]

    print(str(genesDF))

    print("ME calculation...")
    df3 <- df2[genesDF$Class != "Neither", ]
    genesDF <- genesDF[genesDF$Class != "Neither", ]
    meM <- moduleEigengenes(t(df3), genesDF$Class)$eigengenes

    print("Formatting...")
    meM$CellID <- colnames(df3)
    df4 <- Marker_Expression_Flag(meM, mk_exp_DF)
    df4 <- df4[ ,c(1:6)]
    df4 <- melt(df4)

    print("T-test...")
    pvals_DF <- sapply(split(df4, df4$variable), function(df1) {
      combinations <- combn(unique(df1$Cell_Subset), 2)
      pvals <- apply(combinations, 2, function(y){
        results <- t.test(df1$value[df1$Cell_Subset == y[1]]
          , df1$value[df1$Cell_Subset == y[2]])
        return(results$p.value)
      })
      names(pvals) <- paste(combinations[1, ], "vs", combinations[2, ], "\npvalue")
      return(round(pvals, 3))
    })
    pvals_DF <- as.data.frame(pvals_DF)

    print("Plotting p-value table...")
    # Plot p-value table
    tt <- ttheme_default(base_size = 10,
                         padding = unit(c(2, 4), "mm"))
    ggt <- tableGrob(pvals_DF, theme = tt)

    # Format for ggplot
    df4$Cell_Subset <- factor(df4$Cell_Subset
      , levels = c("RG", "RG IP", "IP", "RG Neuron", "IP Neuron", "Neuron"))

    print("Plotting ME boxplot...")
    gg <- ggplot(df4, aes(x = Cell_Subset, y = value, fill = variable)) +
      geom_boxplot() +
      ylab("ME value") +
      xlab("Cell subset") +
      ggtitle(paste0("DE log2 fold change cutoff: ", name
        , "\nNumber of DE genes: ", nrow(df3)))

    gg <- plot_grid(gg, ggt, ncol = 2, rel_widths = c(1,0.5))
    return(gg)
  })

  return(ggL)
}

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_RG_v_IP_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "IP", Class2_label = "RG", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE RG vs IP genes in cluster 8 RG+, IP+, RG+IP+, Neuron-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPpNn_Cluster8.pdf"), height = 9)

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_RG_v_IP_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "IP", Class2_label = "RG", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE RG vs IP genes in cluster 10 RG+, IP+, RG+IP+, Neuron-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPpNn_Cluster10.pdf"), height = 9)

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,8)])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_RG_v_Ne_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "Neuron", Class2_label = "RG", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE RG vs Neuron genes in cluster 0,8 RG+, Neuron+, RG+Neuron+, IP-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPnNp_Cluster08.pdf"), height = 9)

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,10)])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_RG_v_Ne_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "Neuron", Class2_label = "RG", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE RG vs Neuron genes in cluster 0,10 RG+, Neuron+, RG+Neuron+, IP-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPnNp_Cluster010.pdf"), height = 9)

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,8)])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_IP_v_Ne_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "Neuron", Class2_label = "IP", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE IP vs Neuron genes in cluster 0,8 IP+, Neuron+, IP+Neuron+, RG-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGnIPpNp_Cluster08.pdf"), height = 9)

mk_exp_DF <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(mk_exp_DF)[
  mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,10)])
ggL <- ME_CellType_Genes(
  cellIDs = cellIDs, deDF = de_IP_v_Ne_DF, mk_exp_DF = mk_exp_DF
  , Class1_label = "Neuron", Class2_label = "IP", foldChange = c(0, 0.25, 0.5)
)
Plot_Grid(ggL, ncol = 1, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nME of DE IP vs Neuron genes in cluster 0,10 IP+, Neuron+, IP+Neuron+, RG-")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_RGnIPpNp_Cluster010.pdf"), height = 9)


## Average Expression
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[
#   df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
# cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
# df2 <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]
# df2 <- melt(df2)
# df2$DE <- "Neither"
# df2$DE[df2$Var1 %in% row.names(deDF)[
#   deDF$Log_FC_C2vC1 > 1 & deDF$FDR < 0.05]] <- "RG"
# df2$DE[df2$Var1 %in% row.names(deDF)[
#   deDF$Log_FC_C2vC1 < -1 & deDF$FDR < 0.05]] <- "IP"
# df2$CellID <- df2$Var2
# df2 <- Marker_Expression_Flag(df2, df1)
# aggregate(value~DE+Cell_Subset, df2, mean)
# ggplot(df2, aes(x = Cell_Subset, y = value, fill = DE)) +
#   geom_boxplot() +
#   geom_violin() +
#   ylim(-0.5, 2)
# ggsave(paste0(outGraph, "DE_boxplot.png"))
################################################################################

### Correlation of cluster mean expression profile and RG IP cells
print("### Correlation of cluster EG and RG IP cells")

df1 <- Average_MarkersExp_Per_Cell(
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

### tSNE of subsets of cells

tSNE_CCgenes_Markers <- function(cellIDs, title, df1) {

  # Subset to cells of interest
  ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
  ssCentSO@raw.data <- ssCentSO@raw.data[
    ,colnames(ssCentSO@raw.data) %in% cellIDs]
  ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

  # ssCentSO@data <- ssNoCentExM
  # ssCentSO <- ScaleData(ssCentSO)

  # tSNE
  ssCentSO <- RunTSNE(ssCentSO, dims.use = 1:8, do.fast = TRUE, perplexity = 20)

  # Format Macosko cell cycle genes table
  genesGroupDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(genesGroupDF) <- c("Grouping", "Gene.Symbol")
  genesGroupDF$Gene.Symbol <- gsub(" *", "", genesGroupDF$Gene.Symbol)
  genesGroupDF <- genesGroupDF[! genesGroupDF$Gene.Symbol == "", ]
  genesGroupDF <- genesGroupDF[! is.na(genesGroupDF$Grouping), ]
  genesGroupDF$Grouping <- gsub(" *", "", genesGroupDF$Grouping)
  genesGroupDF$Grouping <- factor(genesGroupDF$Grouping, levels = unique(genesGroupDF$Grouping))

  # Gene groups to plot - add known markers to cell cycle genes table
  genesGroupDF <- rbind(genesGroupDF
    , kmDF[kmDF$Grouping %in% c("RG", "IP", "Neuron", "vRG", "oRG") ,c(3,2)])

  # Feature plot
  ggL <- FeaturePlot(
    genes = genesGroupDF$Gene.Symbol
    , tsneDF = as.data.frame(ssCentSO@dr$tsne@cell.embeddings)
    , seuratO = ssCentSO
    , exM = ssCentSO@scale.data
    , centScale = TRUE
    , limLow = -1.5, limHigh = 1.5
    , geneGrouping = genesGroupDF$Grouping)
  ggL <- lapply(ggL, function(gg){
    gg + geom_point(size = 0.5)
  })

  df2 <- as.data.frame(ssCentSO@dr$tsne@cell.embeddings)
  df2$CellID <- row.names(df2)
  df2 <- Marker_Expression_Flag(df2, df1)
  print(str(df2))
  gg <- ggplot(df2, aes(x = tSNE_1, y = tSNE_2, color = Cell_Subset)) +
    geom_point(size = 0.5)
    # xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
    # ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
  ggL <- append(ggL, list(gg))
  # gg <- ggplot(df2, aes(x = tSNE_1, y = tSNE_2, color = Cell_Subset_075)) +
  #   geom_point()
  #   # xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
  #   # ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
  # ggL <- append(ggL, list(gg))

  pg <- Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v'
    , axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\n", title
      , "\nExpression of CC phase genes and cell type marker genes"
      , "\nNormalized centered scaled expression"
      , "\nGene list from Macosko et al. 2016")
    )
    return(pg)
}

# Identify RG+ and / or IP+ but Neuron- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
# tSNE and plot
tSNE_CCgenes_Markers(cellIDs = cellIDs
  , title = "Subset cluster 8 RG+ IP+ Neuron- (>0.5, >0.5, <0.25)"
  , df1 = df1)
ggsave(paste0(outGraph, "RG_IP_Neuronn_Cluster8_FeaturePlot.png")
  , width = 12, height = 12, limitsize = FALSE)

# Identify RG+ and / or IP+ but Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
# tSNE and plot
tSNE_CCgenes_Markers(cellIDs = cellIDs
  , title = "Subset cluster 8 RG+ IP+ Neuron- (>0.5, >0.5, <0.25)"
  , df1 = df1)
ggsave(paste0(outGraph, "RG_IP_Neuronn_Cluster10_FeaturePlot.png")
  , width = 12, height = 9, limitsize = FALSE)
################################################################################

### Monocle

# Monocle
# load("../analysis/analyzed_data/Monocle/Monocle_PC1-40_monocleO.Robj")

Monocle_Run <- function(cellIDs, max_components) {

  print("Monocle: running...")

  # Subset to cells of interest
  exDF <- centSO@raw.data[ ,colnames(centSO@raw.data) %in% cellIDs]
  metDF <- centSO@meta.data[row.names(centSO@meta.data) %in% cellIDs, ]
  # Sort order of cells in exDF to match metadata
  exDF <- exDF[ ,match(row.names(metDF), colnames(exDF))]

  # Format data for monocle
  feature_data = data.frame(gene_short_name = rownames(exDF))
  rownames(feature_data) = feature_data$gene_short_name
  pd <- new("AnnotatedDataFrame", data = metDF)
  fd <- new("AnnotatedDataFrame", data = feature_data)

  # Initialize monocle object
  ssMo <- newCellDataSet(cellData = as(as.matrix(exDF), "sparseMatrix"),
    phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size())

  print("Monocle: Metadata dimensions...")
  print(dim(pData(ssMo)))

  # Filter genes and cells
  # retains all
  ssMo <- detectGenes(ssMo, min_expr = 0.1)
  print(dim(fData(ssMo)))
  # genes expressed in at least 10 cells
  expressed_genes <- row.names(subset(fData(ssMo), num_cells_expressed >= 10))
  ssMo <- ssMo[expressed_genes, ]
  print(dim(fData(ssMo)))
  pData(ssMo)$Total_mRNAs <- Matrix::colSums(exprs(ssMo))

  print("Monocle: Size factors and dispersion")
  # Size factors and dispersion
  ssMo <- estimateSizeFactors(ssMo)
  ssMo <- estimateDispersions(ssMo)

  print("Monocle: Dispersed genes to use for pseudotime ordering")
  # Dispersed genes to use for pseudotime ordering
  disp_table <- dispersionTable(ssMo)
  ordering_genes <- subset(disp_table
    , mean_expression >= 0.01 & dispersion_empirical >= 0.75 * dispersion_fit)$gene_id
  print(length(expressed_genes))
  ssMo_filtered <- setOrderingFilter(ssMo, ordering_genes)
  # Plot
  # png(paste0(outGraph, "Monocle_OrderingGenesDispersion_Cluster.png"))
  #   # , paste0(clid, collapse = '-'), ".png"))
  # plot_ordering_genes(ssMo_filtered)
  # dev.off()

  # print("Variance explained by each PC")
  # # Variance explained by each PC
  # png(paste0(outGraph, "Monocle_PCA_VarianceExplained_Cluster"
  #   , title, ".png"))
  # plot_pc_variance_explained(ssMo_filtered, verbose = TRUE
  #   , use_existing_pc_variance = TRUE, return_all = FALSE)
  # dev.off()

  print("Monocle: Reduce data dimensionality")
  # Reduce data dimensionality
  # Use number of genes expressed or total mRNAs?
  ssMo_filtered <- reduceDimension(ssMo_filtered
    , max_components = max_components
    , residualModelFormulaStr = "~individual + librarylab + Total_mRNAs"
    , verbose = TRUE)

  print("Monocle: Order cells along trajectory")
  # Order cells along trajectory
  ssMo_filtered <- orderCells(ssMo_filtered)

  return(ssMo_filtered)

  ## Differential gene test as function of pseudotime
}

Monocle_Plot_Trajectory <- function(ssMo_filtered, title){
  plot_cell_trajectory(ssMo_filtered, 1, 2, color_by = "Pseudotime"
    , cell_size = 1) +
    scale_color_viridis() +
    theme(legend.position = "right") +
    ggtitle(title)
    # ggtitle(paste0("Cluster: ", cluster))
  ggsave(paste0(outGraph, "Monocle_Trajectory_Pseudotime_Comp12", title, ".png")
    , width = 7, height = 7, limitsize = FALSE)
}

Format_Macosko_Cell_Cycle_Genes_Table <- function(ccDF) {
  # Subset to genes
  # Format Macosko cell cycle genes table
  genesGroupDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(genesGroupDF) <- c("Grouping", "Gene.Symbol")
  genesGroupDF$Gene.Symbol <- gsub(" *", "", genesGroupDF$Gene.Symbol)
  genesGroupDF <- genesGroupDF[! genesGroupDF$Gene.Symbol == "", ]
  genesGroupDF <- genesGroupDF[! is.na(genesGroupDF$Grouping), ]
  genesGroupDF$Grouping <- gsub(" *", "", genesGroupDF$Grouping)
  genesGroupDF$Grouping <- factor(genesGroupDF$Grouping,
     levels = unique(genesGroupDF$Grouping))
  return(genesGroupDF)
}

Format_Genes_Table_To_Plot <- function(ccDF, kmDF, deDF) {
  # Add known markers to cell cycle genes table
  genesGroupDF <- rbind(ccDF
    , kmDF[kmDF$Grouping %in% c("RG", "IP", "Neuron") ,c(3,2)])

  # Add RG IP DE genes
  genesGroupDF <- rbind(genesGroupDF,
    data.frame(Gene.Symbol = row.names(deDF)[
      deDF$Log_FC_C2vC1 > 0 & deDF$FDR < 0.05]
      , Grouping = "RGvsIP DE Genes - RG")
    )
  genesGroupDF <- rbind(genesGroupDF,
    data.frame(Gene.Symbol = row.names(deDF)[
      deDF$Log_FC_C2vC1 < 0 & deDF$FDR < 0.05]
      , Grouping = "RGvsIP DE Genes - IP")
    )

  return(genesGroupDF)
}

Format_GenesGroups_Pseudotime_For_GGplot <- function(genesGroupDF, ssMo_filtered){
  genes <- genesGroupDF$Gene.Symbol
  df1 <- melt(noCentExM[row.names(noCentExM) %in% genes
    , colnames(noCentExM) %in% pData(ssMo_filtered)$CELL])
  df1$Pseudotime <- pData(ssMo_filtered)$Pseudotime[
    match(df1$Var2, pData(ssMo_filtered)$CELL)]
  df1$Grouping <- genesGroupDF$Grouping[
    match(df1$Var1, genesGroupDF$Gene.Symbol)]
  return(df1)
}


GGplot_GenesGroups_Pseudotime <- function(ggDF, title) {
  ggplot(ggDF, aes(x = Pseudotime, y = value, color = Grouping)) +
    geom_smooth() +
    ggtitle(title)
  ggsave(paste0(outGraph, "Monocle_MarkerExpr_Pseudotime", title, ".png"))
}

Plot_Monocle_Analysis <- function(mo, ccDF, kmDF, deDF, title){
  ccDF <- Format_Macosko_Cell_Cycle_Genes_Table(ccDF)
  genesGroupDF <- Format_Genes_Table_To_Plot(ccDF, kmDF, deDF)
  ggDF <- Format_GenesGroups_Pseudotime_For_GGplot(
    genesGroupDF, ssMo_filtered = mo
  )
  GGplot_GenesGroups_Pseudotime(ggDF, title)
}

# Identify RG+ and / or IP+ but Neuron- cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = 4)
Plot_Monocle_Analysis(
  mo = ssMO
  , ccDF = ccDF
  , kmDF = kmDF
  , deDF = de_RG_v_IP_DF
  , title = "Cluster810_RGp_IPp_Nn_Comp4")
)



# Identify RG+ and / or IP- but Neuron+ cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
Monocle_Run(cellIDs = cellIDs, title = "Cluster810_RGp_IPn_Np")

# Identify RG+ and / or IP+ but Neuron- cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
Monocle_Run(cellIDs = cellIDs, title = "Cluster810_RGp_IPp_Nn")

# Identify RG+ and / or RG+IP+ but Neuron- cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25 & df1$RG > 0.5]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
Monocle_Run(cellIDs = cellIDs, title = "Cluster810_RGp_RGpIPp_Nn")

# Identify RG+ and / or RG+Neuron+ but IP- cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25 & df1$RG > 0.5]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
Monocle_Run(cellIDs = cellIDs, title = "Cluster810_RGp_IPn_RGpNp")

# Identify IP+ and / or IP+Neuron+ but RG- cluster 8, 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$IP > 0.5 & df1$RG < 0.25 | df1$Neuron > 0.5 & df1$RG < 0.25 & df1$IP > 0.5]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident %in% c(8,10)]))
# Run monocle
Monocle_Run(cellIDs = cellIDs, title = "Cluster810_RGn_IPp_IPpNp")


# print("Differential gene test as function of pseudotime")
#
# ptDE <- differentialGeneTest(ssMo_filtered
#   , fullModelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
#   , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")


#
#
# diff_test_res <- ptDE
#
# # Subset to genes and by qval
# genes <- kmDF$Gene.Symbol[
#   kmDF$Grouping %in% c("RG", "IP", "Neuron", "vRG", "oRG")]
#
# diff_test_res <- diff_test_res[diff_test_res$gene_short_name %in% genes, ]
#
# genes <- row.names(diff_test_res)
#
# png(paste0(outGraph, "Monocle_DEpseudoTime_Heatmap.png")
#   , width = 9, height = 6, units = "in", res = 300)
# tryCatch(plot_pseudotime_heatmap(ssMo_filtered[genes, ]
#   , num_clusters = 3
#   , cores = 1
#   , show_rownames = T)
#   , error = function(cond) {
#     message("Error message:")
#     message(cond)
#     # Choose a return value in case of error
#     return(NULL)
#   }
# )
# dev.off()
################################################################################

# ### Heatmap and heirarchical clustering of RG and IP cells
# print("### Heatmap and heirarchical clustering of RG and IP cells")
#
# ## Heatmap and heirarchical clustering of RG and IP G2M cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.5 | df1$IP > 0.5 & df1$Neuron < 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "G2M"])
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-2, 2, by = 0.1)
#
# png(paste0(outGraph, "RG_IP_G2M_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
# )
# dev.off()
#
#
# ## Heatmap and heirarchical clustering of RG cluster 10 cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
# df1 <- df1[row.names(df1) %in% cellIDs, ]
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# # annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# # annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-3, 3, by = 0.1)
#
# png(paste0(outGraph, "RG_Cluster10_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
#   , treeheight_col = 150
#   , cutree_cols = 10
# )
# dev.off()
#
#
# ## Heatmap and heirarchical clustering of RG cluster 8 cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
# df1 <- df1[row.names(df1) %in% cellIDs, ]
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# # annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# # annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-3, 3, by = 0.1)
#
# png(paste0(outGraph, "RG_Cluster8_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
#   , treeheight_col = 150
#   , cutree_cols = 10
# )
# dev.off()
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
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IP_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or IP+ but Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IP_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuron_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuron_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ IP- Neuron- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+  IP- Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]


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

# Subset to top expressed genes
exLM <- lapply(exLM, function(exM){
  mns <- rowMeans(exM)
  genes <- names(sort(mns, decreasing = TRUE))[1:10000]
  # genes <- names(sort(mns, decreasing = TRUE))[1:1000]
  exM <- exM[row.names(exM) %in% genes, ]
  return(exM)
})
# Subset to high variance genes
exLM <- lapply(exLM, function(exM){
  varDF <- data.frame(apply(exM, 1, var))
  # ggplot(varDF, aes(x = varDF[,1])) +
  #   geom_histogram(binwidth = 0.1)
  # ggsave(paste0(outGraph, "VarHist.png"))
  exM <- exM[row.names(exM) %in% row.names(varDF)[varDF > 0.3], ]
  return(exM)
})

# PCA
pcaL <- lapply(exLM, function(exM){
  pca <- prcomp(t(exM), center = FALSE)
  print(head((pca$sdev)^2 / sum(pca$sdev^2)*100))
  return(pca)
})

# Plot PCA

outGraphPCA <- paste0(
  dirname(outGraph)
  , "/PCA_Top10000_Var03_NoCent/"
  , basename(outGraph)
)
dir.create(dirname(outGraphPCA), recursive = TRUE)

# Plot PCA loadings
lapply(names(pcaL), function(name){
  pca <- pcaL[[name]]
  Prcomp_Loadings_Plot(pca = pca, nGenes = c(1:20), nPCs = c(1:8)
    , title = paste0(graphCodeTitle, "\n\n", name, "\nGenes with highest PC loadings"))
  ggsave(paste0(outGraphPCA, name, "_PCAloadings.pdf"), width = 13
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
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel05.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset_025")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.25 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel025.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset_075")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.75 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel075.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG_oRG_Subset")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG+ and / or oRG+"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRGoRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "S_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by S phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Sscore.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "G2M_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by G2M phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_G2Mscore.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "oRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_oRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "oRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_oRGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "RG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_RG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "RG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_RGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "IP", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_IP.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Neuron.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "IP", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_IP025scale.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Neuron025scale.png"), width = 14, height = 12)

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
  # ggsave(paste0(outGraphPCA, name, "_PCA_CorCluster7.png"), width = 14, height = 12)
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
  # ggsave(paste0(outGraphPCA, name, "_PCA_CorCluster9.png"), width = 14, height = 12)

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

### DE RG vs IP vs IP+ RG+ analysis

## Heatmap of RG vs IP DE genes in IP+ and Neuron- or RG+ and Neuron- cells

# S phase cells

# Subset expression matrix to DE genes and IP+ and Neuron- or RG+ and Neuron- cells
# DE genes
genes <- row.names(deDF)[deDF$FDR < 0.05]
# IP+ and Neuron- or RG+ and Neuron- cells
df1 <- Average_MarkersExp_Per_Cell(
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
df1 <- Average_MarkersExp_Per_Cell(
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
df1 <- Average_MarkersExp_Per_Cell(
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
