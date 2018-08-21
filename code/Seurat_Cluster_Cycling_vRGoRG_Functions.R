# Damon Polioudakis
# 2017-11-01
# Separate cycling vRG and oRG

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
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
  # df$TYPE[apply((df[ ,c("Neuron", "RG", "IP")] > highThreshold), 1, all)] <- "Neuron+ RG+ IP+"

  # Controls
  df$TYPE[df[ ,c("Endothelial")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Endothelial+ IP+"
  df$TYPE[df[ ,c("Interneuron")] > highThreshold &
      df[ ,c("IP")] > highThreshold] <- "Interneuron+ IP+"

  df$TYPE <- factor(df$TYPE, levels = c("Neuron+", "IP+", "RG+"
    , "RG+ vRG- oRG-", "IP+ RG+", "Neuron+ IP+", "Neuron+ RG+"
    , "Endothelial+ IP+", "Interneuron+ IP+", "Neuron+ RG+ IP+"))
  df$CLUSTER <- factor(as.character(df$CLUSTER)
    , levels = sort(unique(as.numeric(as.character(df$CLUSTER)))))

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
    , vRG_PollenS3 = Mean_Gene_Group_Expression(
      exM = exM, grouping = "vRG-PollenS3")
    , oRG = Mean_Gene_Group_Expression(exM = exM, grouping = "oRG")
    , oRG_PollenS3 = Mean_Gene_Group_Expression(
      exM = exM, grouping = "oRG-PollenS3")
    , RG = Mean_Gene_Group_Expression(exM = exM, grouping = "RG")
    , RG_PollenS3 = Mean_Gene_Group_Expression(
      exM = exM, grouping = "RG-PollenS3")
    , IP = Mean_Gene_Group_Expression(exM = exM, grouping = "IP")
    , Endothelial = Mean_Gene_Group_Expression(
      exM = exM, grouping = "Endothelial Cell")
    , Neuron = Mean_Gene_Group_Expression(exM = exM, grouping = "Neuron")
    , Interneuron = Mean_Gene_Group_Expression(
      exM = exM, grouping = "GABAergic interneuron")
  )

  idx <- match(row.names(mnExDF), row.names(seuratO@meta.data))
  mnExDF$PHASE <- seuratO@meta.data$Phase[idx]
  mnExDF$CLUSTER <- seuratO@ident[idx]
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

  df <- Average_MarkersExp_Per_Cell(exM, seuratO)
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
