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
require(ggdendro)
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

# Iossifov ASD - from Luis metaMat
# From the supplmental tables of Iossifov et al., 2014 - note that some gene
# symbols have the excel conversion error in this data, but it will get removed
# when we overlap things
iosDF <- read.csv(
  "../source/metaMat/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/AllDeNovoByGene.csv")

# de Rubeis ASD - from Luis metaMat
rubDF <- read.table(
  "../source/metaMat/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/deRubeis_mutationlist.txt"
  , header = TRUE, sep = "\t")

# TADA Sanders 2015 = from Luis metaMat
tadaDF <- read.csv(
  "../source/metaMat/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv")

# Known marker Luis table
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv")

## Variables
graphCodeTitle <- "Human_Specific_ASD_Expression.R"
outGraph <- "../analysis/graphs/Human_Specific_ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_ASD_Expression_"
outTable <- "../analysis/tables/Human_Specific_ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_ASD_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Human_Specific_ASD_Expression_"
# outGraph <- "../analysis/graphs/Human_Specific_ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Human_Specific_ASD_Expression_"

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

Seurat_Heatmap_By_Cluster_Hclust_Genes <- function(genes) {
  
  # Subset expression matrix
  exM <- centSO@scale.data
  exM <- exM[row.names(exM) %in% genes, ]
  
  # Obtain the dendrogram
  dend <- as.dendrogram(hclust(dist(exM)))
  dend_data <- dendro_data(dend)
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  # Use the dendrogram label data to position the gene labels
  gene_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, gene = as.character(label), height = 1))
  # Limits for the vertical axes
  gene_axis_limits <- with(
    gene_pos_table,
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1
  # Facet dendrogram so it lines up with faceted heatmaps
  segment_data$Facet <- ""
  
  # Dendrogram plot
  plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_x_reverse(expand = c(0, 0.5)) + 
    scale_y_continuous(breaks = gene_pos_table$y_center, 
      labels = gene_pos_table$gene, 
      limits = gene_axis_limits,
      expand = c(0, 0)) + 
    facet_wrap(~Facet) +
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) +
    theme(strip.background = element_blank())
  
  # Heatmap plot
  geneGroupDF <- data.frame(GENE = row.names(exM), GROUP = "")
  ggL <- lapply(c(0:17), function(cluster){
    tryCatch(
      Heatmap_By_Cluster(
        geneGroupDF = geneGroupDF
        , exprM = as.matrix(centSO@scale.data)
        , seuratO = centSO
        , clusters = cluster
        , lowerLimit = -1.5
        , upperLimit = 1.5
        , geneOrder = gene_pos_table$gene
      )
      , error = function(e) NULL)
  })
  # Remove nulls from ggplot list
  ggL <- ggL[! sapply(ggL, is.null)]
  # Extract legend
  legend <- get_legend(ggL[[1]])
  # Format - remove axis labels
  # ggL[[1]] <- ggL[[1]] + theme(
  #   axis.title.x = element_blank()
  #   , legend.position = "none"
  #   , strip.text.y = element_blank())
  ggL[1:length(ggL)] <- lapply(ggL[1:length(ggL)], function(gg) {
    gg + theme(
      strip.text.y = element_blank()
      , legend.position = "none"
      , axis.title.y = element_blank()
      , axis.text.y = element_blank()
      , axis.ticks.y = element_blank()
      , axis.title.x = element_blank()
      # margin: top, right, bottom, and left
      , plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")
    )
  })
  
  # Combine individual heatmaps and dendrogram
  rel_widths <- as.vector(log((table(centSO@ident) + 1), 10)) + 1
  rel_widths <- c(20, rel_widths, 1)
  # Combine
  pg <- plot_grid(plotlist = append(list(plt_dendr), ggL), ncol = 19
    , rel_widths = rel_widths, align = 'h', axis = 't')
  
  return(pg)
}
################################################################################

### Format

## Iossifov ASD
iosDF <- iosDF[!duplicated(iosDF[,1]),]
rownames(iosDF) <- iosDF[,"gene"]
exomeLength <- iosDF[,"codingLenInTarget"] ## For later use in making the covariate matrix
names(exomeLength) <- rownames(iosDF)
iosDF <- iosDF[,c(5:9,12,17:18,23:24,29)] ## Only keeping the desired columns
iosDF[iosDF>1] <- 1 ## If there is more than one count, just count it as 1 - we are using "yes" vs "no" criteria
iosDF = iosDF[,c(6,7)]
iosDF$gene = rownames(iosDF)

## de Rubeis ASD
rownames(rubDF) = rubDF$Gene
rubDF = rubDF[,3:7]
colnames(rubDF) = paste0("Rubeis_", colnames(rubDF))
rubDF = rubDF
rubDF$gene = rownames(rubDF)

## combined de novo proband LOF iossifov + de Rubeis
asdDF = merge(rubDF, iosDF, by = intersect("gene", "gene"))
asdDF$combined_dn_prob_LOF = asdDF$dnv_LGDs_prb + asdDF$Rubeis_dn.LoF

## TADA Sanders_TADA0.1_exomedel
tadaDF <- tadaDF[,1:21]
tada <- unique(tadaDF[tadaDF$tadaFdrAscSscExomeSscAgpSmallDel<0.1, "RefSeqGeneName"])
################################################################################

### ASD combined de novo proband LOF iossifov + de Rubeis

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- asdDF$gene[asdDF$combined_dn_prob_LOF > 0]

# Genes DE >0.4 in cluster 4 or cluster 14
df1 <- rbind(
  data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 3]
    , Cluster = 3)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 4]
    , Cluster = 4)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 5]
    , Cluster = 5)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 6]
    , Cluster = 6)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 14]
    , Cluster = 14)
)

# Subset to ASD
df1 <- df1[df1$Gene %in% genes, ]

# Feature plot
ggL <- FeaturePlot(
  genes = df1$Gene
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
  )
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.05, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of ASD combined de novo proband LOF iossifov + de Rubeis genes DE in clusters of interest"
    , "\nClusters: 3, 4, 5, 6, 14"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
  )
ggsave(paste0(
  outGraph, "DE_IossifovRubeis_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 12, height = 50, limitsize = FALSE)

# Heatmap
# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = df1$Gene, GROUP = df1$Cluster)
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
    , "\n\nExpression of ASD combined de novo proband LOF iossifov + de Rubeis genes genes DE in clusters of interest"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
)
ggsave(paste0(
  outGraph, "DE_IossifovRubeis_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 16, limitsize = FALSE)

# Heatmap + hclust of genes
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(df1$Gene)
# Title
title = paste0(graphCodeTitle
  , "\n\nExpression of ASD combined de novo proband LOF iossifov + de Rubeis genes DE in clusters of interest"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters:"
  , "\n> 0.4 log fold change cluster vs all other cells"
  , "\n")
# now add the title
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(
  outGraph, "DE_IossifovRubeis_HeatmapDend_NormalizedCenteredScaled.png")
  , width = 12, height = 13, limitsize = FALSE)

# Heatmap + hclust of ASD combined de novo proband LOF iossifov + de Rubeis genes
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(genes)
# Title
title = paste0(graphCodeTitle
  , "\n\nExpression of ASD combined de novo proband LOF iossifov + de Rubeis genes"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters:"
  , "\n> 0.4 log fold change cluster vs all other cells"
  , "\n")
# now add the title
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.03, 1))
# Save
ggsave(paste0(
  outGraph, "IossifovRubeis_HeatmapDend_NormalizedCenteredScaled.png")
  , width = 12, height = 50, limitsize = FALSE)
################################################################################

### Human specific ASD combined de novo proband LOF iossifov + de Rubeis

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- asdDF$gene[asdDF$combined_dn_prob_LOF > 0]

# Allen human specific
genes <- genes[genes %in% hsDF$Gene[hsDF$Set == "Human-specific"]]

# Feature plot
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific + ASD combined de novo proband LOF iossifov + de Rubeis genes"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(
  outGraph, "HumanSpecific_IossifovRubeis_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 12, height = 6, limitsize = FALSE)

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
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.4, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen human specific + ASD combined de novo proband LOF iossifov + de Rubeis genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(
  outGraph, "HumanSpecific_IossifovRubeis_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 5, limitsize = FALSE)
################################################################################

### ASD TADA

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- tada

# Genes DE >0.4 in cluster 4 or cluster 14
df1 <- rbind(
  data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 3]
    , Cluster = 3)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 4]
    , Cluster = 4)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 5]
    , Cluster = 5)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 6]
    , Cluster = 6)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 14]
    , Cluster = 14)
)

# Subset to ASD
df1 <- df1[df1$Gene %in% genes, ]

# Feature plot
ggL <- FeaturePlot(
  genes = df1$Gene
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of ASD TADA genes DE in clusters of interest"
    , "\nClusters: 3, 4, 5, 6, 14"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
)
ggsave(paste0(
  outGraph, "DE_TADA_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 12, height = 20, limitsize = FALSE)

# Heatmap
# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = df1$Gene, GROUP = df1$Cluster)
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
    , "\n\nExpression of ASD TADA genes genes DE in clusters of interest"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
)
ggsave(paste0(
  outGraph, "DE_TADA_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 14, limitsize = FALSE)

# Heatmap + hclust of genes
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(df1$Gene)
# Title
title = paste0(graphCodeTitle
  , "\n\nExpression of ASD TADA genes DE in clusters of interest"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters:"
  , "\n> 0.4 log fold change cluster vs all other cells"
  , "\n")
# now add the title
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(
  outGraph, "DE_TADA_HeatmapDend_NormalizedCenteredScaled.png")
  , width = 12, height = 12, limitsize = FALSE)

# Heatmap + hclust of ASD combined de novo proband LOF iossifov + de Rubeis genes
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(genes)
# Title
title = paste0(graphCodeTitle
  , "\n\nExpression of ASD TADA genes"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters:"
  , "\n> 0.4 log fold change cluster vs all other cells"
  , "\n")
# now add the title
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.12, 1))
# Save
ggsave(paste0(
  outGraph, "TADA_HeatmapDend_NormalizedCenteredScaled.png")
  , width = 12, height = 14, limitsize = FALSE)
################################################################################

### Human specific or primate specific + TADA

# ASD TADA
genes <- as.character(tada)

# Allen human specific
# No intersection of TADA and Allen human specific
table(genes %in% as.character(hsDF$Gene[hsDF$Set == "Human-specific"]))
# FALSE
# 66

# Allen primate specific
genes <- genes[genes %in% as.character(hsDF$Gene[hsDF$Set == "Primate-specific"])]

# Feature plot
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = NULL
  , centScale = TRUE
)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen primate specific + ASD TADA genes"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(
  outGraph, "PrimateSpecific_TADA_FeaturePlot_NormalizedCenteredScaled.png")
  , width = 12, height = 6, limitsize = FALSE)

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
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.4, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of Allen primate specific + ASD TADA genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(
  outGraph, "PrimateSpecific_TADA_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 5, limitsize = FALSE)
################################################################################