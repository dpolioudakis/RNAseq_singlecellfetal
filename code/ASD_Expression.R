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
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Cluster DE table
deDF <- read.table(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Nowakowski cluster DE table
nowakowski_DE_DF <- read.csv(
  "../nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)

# Nowakowski expression matrix
nowakowski_ex_DF <- read.table(
  "../nowakowski_2017/geneMatrix.tsv", header = TRUE
)
# Nowakowski metadata
nowakowski_Mt_DF <- read.table(
  "../nowakowski_2017/metaData.txt", header = TRUE, fill = TRUE
)

# Lake 2017 cluster DE tables
lake_DE_DF <- read.csv(
  "../source/Lake_2018_TS3_Cluster_DE_Genes.csv"
  , skip = 4, header = TRUE, fill = TRUE
)

# Lake 2017 expression matrix
lake_ex_DF <- read.table(
  "../lake_2017/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", header = TRUE
)

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

ihartDF <- read.csv("../source/metaMat/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv")

# Known marker Luis table
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv")

## Variables
graphCodeTitle <- "ASD_Expression.R"
outGraph <- "../analysis/graphs/ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/ASD_Expression_WardD2Link_"
outTable <- "../analysis/tables/ASD_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/ASD_Expression_WardD2Link_"

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
  clsNo <- clusterOrder[! c(0:17) %in% okayClusters]
  # Gene is > X FC in cluster
  genes1 <- deDF$Gene[deDF$Log2_Fold_Change > fcHigh & deDF$Cluster == cluster]
  # Genes in clusters genes cannot be DE in > 0.3
  genes2 <- deDF$Gene[deDF$Log2_Fold_Change > fcLow & deDF$Cluster %in% clsNo]
  # Check
  print(table(genes1 %in% genes2))
  # Remove genes in clusters genes cannot be DE in
  genes1 <- genes1[! genes1 %in% genes2]
  # Filter DE DF
  utdeDF <- deDF[deDF$Gene %in% genes1, ]
  utdeDF <- utdeDF[utdeDF$Cluster == cluster, ]
  return(utdeDF)
}

Combine_DE_and_Expression <- function(deDF, exDF) {
  # Column for setting order of genes
  utdeDF$ORDER <- seq(1, nrow(utdeDF))
  # Merge with expression data frame
  ggDF <- merge(utdeDF[c("Cluster", "Cluster", "ORDER")], exDF
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$Cluster <- as.factor(ggDF$Cluster)
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

# hclust + heatmap
Seurat_Heatmap_By_Cluster_Hclust_Genes_SetColWidths <- function(
  genes, clusterOrder) {

  print("Seurat_Heatmap_By_Cluster_Hclust_Genes_SetColWidths")

  # Subset expression matrix
  exM <- centSO@scale.data
  exM <- exM[row.names(exM) %in% genes, ]

  # Obtain the dendrogram
  dend <- as.dendrogram(hclust(d = dist(exM), method = "ward.D2"))
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
  geneGroupDF <- data.frame(GENE = gene_pos_table$gene, GROUP = "")
  gg <- Plot_Marker_Genes_Heatmap_SetColWidths(geneGroupDF = geneGroupDF)

  # Combine heatmap and dendrogram
  pg <- plot_grid(plotlist = list(gg_dendr, gg)
    , ncol = 2
    , rel_widths = c(0.2,1), align = 'h', axis = 't')

  return(pg)
}

Seurat_Heatmap_By_Cluster_Hclust_Genes <- function(genes, clusterOrder) {

  # Subset expression matrix
  exM <- centSO@scale.data
  exM <- exM[row.names(exM) %in% genes, ]

  # Obtain the dendrogram
  dend <- as.dendrogram(hclust(d = dist(exM), method = "ward.D2"))
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
  geneGroupDF <- data.frame(Cluster = gene_pos_table$gene, GROUP = "")
  ggL <- lapply(clusterOrder, function(cluster){
    tryCatch(
      Heatmap_By_Cluster(
        geneGroupDF = geneGroupDF
        , exprM = as.matrix(centSO@scale.data)
        , seuratO = centSO
        , clusters = cluster
        , lowerLimit = -1.5
        , upperLimit = 1.5
        , geneOrder = TRUE
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
  rel_widths <- data.frame(log((table(centSO@ident) + 1), 5))
  rel_widths <- rel_widths[match(clusterOrder, rel_widths$Var1), ]
  rel_widths <- as.vector(rel_widths$Freq) + 1
  rel_widths <- c(20, rel_widths, 1)
  # Combine
  pg <- plot_grid(plotlist = append(list(plt_dendr), ggL), ncol = length(rel_widths)
    , rel_widths = rel_widths, align = 'h', axis = 't')

  return(pg)
}

# Expression heatmap, cells ordered by cluster
Plot_Marker_Genes_Heatmap_SetColWidths <- function(
  geneGroupDF
  , exprM = centSO@scale.data
  , seuratO = centSO
  , clusters = c(0:15)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = TRUE
  , clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  ) {

  print("Plot_Marker_Genes_Heatmap_SetColWidths")

  # Heatmap plot
  # Normalized, mean centered and scaled
  ggDF <- Heatmap_By_Cluster_Format_Data(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , seuratO = seuratO
    , clusters = clusters
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneOrder
    , clusterOrder = clusterOrder
  )

  print("Heatmap_By_Cluster: plotting...")
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = GENE_GROUP, fill = value)) +
    geom_tile() +
    facet_grid(GROUP~SEURAT_CLUSTERS, space = "free_y", scales = "free"
      , drop = TRUE) +
    # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(lowerLimit, upperLimit)) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(strip.background = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10)) +
    ylab("Genes") +
    xlab("Cells ordered by cluster")
  # gg <- gg + ...
  return(gg)
}
################################################################################

### Format

## Lake
lake_DE_DF <- lake_DE_DF[ ,1:7]

## Nowakowski
rownames(nowakowski_ex_DF) <- nowakowski_ex_DF$geneId
nowakowski_ex_DF <- nowakowski_ex_DF[ ,-1]

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

## iHART
ihartDF$HGNC.gene.symbol <- gsub("\"", "", ihartDF$HGNC.gene.symbol)
################################################################################

### Data processing

# Read depth normalize Lake 2017 (ln(transcripts-per-10,000 +1))
for(i in 1:ncol(lake_ex_DF)){
  read_depth <- sum(lake_ex_DF[ ,i])
  lake_ex_DF[ ,i] <- lake_ex_DF[ ,i] / read_depth
}
lake_ex_DF <- (lake_ex_DF*10000)+1
lake_ex_DF <- log(lake_ex_DF)
################################################################################

### ASD combined de novo proband LOF iossifov + de Rubeis

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- asdDF$gene[asdDF$combined_dn_prob_LOF > 0]

# Genes DE >0.4 in cluster 4 or cluster 13
df1 <- rbind(
  data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 3]
    , Cluster = 3)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 4]
    , Cluster = 4)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 5]
    , Cluster = 5)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 6]
    , Cluster = 6)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 13]
    , Cluster = 13)
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
    , "\nClusters: 3, 4, 5, 6, 13"
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
geneGroupDF <- data.frame(Cluster = df1$Gene, GROUP = df1$Cluster)
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
geneGroupDF <- data.frame(Cluster = genes, GROUP = "")
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

# ASD TADA genes
genes <- tada

# Genes DE >0.4 in cluster 4 or cluster 13
df1 <- rbind(
  data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 3]
    , Cluster = 3)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 4]
    , Cluster = 4)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 5]
    , Cluster = 5)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 6]
    , Cluster = 6)
  , data.frame(Gene = deDF$Gene[deDF$Log2_Fold_Change > 0.4 & deDF$Cluster == 13]
    , Cluster = 13)
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
    , "\nClusters: 3, 4, 5, 6, 13"
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
geneGroupDF <- data.frame(Cluster = df1$Gene, GROUP = df1$Cluster)
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

# Heatmap + hclust of ASD TADA genes
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(genes = genes
  , clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
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

## Heatmap of ASD TADA genes with fixed column widths
genes_DF <- deDF[deDF$Cluster %in% c(0,1,3,4,13,5,6), ]
genes_DF <- genes_DF[order(-genes_DF$Log2_Fold_Change), ]
genes <- intersect(genes_DF$Gene, tada)
genes <- rev(genes)
geneGroupDF <- data.frame(
  GENE = genes
  , GROUP = NA
)
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(geneGroupDF = geneGroupDF)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_HeatmapSetColWidths_NormalizedCenteredScaled_paper.png"
  )
  , width = 12, height = 6
)
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
geneGroupDF <- data.frame(Cluster = genes, GROUP = "")
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

### iHART

# iHART novel 17
# Heatmap
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.17novel" == 1]
geneGroupDF <- data.frame(GENE = genes, GROUP = "")
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(geneGroupDF = geneGroupDF)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD iHART novel 17 genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "iHART17novel_HeatmapSetColWidths_NormalizedCenteredScaled.png"
  )
  , width = 12, height = 6
)
# Hclust + heatmap
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes_SetColWidths(
  genes,
  clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
title <- paste0(
  graphCodeTitle
    , "\n\nExpression of ASD iHART novel 17 genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n"
)
title <- ggdraw() + draw_label(title)
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph
    , "iHART17novel_Hclust_HeatmapSetColWidths_NormalizedCenteredScaled.png"
  )
  , width = 12, height = 6
)

# iHART 69
# Hclust + heatmap
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.69" == 1]
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes_SetColWidths(
  genes,
  clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
title <- paste0(
  graphCodeTitle
    , "\n\nExpression of ASD iHART 69 genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n"
)
title <- ggdraw() + draw_label(title)
plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
ggsave(paste0(outGraph
    , "iHART69_Hclust_HeatmapSetColWidths_NormalizedCenteredScaled.png"
  )
  , width = 12, height = 11
)

## iHART 69 pie chart
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.69" == 1]
# DE
subset_DE_DF <- deDF[
  deDF$Log2_Fold_Change > 0.2, c("Gene", "Log2_Fold_Change", "Cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c(7, 9, 8, 10, 2)
  , "Glutamatergic" = c(0, 1, 4, 3, 13)
  , "GABAergic" = c(5, 6)
  , "Glia or support cells" = c(11, 12, 14, 15)
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = noCentExM, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = noCentExM, ex_cutoff = 0.5, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART69_CellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART 69 classified by cell type expression in human fetal brain"
    , "\nDrop-seq human fetal brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART69_CellTypes_Piechart.pdf"))


## iHART 16 novel pie chart
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.17novel" == 1]
genes <- genes[! genes %in% "CACNA2D3"]
# DE
subset_DE_DF <- deDF[
  deDF$Log2_Fold_Change > 0.2, c("Gene", "Log2_Fold_Change", "Cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c(7, 9, 8, 10, 2)
  , "Glutamatergic" = c(0, 1, 4, 3, 13)
  , "GABAergic" = c(5, 6)
  , "Glia or support cells" = c(11, 12, 14, 15)
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = noCentExM, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = noCentExM, ex_cutoff = 0.5, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART16novel_CellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART novel 16 classified by cell type expression in human fetal brain"
    , "\nDrop-seq human fetal brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART16novel_CellTypes_Piechart.pdf"))


## iHART 69 pie chart with Nowakowski dataset
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.69" == 1]
# Check expression levels of iHART genes
subset_ex_M <- nowakowski_ex_DF[rownames(nowakowski_ex_DF) %in% genes, ]
gg_DF <- subset_ex_M
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$Cluster <- nowakowski_Mt_DF$WGCNAcluster[
  match(gg_DF$variable, nowakowski_Mt_DF$X_id)]
gg_DF <- gg_DF[! is.na(gg_DF$Cluster), ]
gg_DF <- gg_DF[! is.na(gg_DF$Gene), ]
ggplot(gg_DF, aes(x = Cluster, y = value)) +
  geom_jitter(size = 0.01) +
  facet_wrap(~Gene)
ggsave(paste0(outGraph, "iHART69_NowkowskiClusters_JitterPlot.png")
  , width = 13, height = 15, dpi = 150)
# Nowakowski DE not very cell type specifc, set DE filter higher
subset_DE_DF <- nowakowski_DE_DF[
  nowakowski_DE_DF$avg_diff > 0.75, c("gene", "avg_diff", "cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c("Astrocyte", "Choroid", "IPC-div1", "RG-div1"
    , "RG-div2", "RG-early", "Endothelial", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3"
    , "IPC-div2", "Microglia", "OPC", "oRG", "tRG", "vRG"
    )
  , "NA" = c("MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "MGE-RG1", "MGE-RG2"
    , "nIN1", "nIN2", "nIN3", "nIN4", "nIN5", "Glyc", "U1", "U2", "U3", "U4"
    , "Mural"
    )
  , "Glutamatergic" = c("EN-PFC2", "EN-PFC3", "EN-V1-2", "EN-PFC1", "EN-V1-1"
    , "EN-V1-3", "nEN-early2", "nEN-early1", "nEN-late"
    )
  , "GABAergic" = c("IN-CTX-CGE", "IN-CTX-CGE2", "IN-CTX-MGE2", "IN-CTX-MGE1"
    , "IN-STR"
    )
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = nowakowski_ex_DF, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = nowakowski_ex_DF, ex_cutoff = 5, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART69_NowkowskiCellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART 69 classified by cell type expression in human fetal brain"
    , "\nNowkowski human fetal brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART69_NowkowskiCellTypes_Piechart.pdf"))

## iHART 16 novel pie chart with Nowakowski dataset
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.17novel" == 1]
genes <- genes[! genes %in% "CACNA2D3"]
# DE filter higher
subset_DE_DF <- nowakowski_DE_DF[
  nowakowski_DE_DF$avg_diff > 1.5, c("gene", "avg_diff", "cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c("Astrocyte", "Choroid", "IPC-div1", "RG-div1"
    , "RG-div2", "RG-early", "Endothelial", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3"
    , "IPC-div2", "Microglia", "OPC", "oRG", "tRG", "vRG"
    )
  , "NA" = c("MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "MGE-RG1", "MGE-RG2"
    , "nIN1", "nIN2", "nIN3", "nIN4", "nIN5", "Glyc", "U1", "U2", "U3", "U4"
    , "Mural"
    )
  , "Glutamatergic" = c("EN-PFC2", "EN-PFC3", "EN-V1-2", "EN-PFC1", "EN-V1-1"
    , "EN-V1-3", "nEN-early2", "nEN-early1", "nEN-late"
    )
  , "GABAergic" = c("IN-CTX-CGE", "IN-CTX-CGE2", "IN-CTX-MGE2", "IN-CTX-MGE1"
    , "IN-STR"
    )
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = nowakowski_ex_DF, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = nowakowski_ex_DF, ex_cutoff = 5, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART16novel_NowkowskiCellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART 16 novel classified by cell type expression in human fetal brain"
    , "\nNowkowski human fetal brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART16novel_NowkowskiCellTypes_Piechart.pdf"))


## iHART 69 pie chart with Lake adult dataset
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.69" == 1]
# Check expression levels of iHART genes
subset_ex_M <- lake_ex_DF[
  rownames(lake_ex_DF) %in% genes
  , ]
gg_DF <- subset_ex_M
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF <- gg_DF[! is.na(gg_DF$Gene), ]
ggplot(gg_DF, aes(x = Gene, y = value)) +
  geom_jitter(size = 0.01) +
ggsave(paste0(outGraph, "iHART69_Lake_JitterPlot.png")
  , width = 7, height = 7, dpi = 150)
# DE filter higher
subset_DE_DF <- lake_DE_DF[
  lake_DE_DF$"Average.Difference..log.fold.change." > 0.25
  , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c("End", "Per", "Ast", "Ast_Cer", "Oli", "OPC"
    , "OPC_Cer", "Mic"
    )
  , "NA" = c("Gran"
    )
  , "Glutamatergic" = c("Ex1", "Ex2", "Ex3a", "Ex3b", "Ex3c", "Ex3d", "Ex3e"
    , "Ex4", "Ex5a", "Ex5b", "Ex6a", "Ex6b", "Ex8"
    )
  , "GABAergic" = c("In1a", "In1b", "In1c", "In2", "In3", "In4a", "In4b"
    , "In6a", "In6b", "In7", "In8", "Purk1", "Purk2"
    )
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = lake_ex_DF, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = lake_ex_DF, ex_cutoff = 1, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART69_LakeCellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART 69 classified by cell type expression in human adult brain"
    , "\nLake human adult brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART69_LakeCellTypes_Piechart.pdf"))



## iHART 16 novel pie chart with Lake adult dataset
# iHART genes
genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.17novel" == 1]
genes <- genes[! genes %in% "CACNA2D3"]
# DE filter higher
subset_DE_DF <- lake_DE_DF[
  lake_DE_DF$"Average.Difference..log.fold.change." > 0.25
  , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Glia or support cells" = c("End", "Per", "Ast", "Ast_Cer", "Oli", "OPC"
    , "OPC_Cer", "Mic"
    )
  , "NA" = c("Gran"
    )
  , "Glutamatergic" = c("Ex1", "Ex2", "Ex3a", "Ex3b", "Ex3c", "Ex3d", "Ex3e"
    , "Ex4", "Ex5a", "Ex5b", "Ex6a", "Ex6b", "Ex8"
    )
  , "GABAergic" = c("In1a", "In1b", "In1c", "In2", "In3", "In4a", "In4b"
    , "In6a", "In6b", "In7", "In8", "Purk1", "Purk2"
    )
)
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# Glutamatergic enrichment of ASD genes
enrich_log2_OR_pval_L <- Calculate_Enrichment_Log2_Odds_Ratio(
  exM = lake_ex_DF, de_DF = subset_DE_DF, genes = genes)
# Number of genes per cell class
number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
  subset_DE_DF, exM = lake_ex_DF, ex_cutoff = 1, genes = genes)
# Percent
percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
# Output table of percents
write.csv(percent_classes_DF
  , file = paste0(outTable, "iHART16novel_LakeCellTypes_Piechart.csv")
)
# Plot
ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggplot_set_theme_publication_nolabels +
  ggtitle(paste0(graphCodeTitle
    , "\n\niHART 16 novel classified by cell type expression in human adult brain"
    , "\nLake human adult brain dataset"
    , "\nGlutamatergic enrichment:"
    , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
    , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
ggsave(paste0(outGraph, "iHART16novel_LakeCellTypes_Piechart.pdf"))


Classify_Cells_By_Type <- function(subset_DE_DF, class_cluster_idx){
  print("Classify_Cells_By_Type")
  # Example class_cluster_idx input
  # $`Glia or support cells`
  # [1] "End"     "Per"     "Ast"     "Ast_Cer" "Oli"     "OPC"     "OPC_Cer"
  # [8] "Mic"
  # $`NA`
  # [1] "Gran"
  class_cluster_idx <- melt(class_cluster_idx)
  subset_DE_DF$Class <- class_cluster_idx$L1[
    match(subset_DE_DF$Cluster, class_cluster_idx$value)]
  return(subset_DE_DF)
}

Calculate_Enrichment_Log2_Odds_Ratio <- function(
  exM, de_DF, genes
  ){
    print("Calculate_Enrichment_Log2_Odds_Ratio")
    # Setup gene 0 1 matrix for glm
    gene_M <- matrix(0, nrow(exM), 2)
    rownames(gene_M) <- rownames(exM)
    idx <- rownames(gene_M) %in% genes
    gene_M[idx, 1] <- 1
    ex_neuron_genes <- subset_DE_DF$Gene[subset_DE_DF$Class == "Glutamatergic"]
    idx <- rownames(gene_M) %in% ex_neuron_genes
    gene_M[idx, 2] <- 1
    # Calculate enrichment with GLM
    dat1 <- as.numeric(gene_M[,1])
    dat2 <- as.numeric(gene_M[,2])
    glm.out <- glm(dat1~dat2,family=binomial)
    enrichP <- summary(glm.out)$coefficients[2,4]
    enrichOR <- summary(glm.out)$coefficients[2,1]
    # log2 odds ratio
    enrich_log2_OR <- log2(exp(enrichOR))
    enrich_log2_OR_pval_L <- list(
      Enrichment_Log2_Odds_Ratio = enrich_log2_OR
      , Enrichment_Pvalue = enrichP
    )
    return(enrich_log2_OR_pval_L)
}

Sum_Classes <- function(classes_to_find, class_L){
  print("Sum_Classes")
  sum(
    sapply(class_L, function(classes_list){
      all(classes_list %in% classes_to_find) & all(classes_to_find %in% classes_list)
    })
  )
}

Calculate_Number_Of_Genes_Per_Cell_Class <- function(
  subset_DE_DF, exM, ex_cutoff, genes){
    # browser()
  print("Calculate_Number_Of_Genes_Per_Cell_Class")
  # Subset DE to iHART
  subset_DE_DF <- subset_DE_DF[
    subset_DE_DF$Gene %in% genes, ]
  # Number of genes with no cluster enrichment that are broadly expressed
  subset_ex_M <- exM[
    rownames(exM) %in% genes[! genes %in% subset_DE_DF$Gene]
    , ]
  number_broadly_expressed <- sum(rowMeans(subset_ex_M) > ex_cutoff)
  # Genes that are NA category
  # (for counting genes not detected in cell types of interest)
  na_class_genes <- subset_DE_DF$Gene[subset_DE_DF$Class == "NA"]
  na_class_genes <- na_class_genes[
    ! na_class_genes %in% subset_DE_DF$Gene[subset_DE_DF$Class != "NA"]]
  # Now remove NA class genes
  subset_DE_DF <- subset_DE_DF[
    subset_DE_DF$Class != "NA", ]
  subset_DE_DF <- subset_DE_DF[
    ! is.na(subset_DE_DF$Class), ]
  class_L <- lapply(subset_DE_DF$Gene, function(Gene){
    class <- subset_DE_DF$Class[subset_DE_DF$Gene == Gene]
    class[! duplicated(class)]
  })
  names(class_L) <- subset_DE_DF$Gene
  class_L <- class_L[! duplicated(names(class_L))]
  # Quantify genes by class
  number_classes_DF <- data.frame(c(
    Glutamatergic = Sum_Classes(class_L = class_L
      , classes_to_find = "Glutamatergic")
    , GABAergic = Sum_Classes(class_L = class_L
      , classes_to_find = "GABAergic")
    , Glia_or_support_cells = Sum_Classes(class_L = class_L
      , classes_to_find = "Glia or support cells")
    , Neuron = Sum_Classes(class_L = class_L
      , classes_to_find = c("Glutamatergic", "GABAergic")
      )
    , Broadly_Expressed = sum(
        Sum_Classes(class_L = class_L
          , classes_to_find = c("Glutamatergic", "Glia or support cells")
        )
        , Sum_Classes(class_L = class_L
          , classes_to_find = c("GABAergic", "Glia or support cells")
        )
        , Sum_Classes(class_L = class_L
          , classes_to_find = c(
            "Glutamatergic", "GABAergic", "Glia or support cells")
        )
      ) + number_broadly_expressed
    , Not_Detected = sum(! genes %in% subset_DE_DF$Gene)
      - number_broadly_expressed
  ))
  return(number_classes_DF)
}

Percent_Of_Genes_Per_Cell_Class <- function(number_classes_DF){
  print("Percent_Of_Genes_Per_Cell_Class")
  percent_classes_DF <- data.frame(
    Percent = round(number_classes_DF[,1] / sum(number_classes_DF[,1]) * 100, 1)
    , Class = rownames(number_classes_DF)
  )
}
