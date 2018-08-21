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

# Lake 2017 cell annotations
lake_cell_anno_DF <- read.csv(
  "../lake_2017/Lake_2017_TS2_Cell_Annotations.csv", header = TRUE
)

# Epilepsy risk genes
epilepsy_DF <- read.table(
  "../source/Gene_Lists/High-Confidence_Epilepsy_Risk_Genes_Ruzzo_2018-05-11.txt", fill = TRUE, header = TRUE, sep = "\t"
)

# TADA Sanders 2015 = from Luis metaMat
tadaDF <- read.csv(
  "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv")
# iHART from Ruzzo et al
ihartDF <- read.csv("../source/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv")

# ID risk genes
# de novo ID: NEJM + Lancet
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_deLigt_NEJM.txt", header=T, sep="\t", fill=T)
nejm <- as.character(listMat$Gene[listMat$nature_of_mutation=="D"])
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_Rauch_Lancet.txt", header=T, sep="\t", fill=T)
lancet <- as.character(listMat$Gene_symbol[
  listMat$Type %in% c("frameshift", "nonsense", "splice")])
id_genes <- as.character(unique(c(lancet, nejm)))

## Variables
graphCodeTitle <- "Epilepsy_Expression.R"
outGraph <- "../analysis/graphs/Epilepsy_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Epilepsy_Expression_"
outTable <- "../analysis/tables/Epilepsy_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Epilepsy_Expression_"

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
  geneGroupDF <- data.frame(Gene = gene_pos_table$gene, Group = "")
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
  geneGroupDF <- data.frame(Cluster = gene_pos_table$gene, Group = "")
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
################################################################################

### Format

## Dropseq DE table
deDF <- deDF[deDF$Cluster != 16, ]

## Lake
lake_DE_DF <- lake_DE_DF[ ,1:7]
# Subset to clusters in frontal cortex dataset
lake_DE_DF <- lake_DE_DF[
  lake_DE_DF$Cluster %in% unique(gsub("_.*", "", colnames(lake_ex_DF)))
  , ]

## Nowakowski
rownames(nowakowski_ex_DF) <- nowakowski_ex_DF$geneId
nowakowski_ex_DF <- nowakowski_ex_DF[ ,-1]

## TADA Sanders_TADA0.1_exomedel
tadaDF <- tadaDF[,1:21]
tada <- unique(tadaDF[tadaDF$tadaFdrAscSscExomeSscAgpSmallDel<0.1, "RefSeqGeneName"])
################################################################################

### Data processing

# Read depth normalize Lake 2017 (ln(transcripts-per-10,000 +1))
for(i in 1:ncol(lake_ex_DF)){
  read_depth <- sum(lake_ex_DF[ ,i])
  lake_ex_DF[ ,i] <- lake_ex_DF[ ,i] / read_depth
}
lake_ex_DF <- (lake_ex_DF*10000)+1
lake_ex_DF <- log(lake_ex_DF)

# Combine Lake and drop-seq expression matrices
lake_dropseq_ex_DF <- merge(lake_ex_DF, noCentExM
  , by = "row.names", all = TRUE)

# Combine Lake and drop-seq expression matrices
# (centered and scaled separately)
cent_lake_dropseq_ex_DF <- merge(t(scale(t(lake_ex_DF))), centSO@scale.data
  , by = "row.names", all = TRUE)
rownames(cent_lake_dropseq_ex_DF) <- cent_lake_dropseq_ex_DF[ ,1]
cent_lake_dropseq_ex_DF <- cent_lake_dropseq_ex_DF[ ,-1]
################################################################################

### Heatmaps

# Heatmap + hclust
pg <- Seurat_Heatmap_By_Cluster_Hclust_Genes(genes = epilepsy_DF$Gene
  , clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
# Title
title = paste0(graphCodeTitle
  , "\n\nExpression of Epilepsy Ruzzo genes"
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
  outGraph, "Epilepsy_HeatmapDend_NormalizedCenteredScaled.png")
  , width = 12, height = 16, limitsize = FALSE)


## Heatmap sorted by cell type enrichment

# Drop-seq
genes_DF <- deDF[! deDF$Cluster %in% 16, ]
# Take cluster with highest enrichment for each gene
genes_DF <- genes_DF[order(genes_DF$Gene, -genes_DF$Log2_Fold_Change), ]
genes_DF <- genes_DF[! duplicated(genes_DF$Gene), ]
# Order by cluster and enrichment
genes_DF$Cluster <- factor(genes_DF$Cluster,
   levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
genes_DF <- genes_DF[order(genes_DF$Cluster, -genes_DF$Log2_Fold_Change), ]
# Subset to Epilepsy
genes_DF <- genes_DF[genes_DF$Gene %in% epilepsy_DF$Gene, ]
genes <- rev(genes_DF$Gene)
genes <- c(as.character(epilepsy_DF$Gene)[! epilepsy_DF$Gene %in% genes], as.character(genes))
# # Remove genes not in dataset
# genes <- genes[genes %in% rownames(centSO@scale.data)]
genes <- genes[! is.na(genes)]
# Plot
geneGroupDF <- data.frame(
  Gene = genes
  , Group = ""
)
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of Epilepsy Ruzzo genes in human fetal cortex"
    , "\nSorted by cluster enrichment"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "Epilepsy_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 7, height = 16
)

# Lake
cellID_clusterID <- gsub("_.*", "", colnames(lake_ex_DF))
names(cellID_clusterID) <- colnames(lake_ex_DF)
clusterOrder <- c("Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
  , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
  , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = t(scale(t(lake_ex_DF)))
  , cellID_clusterID = cellID_clusterID
  , clusterOrder = clusterOrder
  , clusters = clusterOrder
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of Epilepsy Ruzzo genes in human adult cortex"
    , "\nSorted by cluster enrichment"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "Epilepsy_Lake_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 11, height = 16
)

# Dropseq and Lake
cellID_clusterID <- c(gsub("_.*", "", colnames(lake_ex_DF))
  , as.numeric(as.character(centSO@ident)))
names(cellID_clusterID) <- c(colnames(lake_ex_DF), names(centSO@ident))
clusterOrder <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15
  , "Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
  , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
  , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = cent_lake_dropseq_ex_DF
  , cellID_clusterID = cellID_clusterID
  , clusterOrder = clusterOrder
  , clusters = clusterOrder
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of Epilepsy Ruzzo genes"
    , "\nHuman fetal and adult sorted by cluster enrichment"
    , "\nLake et al. frontal cortex"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "Epilepsy_LakeDropseq_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 16, height = 18
)


## Subset to high confidence epilepsy genes

# Subset
genes_highconf <- epilepsy_DF$Gene[
  epilepsy_DF$Classification == "High-confidence"]
geneGroupDF <- geneGroupDF[geneGroupDF$Gene %in% genes_highconf, ]
# Write out gene list to copy paste into illustrator
write.table(rev(geneGroupDF$Gene)
  , file = paste0(outTable, "EpilepsyHighConf.txt")
  , quote = FALSE, row.names = FALSE)
# Plot
# Dropseq and Lake
cellID_clusterID <- c(gsub("_.*", "", colnames(lake_ex_DF))
  , as.numeric(as.character(centSO@ident)))
names(cellID_clusterID) <- c(colnames(lake_ex_DF), names(centSO@ident))
clusterOrder <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15
  , "Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
  , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
  , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = cent_lake_dropseq_ex_DF
  , cellID_clusterID = cellID_clusterID
  , clusterOrder = clusterOrder
  , clusters = clusterOrder
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of Epilepsy Ruzzo high confidence genes"
    , "\nHuman fetal and adult sorted by cluster enrichment"
    , "\nLake et al. frontal cortex"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "EpilepsyHighConf_LakeDropseq_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 16, height = 22
)

## Intersection ASD, epilepsy, ID

# Initialize matrix with high conf epilepsy genes
genes_highconf <- epilepsy_DF$Gene[
  epilepsy_DF$Classification == "High-confidence"]
geneGroupDF <- geneGroupDF[geneGroupDF$Gene %in% genes_highconf, ]
gene_binary_M <- matrix(NA, length(geneGroupDF$Gene), 0)
rownames(gene_binary_M) <- geneGroupDF$Gene
# Add genes lists
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M, tada, "ASD")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M, geneGroupDF$Gene, "Epilepsy")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M, id_genes, "ID")
# Format
gg_DF <- melt(gene_binary_M)
gg_DF$Color <- with(gg_DF
  , ifelse(Var2 == "ASD" & value == 1, "ASD"
  , ifelse(Var2 == "Epilepsy" & value == 1, "Epilepsy"
  , ifelse(Var2 == "ID" & value == 1, "ID"
    , NA)))
)
# Plot
ggplot(gg_DF, aes(x = Var2, y = Var1, fill = Color)) +
  geom_tile(width=0.7, height=0.7) +
  scale_fill_discrete(na.value = "lightgrey") +
  ggplot_set_theme_publication +
  ggtitle(paste0(
    graphCodeTitle
      , "\n\nIntersection of ASD, epilepsy, ID gene lists"
      , "\n")
  )
ggsave(paste0(
  outGraph, "EpilepsyHighConf_EnrichmentSorted_DiseaseIntersection.pdf")
  , width = 4, height = 12)
################################################################################
