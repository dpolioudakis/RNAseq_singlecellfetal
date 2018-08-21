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
require(ggpubr)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

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
  "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/AllDeNovoByGene.csv")

# de Rubeis ASD - from Luis metaMat
rubDF <- read.table(
  "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/deRubeis_mutationlist.txt"
  , header = TRUE, sep = "\t")

# TADA Sanders 2015 = from Luis metaMat
tadaDF <- read.csv(
  "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv")

ihartDF <- read.csv("../source/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv")

# Epilepsy risk genes
epilepsy_DF <- read.table(
  "../source/Gene_Lists/High-Confidence_Epilepsy_Risk_Genes_Ruzzo_2018-05-11.txt", fill = TRUE, header = TRUE, sep = "\t"
)

# ID risk genes
# de novo ID: NEJM + Lancet
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_deLigt_NEJM.txt", header=T, sep="\t", fill=T)
nejm <- as.character(listMat$Gene[listMat$nature_of_mutation=="D"])
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_Rauch_Lancet.txt", header=T, sep="\t", fill=T)
lancet <- as.character(listMat$Gene_symbol[
  listMat$Type %in% c("frameshift", "nonsense", "splice")])
id_genes <- as.character(unique(c(lancet, nejm)))

# Known marker Luis table
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv")

# Biomart gene info
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# BrainSpan developmental transcriptome (FPKM)
bsDF <- read.csv(
  "../source/BrainSpan_DevTranscriptome/genes_matrix_csv/expression_matrix.csv", header = FALSE)
bs_rnames <- read.csv(
  "../source/BrainSpan_DevTranscriptome/genes_matrix_csv/rows_metadata.csv")
bsMtDF <- read.csv(
  "../source/BrainSpan_DevTranscriptome/genes_matrix_csv/columns_metadata.csv")

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

## BrainSpan
Calculate_BrainSpan_Expression <- function(genes){
  # Subset to cortex, remove visual and cerebellar
  subset_BsMt_DF <- bsMtDF[grep("cortex", bsMtDF$structure_name), ]
  subset_BsMt_DF <- subset_BsMt_DF[
    subset_BsMt_DF$structure_name != "cerebellar cortex" &
    subset_BsMt_DF$structure_name != "primary visual cortex (striate cortex, area V1/17)", ]
  ssBsDF <- bsDF[ ,subset_BsMt_DF$column_num]
  bsDF[bs_rnames$gene_symbol == "TCF7L2", ]
  # Log2 transform
  ssBsDF <- log(ssBsDF, 2)
  # Subset to genes of interest
  idx <- match(genes, bs_rnames$gene_symbol)
  idx <- idx[! is.na(idx)]
  df1 <- bs_rnames[idx, ]
  ssBsDF <- ssBsDF[df1$row_num, ]
  row.names(ssBsDF) <- df1$gene_symbol
  # Genes as columns
  ssBsDF <- as.data.frame(t(ssBsDF))
  # Remove outliers by stdev
  ssBsDF <- as.data.frame(apply(ssBsDF, 2, function(x) {
    Remove_Outliers_By_SD(x, nStdev = 2.5)
    }))
  # Add age
  ssBsDF$Age <- factor(subset_BsMt_DF$age, levels = unique(subset_BsMt_DF$age))
  # Add donor ID
  ssBsDF$Donor <- factor(subset_BsMt_DF$donor_id
    , levels = unique(subset_BsMt_DF$donor_id))
  # Add structure ID
  ssBsDF$Structure <- factor(subset_BsMt_DF$structure_acronym
    , levels = unique(subset_BsMt_DF$structure_acronym))
  # Format
  ssBsDF <- melt(ssBsDF)
  return(ssBsDF)
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
geneGroupDF <- data.frame(Cluster = df1$Gene, Group = df1$Cluster)
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
geneGroupDF <- data.frame(Cluster = genes, Group = "")
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

### ASD TADA featureplots and heatmaps

# ASD TADA genes
genes <- tada

# Check expression levels of TADA genes
subset_ex_M <- lake_ex_DF[
  rownames(lake_ex_DF) %in% genes
  , ]
gg_DF <- subset_ex_M
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF <- gg_DF[! is.na(gg_DF$Gene), ]
ggplot(gg_DF, aes(x = Gene, y = value)) +
  geom_jitter(size = 0.01) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste0(outGraph, "TADA_Lake_JitterPlot.png")
  , width = 7, height = 7, dpi = 150)

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
geneGroupDF <- data.frame(Cluster = df1$Gene, Group = df1$Cluster)
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
################################################################################

### ASD TADA sorted by enrichment plots

## Heatmap of ASD TADA genes with fixed column widths

# Drop-seq
genes_DF <- deDF[! deDF$Cluster %in% 16, ]
# Take cluster with highest enrichment for each gene
genes_DF <- genes_DF[order(genes_DF$Gene, -genes_DF$Log2_Fold_Change), ]
genes_DF <- genes_DF[! duplicated(genes_DF$Gene), ]
# Order by cluster and enrichment
genes_DF$Cluster <- factor(genes_DF$Cluster,
   levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
genes_DF <- genes_DF[order(genes_DF$Cluster, -genes_DF$Log2_Fold_Change), ]
# Subset to TADA
genes_DF <- genes_DF[genes_DF$Gene %in% tada, ]
genes <- rev(genes_DF$Gene)
genes <- c(as.character(tada)[! tada %in% genes], as.character(genes))
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
    , "\n\nExpression of ASD TADA genes sorted by cluster enrichment"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
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
    , "\n\nExpression of ASD TADA genes in adult sorted by cluster enrichment"
    , "\nLake et al. frontal cortex and drop-seq fetal centered and scaled together"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_Lake_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
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
    , "\n\nExpression of ASD TADA genes"
    , "\nHuman fetal and adult sorted by cluster enrichment"
    , "\nLake et al. frontal cortex"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
  , "TADA_LakeDropseq_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 16, height = 18
)

# Dropseq and Lake centered and scaled together
cellID_clusterID <- c(gsub("_.*", "", colnames(lake_ex_DF)), centSO@ident)
names(cellID_clusterID) <- c(colnames(lake_ex_DF), names(centSO@ident))
clusterOrder <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15
  , "Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
  , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
  , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = t(scale(t(lake_dropseq_ex_DF)))
  , cellID_clusterID = cellID_clusterID
  , clusterOrder = clusterOrder
  , clusters = clusterOrder
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes in human adult sorted by cluster enrichment"
    , "\nLake et al. frontal cortex"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_LakeDropseq_EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 16, height = 12
)

## Intersection ASD, epilepsy, ID

# Initialize matrix
gene_binary_M <- matrix(NA, length(geneGroupDF$Gene), 0)
rownames(gene_binary_M) <- geneGroupDF$Gene
# Add genes lists
epilepsy_genes <- epilepsy_DF$Gene[
  epilepsy_DF$Classification == "High-confidence"]
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M, geneGroupDF$Gene, "ASD")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M, epilepsy_genes, "Epilepsy")
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
ggsave(paste0(outGraph, "TADA_EnrichmentSorted_DiseaseIntersection.pdf")
  , width = 4, height = 8)

## Check characteristics of genes not cluster enriched in dropseq

# CDS length and GC content
idx <- match(geneGroupDF$Gene, bmDF$hgnc_symbol)
gg_DF <- bmDF[idx, ]
gg_DF$Enriched <- gg_DF$hgnc_symbol %in% deDF$Gene
aggregate(cds_length~Enriched, gg_DF, mean)
aggregate(percentage_gc_content~Enriched, gg_DF, mean)
ggplot(gg_DF, aes(x = Enriched, y = cds_length, fill = Enriched)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,16000)) +
  stat_compare_means(cds_length ~ Enriched, data = gg_DF
    , comparisons = list(c(1,2))
    , method = "t.test", p.adjust.method = "none", label = "p.signif"
    # , label.y = c(limits_L$Maximum*c(1.2, 1.4, 1.6))
  ) +
  ylab("CDS length") +
  ggplot_set_theme_publication +
  ggtitle(paste0(
    graphCodeTitle
      , "\n\nCDS length of ASD TADA genes not enriched in dropseq"
      , "\n")
  )
ggsave(paste0(outGraph, "TADA_ClusterEnriched_CDSlength.pdf")
  , width = 5, height = 5)
ggplot(gg_DF, aes(x = Enriched, y = percentage_gc_content, fill = Enriched)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,100)) +
  stat_compare_means(cds_length ~ Enriched, data = gg_DF
    , comparisons = list(c(1,2))
    , method = "t.test", p.adjust.method = "none", label = "p.signif"
    # , label.y = c(limits_L$Maximum*c(1.2, 1.4, 1.6))
  ) +
  ylab("GC content") +
  ggplot_set_theme_publication +
  ggtitle(paste0(
    graphCodeTitle
      , "\n\nPercent GC content of ASD TADA genes not enriched in dropseq"
      , "\n")
  )
ggsave(paste0(outGraph, "TADA_ClusterEnriched_percentGC.pdf")
  , width = 5, height = 5)

# Expression level in brainspan
ssBsDF <- Calculate_BrainSpan_Expression(gg_DF$hgnc_symbol)
ssBsDF$Enriched <- gg_DF$Enriched[match(ssBsDF$variable, gg_DF$hgnc_symbol)]
ssBsDF <- ssBsDF[ssBsDF$value != "-Inf", ]
ssBsDF <- ssBsDF[! is.na(ssBsDF$Enriched), ]
ssBsDF <- ssBsDF[! is.na(ssBsDF$value), ]
ssBsDF$Pre_Vs_Postnatal <- "Postnatal"
ssBsDF$Pre_Vs_Postnatal[grepl("pcw", ssBsDF$Age)] <- "Prenatal"
ssBsDF$Pre_Vs_Postnatal <- factor(
  ssBsDF$Pre_Vs_Postnatal, levels = c("Prenatal", "Postnatal")
)
ggplot(ssBsDF, aes(x = Enriched, y = value, fill = Enriched)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-3,10)) +
  facet_wrap(~Pre_Vs_Postnatal) +
  stat_compare_means(value ~ Enriched, data = ssBsDF
    , comparisons = list(c(1,2))
    , method = "t.test", p.adjust.method = "none", label = "p.signif"
    , label.y = c(9)
  ) +
  ylab("Normalized expression (log2 FPKM)") +
  ggplot_set_theme_publication +
  ggtitle(paste0(
    graphCodeTitle
      , "\n\nExpression of ASD TADA genes in brainspan not enriched in dropseq"
      , "\n")
  )
ggsave(paste0(outGraph, "TADA_ClusterEnriched_BrainSpan.pdf")
  , width = 5, height = 5)
aggregate(value~Enriched, ssBsDF, mean)
# Plot as fit line
# Mean of expression of genes
ssBsDF <- aggregate(value~Age+Donor+Structure+Enriched, ssBsDF, mean)
# ggplot fit line
ggplot(ssBsDF, aes(x = Age, y = value, group = Enriched, color = Enriched)) +
  # facet_wrap(~Enriched, scales = "free_x", ncol = 2) +
  geom_jitter(size = 0.1, width = 0.2) +
  geom_smooth() +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  # ylim(c(0,5)) +
  xlab("Age") +
  ylab("Normalized expression (log2 FPKM)") +
  ggtitle(paste0(
    graphCodeTitle
      , "\n\nExpression of ASD TADA genes in brainspan not enriched in dropseq"
      , "\n")
  )
ggsave(paste0(outGraph, "TADA_ClusterEnriched_BrainSpan_Fit.png")
  , width = 5, height = 5)

## TADA neuron enriched

# Subset to genes enriched in neuron clusters
genes_DF <- deDF[deDF$Cluster %in% c(0,1,3,4,13,5,6) & deDF$Gene %in% tada, ]
# Take cluster with highest enrichment for each gene
genes_DF <- genes_DF[order(genes_DF$Gene, -genes_DF$Log2_Fold_Change), ]
genes_DF <- genes_DF[! duplicated(genes_DF$Gene), ]
# Order by cluster and enrichment
genes_DF$Cluster <- factor(genes_DF$Cluster,
   levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
genes_DF <- genes_DF[order(genes_DF$Cluster, -genes_DF$Log2_Fold_Change), ]
genes <- rev(genes_DF$Gene)
geneGroupDF <- data.frame(
  Gene = genes
  , Group = ""
)
cellID_clusterID <- centSO@ident
# Write out gene list to copy paste into illustrator
write.table(rev(geneGroupDF$Gene)
  , file = paste0(outTable, "TADA_NeuronEnriched.txt")
  , quote = FALSE, row.names = FALSE)
# Plot
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , cellID_clusterID = cellID_clusterID)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_NeuronEnriched_HeatmapSetColWidths_NormalizedCenteredScaled_paper.png"
  )
  , width = 12, height = 6
)

## TADA progenitor, OPC, and pericyte enriched genes
geneGroupDF <- data.frame(
  Gene = rev(c("ILF2", "TRIO", "SETD5", "TCF7L2", "KAT2B", "SLC6A1"))
  , Group = ""
)
# Dropseq
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA progenitor OPC and pericyte enriched genes"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_ProgenitorOPCandPericyte_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 12, height = 4
)

## TADA genes of interest
geneGroupDF <- data.frame(Gene = rev(c("ILF2", "BCL11A", "CTTNBP2", "DIP2A", "CUL3", "SLC6A1", "PTEN", "MBD5A", "DSCAM")), Group = "")
# Dropseq
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes of interest"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_GenesOfInterest_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 12, height = 6
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
    , "\n\nExpression of ASD TADA genes of interest in human adult"
    , "\nLake et al. frontal cortex"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_Lake_GenesOfInterest_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 12, height = 6
)
################################################################################

### ASD TADA unique to dropseq or lake et al

mean_exp <- rowMeans(noCentExM[rownames(noCentExM) %in% tada, ])
genes_1 <- names(mean_exp)[mean_exp < 0.2]
genes_2 <- tada[! tada %in% deDF$Gene]
genes_not_dropseq <- intersect(genes_1, genes_2)

mean_exp <- rowMeans(lake_ex_DF[rownames(lake_ex_DF) %in% tada, ])
genes_1 <- names(mean_exp)[mean_exp >= 0.5]
genes_2 <- tada[tada %in% lake_DE_DF$Gene]
genes_in_lake <- union(genes_1, genes_2)

genes_unique_lake <- intersect(genes_not_dropseq, genes_in_lake)

mean_exp <- rowMeans(noCentExM[rownames(noCentExM) %in% tada, ])
genes_1 <- names(mean_exp)[mean_exp >= 0.5]
genes_2 <- tada[tada %in% deDF$Gene]
genes_in_dropseq <- union(genes_1, genes_2)

mean_exp <- rowMeans(lake_ex_DF[rownames(lake_ex_DF) %in% tada, ])
genes_1 <- names(mean_exp)[mean_exp < 0.2]
genes_2 <- tada[! tada %in% lake_DE_DF$Gene]
genes_not_lake <- intersect(genes_1, genes_2)

genes_unique_dropseq <- intersect(genes_in_dropseq, genes_not_lake)

## TADA genes in dropseq not in Lake et al
# Plot
geneGroupDF <- data.frame(
  Gene = genes_unique_dropseq
  , Group = ""
)
# Dropseq
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes"
    , "\nGenes detected in fetal drop-seq and not in human adult (Lake et al.)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_UniqueDropseq_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 8, height = 12
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
    , "\n\nExpression of ASD TADA genes adult"
    , "\nGenes detected in fetal drop-seq and not in human adult (Lake et al.)"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_UniqueDropseq_Lake_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 12, height = 12
)

## TADA genes in Lake et al not in dropseq
# Plot
geneGroupDF <- data.frame(
  Gene = genes_unique_lake
  , Group = ""
)
# Dropseq
cellID_clusterID <- centSO@ident
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = geneGroupDF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of ASD TADA genes"
    , "\nGenes detected in human adult (Lake et al.) not in fetal dropseq"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_UniqueLake_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 8, height = 12
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
    , "\n\nExpression of ASD TADA genes in adult"
    , "\nGenes detected in human adult (Lake et al.) not in fetal dropseq"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "TADA_UniqueLake_Lake_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 12, height = 12
)

ssBsDF <- Calculate_BrainSpan_Expression(genes_unique_lake)
# Mean of expression of genes
# ssBsDF <- aggregate(value~Donor+Age+Structure, ssBsDF, mean)
# ggplot fit line
ggplot(ssBsDF, aes(x = Age, y = value, group = 1)) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  geom_jitter(size = 0.1, width = 0.2) +
  geom_smooth(color = "red") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  # ylim(c(0,5)) +
  xlab("Age") +
  ylab("Normalized expression")
ggsave(paste0(outGraph, "TADA_Lake_Unique_BrainSpan_Fit_paper.pdf")
  , width = 12, height = 12)

ssBsDF <- Calculate_BrainSpan_Expression(genes_unique_dropseq)
# ggplot fit line
ggplot(ssBsDF, aes(x = Age, y = value, group = 1)) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  geom_jitter(size = 0.1, width = 0.2) +
  geom_smooth(color = "red") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  # ylim(c(0,5)) +
  xlab("Age") +
  ylab("Normalized expression")
ggsave(paste0(outGraph, "TADA_Dropseq_Unique_BrainSpan_Fit_paper.pdf")
  , width = 12, height = 12)

ssBsDF <- Calculate_BrainSpan_Expression(c("ILF2", "BCL11A", "CTTNBP2", "DIP2A", "CUL3", "SLC6A1", "PTEN", "MBD5A", "DSCAM"))
# ggplot fit line
ggplot(ssBsDF, aes(x = Age, y = value, group = 1)) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  geom_jitter(size = 0.1, width = 0.2) +
  geom_smooth(color = "red") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  # ylim(c(0,5)) +
  xlab("Age") +
  ylab("Normalized expression")
ggsave(paste0(outGraph, "TADA_GenesOfInterest_BrainSpan_Fit_paper.pdf")
  , width = 12, height = 12)
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
geneGroupDF <- data.frame(Cluster = genes, Group = "")
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
geneGroupDF <- data.frame(Gene = genes, Group = "")
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
################################################################################
