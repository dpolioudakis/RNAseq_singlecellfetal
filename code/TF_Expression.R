# Damon Polioudakis
# 2017-05-30
# TF expression

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
require(cowplot)
require(xlsx)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# Log normalized, regressed nUMI and percent mito
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
centSO <- ssCentSO
noCentExM <- ssNoCentExM

# BrainSpan developmental transcriptome
bsDF <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/expression_matrix.csv", header = FALSE)
rnames <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/rows_metadata.csv")
row.names(bsDF) <- rnames$gene_symbol
bsMtDF <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/columns_metadata.csv")

# Miller
load("../neurogenesis/orig.data/LCMDE/AllenLCM.Rdata")
MillerExprRAW = AllenLCM$datExpr
MillerMetaRAW = AllenLCM$datTraits
Zones = read.csv("../neurogenesis/orig.data/LCMDE/LCM_Zones_CPio.csv")
MillerAnnotRAW = read.csv("../neurogenesis/orig.data/LCMDE/annot.csv", row.names = 1)

# Marker gene lists
# deDF <- read.table(
#   "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , header = TRUE)
deDF <- read.table(
  "../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Human TFs
tfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_TF_EnsemblID.txt")

# Human chromatin remodeling factors
crDF <- read.table(
  "../source/AnimalTFDB_Homo_sapiens_chr_remodeling_factor_EnsemblID.txt")

# Human co-factors
cfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_cofactor_EnsemblID.txt")

# biomaRt gene info
# bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh37_Ensembl75.csv"
#   , header = TRUE)
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# He 2007 Supp Table 2 layer expression
heDF <- read.csv("../source/He_NatNeurosci_2017_ST2.csv", header = TRUE)

# Allen Developmental Macaque human specific genes
hsDF <- read.csv("../source/Bakken_2016_AllenDevMacaque_ST10_HumanSpecific.csv"
  , header = TRUE, skip = 1)

# Luis ATAC human enhancer genes
aheCpDF <- read.table("../source/ATAC-HiC-CombProxDist-CP.txt", header = TRUE)
aheGzDF <- read.table("../source/ATAC-HiC-CombProxDist-VZ.txt", header = TRUE)

# DE Vijay pipeline TF enrichment
# JASPAR
# Results
tbejDF <- read.csv("../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/Seurat_ClusterDE_TFenrichment_logFC04_TF_Enrichment_Results.csv"
  , header = TRUE)
# TFs in JASPAR database
jasparDF <- read.csv("../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/Compiled_TFs_in_JASPAR_Database.csv"
  , header = TRUE)
# TRANSFAC
# Results
tbetDF <- read.csv("../analysis/Seurat_ClusterDE_TFenrichment/TRANSFAC/Seurat_ClusterDE_TFenrichment_logFC04_TF_Enrichment_Results.csv"
  , header = TRUE)
# TFs in JASPAR database
jasparDF <- read.csv("../analysis/Seurat_ClusterDE_TFenrichment/TRANSFAC/Compiled_TFs_in_TRANSFAC_Database.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "TF_Expression.R"
outGraph <- "../analysis/graphs/TF_Expression_DS2-11/TF_Expression_DS2-11_"

## Output Directories
outGraphDir <- dirname(outGraph)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

Subset_TFs <- function(deDF, cluster, okayClusters, fcHigh, fcLow) {
  # Clusters gene cannot be DE in
  clsNo <- c(0:17)[! c(0:17) %in% okayClusters]
  # Gene is > X FC in cluster
  genes1 <- deDF$GENE[deDF$LOG_FC > fcHigh & deDF$CLUSTER == cluster]
  # Genes in clusters genes cannot be DE in > 0.3
  genes2 <- deDF$GENE[deDF$LOG_FC > fcLow & deDF$CLUSTER %in% clsNo]
  # Check
  table(genes1 %in% genes2)
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

# Expression heatmap
Heatmap_By_Cluster <- function(
  ggDF, seuratO, clusters, lowerLimit, upperLimit, ggTitle, levels) {

  colnames(ggDF)[1:2] <- c("GENE", "GROUP")
  # Remove blanks
  ggDF <- ggDF[! ggDF$GENE == "", ]
  # Save order to set levels later
  levels <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  ggDF <- melt(ggDF)
  # Add seurat clusters
  idx <- match(ggDF$variable, names(seuratO@ident))
  ggDF$SEURAT_CLUSTERS <- seuratO@ident[idx]
  # Subset clusters
  ggDF <- ggDF[ggDF$SEURAT_CLUSTERS %in% clusters, ]
  ggDF$GENE_GROUP <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  ggDF$GENE_GROUP <- factor(ggDF$GENE_GROUP, levels = levels)
  # Set limits
  ggDF$value[ggDF$value < lowerLimit] <- lowerLimit
  ggDF$value[ggDF$value > upperLimit] <- upperLimit
  # ggplot
  ggplot(ggDF, aes(x = variable, y = GENE_GROUP, fill = value)) +
    geom_tile() +
    facet_grid(GROUP~SEURAT_CLUSTERS, space = "free", scales = "free") +
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
}

# Use plot_grid to combine 3 heatmaps to deal with cell number scaling
Heatmaps_Combined <- function(ggDF, seuratO, clusters1, clusters2, clusters3
  , lowerLimit, upperLimit, title) {
  p1 <- Heatmap_By_Cluster(ggDF, seuratO = centSO, clusters = clusters1
    , lowerLimit = lowerLimit, upperLimit = upperLimit
  )
  p1 <- p1 + theme(
    axis.title.x = element_blank()
    , strip.text.y = element_blank()
    , legend.position = "none"
  )
  p2 <- Heatmap_By_Cluster(ggDF, seuratO = centSO, clusters = clusters2
    , lowerLimit = lowerLimit, upperLimit = upperLimit
  )
  p2 <- p2 + theme(
    strip.text.y = element_blank()
    , legend.position = "none"
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank()
  )
  p3 <- Heatmap_By_Cluster(ggDF, seuratO = centSO, clusters = clusters3
    , lowerLimit = lowerLimit, upperLimit = upperLimit
  )
  p3 <- p3 + theme(
    axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank())
  # plot_grid combine
  pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'h', axis = 'b')
  # # now add the title
  # title <- ggdraw() + draw_label(title)
  # # rel_heights values control title margins
  # pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  return(pg)
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

### Format

## Format and subset Miller

ZonesMiller = MillerMetaRAW$structure_acronym

# Select cortical layers only
MillerMetaRAW$Labels = rep(NA, nrow(MillerMetaRAW))
for (Label in Zones$Label)
  MillerMetaRAW[grep(Label, ZonesMiller), "Labels"] = Label
Zoneindex = which(! is.na(MillerMetaRAW$Labels))
MillerExpr = MillerExprRAW[, Zoneindex]
MillerMeta = MillerMetaRAW[Zoneindex, ]

# Only ENTREZ annotated probes
MillerExpr = MillerExpr[! is.na(MillerAnnotRAW$ENTREZ_ID), ]
MillerAnnot = MillerAnnotRAW[! is.na(MillerAnnotRAW$ENTREZ_ID), ]

# Get maximum expression probe
Genes = unique(MillerAnnot$ENTREZ_ID)
keepind = matrix(nrow = 0, ncol = 0);
for (ii in 1:length(Genes)) {
  genematchind = which(MillerAnnot$ENTREZ_ID == Genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(MillerExpr[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}

MillerExpr = MillerExpr[keepind, ];
rownames(MillerExpr) <- MillerAnnot$ENTREZ_ID[
  match(rownames(MillerExpr), rownames(MillerAnnot))]
# rownames(MillerExpr) = Genes;
MillerAnnot = MillerAnnot[keepind, ];

# Miller
millerZones <- MillerMeta[match(MillerMeta$well_id, colnames(MillerExpr)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)

# ## Convert gene symbol to ensembl
# deDF$ensembl_gene_id <- bmDF$ensembl_gene_id[match(deDF$GENE, bmDF$hgnc_symbol)]

# ## Assign annotated cluster names to clusters
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
# new.cluster.ids <- c(
#   "Excitatory Upper Layer Neuron 1"
#   , "Excitatory Neuron"
#   , "Excitatory Upper Layer Neuron 2"
#   , "Excitatory Deep Layer Neuron"
#   , "Intermediate Progenitors"
#   , "Interneuron"
#   , "Mitotic Progenitors"
#   , "oRG"
#   , "Oligodendrocyte Precursor"
#   , "Endothelial")
# seuratO@ident <- plyr::mapvalues(seuratO@ident, from = current.cluster.ids
#   , to = new.cluster.ids)
# # Stash numerical cluster identities if want to use later
# seuratO <- StashIdent(seuratO, save.name = "Cluster_Numbers")
# 
# mapDF <- data.frame(VALUE = current.cluster.ids, MAPVALUE = new.cluster.ids)

current.cluster.ids <- c(sort(as.numeric(as.character(unique(centSO@ident)))))
new.cluster.ids <- c(sort(as.numeric(as.character(unique(centSO@ident)))))
mapDF <- data.frame(VALUE = current.cluster.ids, MAPVALUE = new.cluster.ids)
################################################################################

### TFs, co-factors, chromatin remodelers

## Genes unique DE for Excitatory (0, 1, 4, 13), Deep layer (3, 14)
## , Interneuron (5, 6), RG (7, 8, 9, 10)

# Data frame of TFs, co-factors, chromatin remodelers
df <- rbind(tfDF, crDF, cfDF)
# 1970
nrow(df)
# Number of DE genes 7839
nrow(deDF)
# Intersect DE gene lists and TFs, co-factors, chromatin remodelers list: 408
length(intersect(deDF$ENSEMBL, df$V1))
df <- deDF[deDF$ENSEMBL %in% df$V1, ]

# RG
# Heatmap
# Normalized, mean centering scaling
ldf <- lapply(c(7, 9), function(cluster) {
  Subset_TFs(deDF = df, cluster = cluster, okayClusters = c(7, 8, 9, 10)
    , fcHigh = 0.4, fcLow = 0.3)
})
utdeDF <- do.call("rbind", ldf)
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = utdeDF, exDF = centSO@scale.data)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed in RG clusters"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster; < 0.3 for other clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(
  outGraph, "TFsCRsCFs_DeTypeUniqueRG_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4 + nrow(ggDF)/6, limitsize = FALSE)

# Excitatory and deep layer
# Heatmap
# Normalized, mean centering scaling
ldf <- lapply(c(0, 1, 4, 13, 3, 14), function(cluster) {
  Subset_TFs(deDF = df, cluster = cluster, okayClusters = c(0, 1, 4, 13, 3, 14)
    , fcHigh = 0.4, fcLow = 0.3)
})
utdeDF <- do.call("rbind", ldf)
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = utdeDF, exDF = centSO@scale.data)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed in excitatory and deep layer clusters"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster; < 0.3 for other clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(
  outGraph, "TFsCRsCFs_DeTypeUniqueExcAndDeep_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4 + nrow(ggDF)/6, limitsize = FALSE)

# Excitatory
# Heatmap
# Normalized, mean centering scaling
ldf <- lapply(c(0, 1, 4, 13), function(cluster) {
  Subset_TFs(deDF = df, cluster = cluster, okayClusters = c(0, 1, 4, 13)
    , fcHigh = 0.4, fcLow = 0.3)
})
utdeDF <- do.call("rbind", ldf)
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = utdeDF, exDF = centSO@scale.data)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed in excitatory clusters"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster; < 0.3 for other clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(
  outGraph, "TFsCRsCFs_DeTypeUniqueExcitatory_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4 + nrow(ggDF)/6, limitsize = FALSE)

# Deep layer
# Heatmap
# Normalized, mean centering scaling
ldf <- lapply(c(3, 14), function(cluster) {
  Subset_TFs(deDF = df, cluster = cluster, okayClusters = c(3, 14)
    , fcHigh = 0.4, fcLow = 0.3)
})
utdeDF <- do.call("rbind", ldf)
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = utdeDF, exDF = centSO@scale.data)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed in excitatory deep layer clusters"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster; < 0.3 for other clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(
  outGraph, "TFsCRsCFs_DeTypeUniqueExcDeepLayer_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4 + nrow(ggDF)/6, limitsize = FALSE)

# Interneuron
# Heatmap
# Normalized, mean centering scaling
ldf <- lapply(c(5, 6), function(cluster) {
  Subset_TFs(deDF = df, cluster = cluster, okayClusters = c(5, 6)
    , fcHigh = 0.4, fcLow = 0.3)
})
utdeDF <- do.call("rbind", ldf)
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = utdeDF, exDF = centSO@scale.data)
# Plot as heatmap
ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed in interneuron clusters"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster; < 0.3 for other clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(
  outGraph, "TFsCRsCFs_DeTypeUniqueInterneuron_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4 + nrow(ggDF)/6, limitsize = FALSE)

## Unique DE genes

# Keep only DE genes that are unique to a cluster
keep <- names(table(deDF$GENE))[table(deDF$GENE) == 1]
udeDF <- deDF[deDF$GENE %in% keep, ]

df <- rbind(tfDF, crDF, cfDF)
# 1970
nrow(df)
nrow(udeDF)
length(intersect(udeDF$ENSEMBL, df$V1))
df <- udeDF[udeDF$ENSEMBL %in% df$V1, ]

# Heatmap
# Normalized, no mean centering scaling
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = df, exDF = noCentExM)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = 0, upperLimit = 3
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n"
  , "\nTFs, Chromatin remodelers, co-factors differentially expressed in a cluster"
  , "\nSeurat normalized expression"
  , "\n"
)
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
ggsave(paste0(outGraph, "TFsCRsCFs_DEunique_Heatmap_Normalized.png")
  , width = 12, height = nrow(ggDF)/6)

# Heatmap
# Normalized, mean centering scaling
# Subset expression matrix by DE genes
ggDF <- Combine_DE_and_Expression(deDF = df, exDF = centSO@scale.data)
# Plot as heatmap
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n"
  , "\nTFs, Chromatin remodelers, co-factors differentially expressed in a cluster"
  , "\nSeurat normalized expression, mean centered, variance scaled"
  , "\n"
)
ggsave(paste0(outGraph, "TFsCRsCFs_DEunique_Heatmap_NormCentScale.png")
  , width = 12, height = (nrow(ggDF)/6))

########################## BROKEN HERE
## Non unique DE genes

# Data frame of TFs, co-factors, chromatin remodelers
df <- rbind(tfDF, crDF, cfDF)
# 1970
nrow(df)
# Number of DE genes 7839
nrow(deDF)
# Number of DE genes > 0.7 log FC: 1431
df2 <- deDF[deDF$LOG_FC > 0.4, ]
nrow(df2)
# Intersect DE gene lists and TFs, co-factors, chromatin remodelers list: 408
length(intersect(df2$ENSEMBL, df$V1))
df <- df2[df2$ENSEMBL %in% df$V1, ]

# Heatmap
# Normalized, mean centering scaling
ggDF <- merge(df[c("GENE", "CLUSTER")], centSO@scale.data
  , by.x = 1, by.y = "row.names", all.x = TRUE)
ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\nDE filters: > 0.4 log fold change in cluster"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "TFsCRsCFs_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = nrow(ggDF)/6, limitsize = FALSE)

# Heatmap
# Normalized, no mean centering scaling
ggDF <- merge(df[c("GENE", "CLUSTER")], noCentExM
  , by.x = 1, by.y = "row.names", all.x = TRUE)
pg <- Heatmaps_Combined(
  ggDF = ggDF, seuratO = centSO, lowerLimit = 0, upperLimit = 3
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17), title = ""
)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of TFs, chromatin remodelers, co-factors differentially expressed"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression"
  , "\nDE filters: > 0.4 log fold change in cluster"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "TFsCRsCFs_Heatmap_Normalized.png")
  , width = 12, height = nrow(ggDF)/6, limitsize = FALSE)
################################################################################

### Brainspan developmental transcriptome

unique(bsMtDF$structure_name)
# Structures with cortex in name:
# [1] "occipital neocortex"
# [2] "primary motor-sensory cortex (samples)"
# [3] "posterior (caudal) superior temporal cortex (area 22c)"
# [4] "anterior (rostral) cingulate (medial prefrontal) cortex"
# [5] "dorsolateral prefrontal cortex"
# [6] "orbital frontal cortex"
# [7] "inferolateral temporal cortex (area TEv, area 20)"
# [8] "ventrolateral prefrontal cortex"
# [9] "parietal neocortex"
# [10] "temporal neocortex"
# [11] "primary auditory cortex (core)"
# [12] "primary visual cortex (striate cortex, area V1/17)"
# [13] "primary motor cortex (area M1, area 4)"
# [14] "posteroventral (inferior) parietal cortex"
# [15] "primary somatosensory cortex (area S1, areas 3,1,2)"
# [16] "cerebellar cortex"

# Subset to cortex, remove visual and cerebellar
df <- bsMtDF[grep("cortex", bsMtDF$structure_name), ]
df[df$structure_name != "cerebellar cortex" &
    df$structure_name != "primary visual cortex (striate cortex, area V1/17)", ]
ssBsDF <- bsDF[ ,df$column_num]
# Subset to genes of interest
idx <- match(c(
  # RG
  "VIM", "HES1", "HOPX", "ITGB5", "CARHSP1", "ZFHX4", "LITAF", "MAFF"
  # Excitatory - not deep layer
  , "TUBB3", "STMN2", "SATB2", "CUX1", "CSRP2"
  # Excitatory deep layer
  , "BCL11B", "SOX5", "TBR1", "LCORL", "ST18", "KAT6B"
  # Excitatory - upper or deep layer
  , "NEUROD6", "MAP2", "YWHAB"
  # Interneuron
  , "DLX1", "DLX5", "CXCR4", "CALB2", "CITED2"), rnames$gene_symbol)
df1 <- rnames[idx, ]
ssBsDF <- ssBsDF[df1$row_num, ]
row.names(ssBsDF) <- df1$gene_symbol
# Genes as columns
ssBsDF <- as.data.frame(t(ssBsDF))
# Remove outliers by stdev
ssBsDF <- as.data.frame(apply(ssBsDF, 2, function(x) {
  Remove_Outliers_By_SD(x, nStdev = 2.5)
  }))
# Add age
ssBsDF$AGE <- factor(df$age, levels = unique(bsMtDF$age))

# mean expression at each time point
ssBsDF <- melt(ssBsDF)
mnDF <- aggregate(value~variable+AGE, data = ssBsDF, mean, na.rm = TRUE)

# Split by gene and plot
ldf <- split(ssBsDF, ssBsDF$variable)
# Duplicate some genes to plot next to genes of interest
idx <- match(c(
  # RG
  "VIM", "HES1", "HOPX", "ITGB5", "CARHSP1", "ZFHX4", "LITAF", "MAFF"
  # Excitatory - not deep layer
  , "TUBB3", "STMN2", "SATB2", "CUX1", "CSRP2"
  # Excitatory deep layer
  , "BCL11B", "SOX5", "TBR1", "LCORL", "ST18", "KAT6B"
  # Excitatory - upper or deep layer
  , "TUBB3", "STMN2", "NEUROD6", "MAP2", "YWHAB"
  # Interneuron
  , "DLX1", "DLX5", "CXCR4", "CALB2", "CITED2"), names(ldf))
ldf <- ldf[idx]
names(ldf)
# Loop through and plot
ggL <- lapply(ldf, function(df) {
  # ggplot
  ggplot(df, aes(x = AGE, y = value, group = 1)) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    geom_jitter(size = 0.1, width = 0.2) +
    stat_summary(fun.y = "mean", color = "red", geom = "line") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # ggtitle(df$variable[1])
    xlab("Age") +
    ylab("Normalized expression")
})
# plot grid
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n"
  , "\nExpression across Brainspan developmental transcriptome"
  , "\nOutliers >2.5 SD from mean expression removed"
  , "\n"
)
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.075, 1))
# Save
ggsave(paste0(outGraph, "BrainSpan.png"), width = 14, height = 30)
################################################################################

### Overlap with Miller zones

miExDF <- MillerExpr
miExDF$ENTREZ <- row.names(miExDF)
miExDF <- melt(miExDF, id.vars = "ENTREZ")
miExDF$LAYER <- MillerMetaRAW$structure_acronym[
  match(miExDF$variable, MillerMetaRAW$well_id)]
# Convert layers to acronyms
miExDF$LAYER <- gsub(".*SG.*", "SG", millerZones)
miExDF$LAYER <- gsub(".*MZ.*", "MZ", millerZones)
miExDF$LAYER <- gsub(".*CP.*", "CP", millerZones)
miExDF$LAYER <- gsub(".*SP.*", "SP", millerZones)
miExDF$LAYER <- gsub(".*IZ.*", "IZ", millerZones)
miExDF$LAYER <- gsub(".*SZ.*", "SZ", millerZones)
miExDF$LAYER <- gsub(".*VZ.*", "VZ", millerZones)
# Add hgnc_symbol
miExDF$hgnc_symbol <- bmDF$hgnc_symbol[match(miExDF$ENTREZ, bmDF$entrezgene)]

## TFs
df <- deDF[deDF$ensembl_gene_id %in% tfDF$V1, ]
df <- miExDF[miExDF$hgnc_symbol %in% df$gene, ]

ldf <- split(df, df$ENTREZ)
llmod <- lapply(ldf, function(df) {
  lm(value ~ LAYER, data = df)
})

v1 <- sapply(llmod, function(lmod) {anova(lmod)$"Pr(>F)"[1] < 0.05})
v2 <- sapply(ldf, function(df) {df$hgnc_symbol[1]})
# OLIG1   FOSB   PBX1 GATAD1  PRDM2
genes <- v2[v1]

ggDF <- miExDF[miExDF$hgnc_symbol %in% genes, ]
ggplot(ggDF, aes(x = LAYER, y = value)) +
  facet_wrap(~hgnc_symbol, ncol = 4) +
  geom_boxplot()
ggsave(paste0(outGraph, "Miller_ANOVA.pdf"), width = 9, height = 9)

ggDF <- df
ggplot(ggDF, aes(x = LAYER, y = value)) +
  facet_wrap(~hgnc_symbol, ncol = 4) +
  geom_boxplot()
ggsave(paste0(outGraph, "Miller.pdf"), width = 9, height = 49)


## TFs, Chromatin remodelers, co-factors
df <- rbind(tfDF, crDF, cfDF)
df <- deDF[deDF$ensembl_gene_id %in% df$V1, ]
df <- miExDF[miExDF$hgnc_symbol %in% df$gene, ]

ldf <- split(df, df$ENTREZ)
llmod <- lapply(ldf, function(df) {
  lm(value ~ LAYER, data = df)
})

v1 <- sapply(llmod, function(lmod) {anova(lmod)$"Pr(>F)"[1] < 0.05})
v2 <- sapply(ldf, function(df) {df$hgnc_symbol[1]})
# OLIG1   FOSB   PBX1 GATAD1  PRDM2
genes <- v2[v1]

ggDF <- miExDF[miExDF$hgnc_symbol %in% genes, ]
ggplot(ggDF, aes(x = LAYER, y = value)) +
  facet_wrap(~hgnc_symbol, ncol = 4) +
  geom_boxplot()
ggsave(paste0(outGraph, "Miller_TFsCRsCFs_ANOVA.pdf"), width = 9, height = 9)

ggDF <- df
ggplot(ggDF, aes(x = LAYER, y = value)) +
  facet_wrap(~hgnc_symbol, ncol = 4) +
  geom_boxplot()
ggsave(paste0(outGraph, "Miller_TFsCRsCFs.pdf"), width = 9, height = 49)
################################################################################

### Overlap with He layer expression markers and human specific

## TFs, Chromatin remodelers, co-factors DE in clusters
df1 <- rbind(tfDF, crDF, cfDF)
df1 <- deDF[deDF$ensembl_gene_id %in% df1$V1, ]

df2 <- heDF
df2$EnsemblID <- gsub("\\..*", "", df2$EnsemblID)
df3 <- merge(df1, df2, by.x = "ensembl_gene_id", by.y = "EnsemblID")

# Stacked bar chart
ggDF <- rbind(
  data.frame(round(table(heDF$Layer.marker.in.human) /
      sum(table(heDF$Layer.marker.in.human))*100, 1), DATASET = "He et al.")
  , data.frame(round(table(df3$Layer.marker.in.human) /
      sum(table(df3$Layer.marker.in.human))*100, 1), DATASET = "Intersect He + Drop-seq DE")
)
colnames(ggDF)[colnames(ggDF) == "Var1"] <- "LAYER"
ggplot(ggDF, aes(x = DATASET, y = Freq)) +
  geom_bar(aes(fill = LAYER), stat = "identity") +
  scale_fill_brewer(type = "qual", palette = 6) +
  ylab("Percentage") +
  xlab("Dataset") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercentage of layer specific genes relative to all layer specific genes"
    , "\nTFs differentially expressed in a cluster"
    , "\nHuman layer specific - He et al. 2017"
    , "\nSeurat normalized expression"
    , "\n"))
ggsave(paste0(outGraph, "TFsCRsCFs_DE_HumanLayer_StackedBar.pdf")
  , width = 6, height = 6)

# Layer markers
# Heatmap
# Normalized, no mean centering scaling
df <- df3[! is.na(df3$Layer.marker.in.human), ]
ldf <- split(df, df$Layer.marker.in.human)
lapply(ldf, function(df) {
  layer <- df$Layer.marker.in.human[1]
  ggDF <- merge(df[c("gene", "cluster")], seuratO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$cluster <- as.factor(ggDF$cluster)
  Heatmap_By_Cluster(ggDF, seuratO = seuratO, mapDF = mapDF, lowerLimit = 0
    , upperLimit = 2
    , ggTitle = paste0(graphCodeTitle
      , "\n"
      , "\nTFs differentially expressed in a cluster"
      , "\nHuman layer specific - He et al. 2017"
      , "\nLayer: ", layer
      , "\nSeurat normalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "TFsCRsCFs_DE_HumanLayer", layer, "_Heatmap_Normalized.png")
    , width = 12, height = (5 + 0.15*nrow(df)))
})

# Species specific
# Heatmap
# Normalized, no mean centering scaling
ldf <- list(HUMAN_SPECIFIC_CHANGED = df3[df3$Human.specific.changed == TRUE, ]
  , CHIMPANZEE_SPECIFIC_CHANGED = df3[df3$Chimpanzee.specific.changed == TRUE, ]
  , MACAQUE_SPECIFIC_CHANGED = df3[df3$Macaque.specific.changed == TRUE, ])
lapply(names(ldf), function(name) {
  df <- ldf[[name]]
  ggDF <- merge(df[c("gene", "cluster")], seuratO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$cluster <- as.factor(ggDF$cluster)
  Heatmap_By_Cluster(ggDF, seuratO = seuratO, mapDF = mapDF, lowerLimit = 0
    , upperLimit = 2
    , ggTitle = paste0(graphCodeTitle
      , "\n"
      , "\nTFs differentially expressed in a cluster"
      , "\nSpecies specific changed - He et al. 2017"
      , "\n", name
      , "\nSeurat normalized expression"
      , "\n")
  )
  ggsave(paste0(outGraph, "TFsCRsCFs_DE_Species_", name, "_Heatmap_Normalized.png")
    , width = 12, height = (5 + 0.15*nrow(df)))
})

## DE in clusters

df2 <- heDF
df2$EnsemblID <- gsub("\\..*", "", df2$EnsemblID)
df3 <- merge(deDF, df2, by.x = "ensembl_gene_id", by.y = "EnsemblID")

# Stacked bar chart
ggDF <- rbind(
  data.frame(round(table(heDF$Layer.marker.in.human) /
      sum(table(heDF$Layer.marker.in.human))*100, 1), DATASET = "He et al.")
  , data.frame(round(table(df3$Layer.marker.in.human) /
      sum(table(df3$Layer.marker.in.human))*100, 1), DATASET = "Intersect He + Drop-seq DE")
)
colnames(ggDF)[colnames(ggDF) == "Var1"] <- "LAYER"
ggplot(ggDF, aes(x = DATASET, y = Freq)) +
  geom_bar(aes(fill = LAYER), stat = "identity") +
  scale_fill_brewer(type = "qual", palette = 6) +
  ylab("Percentage") +
  xlab("Dataset") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercentage of layer specific genes relative to all layer specific genes"
    , "\nGenes differentially expressed in a cluster"
    , "\nHuman layer specific - He et al. 2017"
    , "\nSeurat normalized expression"
    , "\n"))
ggsave(paste0(outGraph, "Genes_DE_HumanLayer_StackedBar.pdf")
  , width = 6, height = 6)

# Layer markers
# Heatmap
# Normalized, no mean centering scaling
df <- df3[! is.na(df3$Layer.marker.in.human), ]
ldf <- split(df, df$Layer.marker.in.human)
lapply(ldf, function(df) {
  layer <- df$Layer.marker.in.human[1]
  ggDF <- merge(df[c("gene", "cluster")], seuratO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$cluster <- as.factor(ggDF$cluster)
  Heatmap_By_Cluster(ggDF, seuratO = seuratO, mapDF = mapDF, lowerLimit = 0
    , upperLimit = 2
    , ggTitle = paste0(graphCodeTitle
      , "\n"
      , "\nGenes differentially expressed in a cluster"
      , "\nHuman layer specific - He et al. 2017"
      , "\nLayer: ", layer
      , "\nSeurat normalized expression"
      , "\n")
  )
  plotHeight <- (5 + 0.15*nrow(df))
  if (plotHeight > 49) {plotHeight <- 49}
  ggsave(paste0(outGraph, "Genes_DE_HumanLayer", layer, "_Heatmap_Normalized.png")
    , width = 12, height = plotHeight)
})

# Species specific
# Heatmap
# Normalized, no mean centering scaling
ldf <- list(HUMAN_SPECIFIC_CHANGED = df3[df3$Human.specific.changed == TRUE, ]
  , CHIMPANZEE_SPECIFIC_CHANGED = df3[df3$Chimpanzee.specific.changed == TRUE, ]
  , MACAQUE_SPECIFIC_CHANGED = df3[df3$Macaque.specific.changed == TRUE, ])
lapply(names(ldf), function(name) {
  df <- ldf[[name]]
  ggDF <- merge(df[c("gene", "cluster")], seuratO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$cluster <- as.factor(ggDF$cluster)
  Heatmap_By_Cluster(ggDF, seuratO = seuratO, mapDF = mapDF, lowerLimit = 0
    , upperLimit = 2
    , ggTitle = paste0(graphCodeTitle
      , "\n"
      , "\nTFs differentially expressed in a cluster"
      , "\nSpecies specific changed - He et al. 2017"
      , "\n", name
      , "\nSeurat normalized expression"
      , "\n")
  )
  plotHeight <- (5 + 0.15*nrow(df))
  if (plotHeight > 49) {plotHeight <- 49}
  ggsave(paste0(outGraph, "Genes_DE_Species_", name, "_Heatmap_Normalized.png")
    , width = 12, height = plotHeight)
})

## All He genes

df2 <- heDF
df2$EnsemblID <- gsub("\\..*", "", df2$EnsemblID)
df3 <- merge(deDF, df2, by.x = "ensembl_gene_id", by.y = "EnsemblID")

# Stacked bar chart
ggDF <- rbind(
  data.frame(round(table(heDF$Layer.marker.in.human) /
      sum(table(heDF$Layer.marker.in.human))*100, 1), DATASET = "He et al.")
  , data.frame(round(table(df3$Layer.marker.in.human) /
      sum(table(df3$Layer.marker.in.human))*100, 1), DATASET = "Intersect He + Drop-seq DE")
)
colnames(ggDF)[colnames(ggDF) == "Var1"] <- "LAYER"
ggplot(ggDF, aes(x = DATASET, y = Freq)) +
  geom_bar(aes(fill = LAYER), stat = "identity") +
  scale_fill_brewer(type = "qual", palette = 6) +
  ylab("Percentage") +
  xlab("Dataset") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercentage of layer specific genes relative to all layer specific genes"
    , "\nGenes differentially expressed in a cluster"
    , "\nHuman layer specific - He et al. 2017"
    , "\nSeurat normalized expression"
    , "\n"))
ggsave(paste0(outGraph, "AllGenes_HumanLayer_StackedBar.pdf")
  , width = 6, height = 6)

# Layer markers
# Heatmap
# Normalized, no mean centering scaling
df <- heDF[! is.na(heDF$Layer.marker.in.human), ]
df$cluster <- NA
ldf <- split(df, df$Layer.marker.in.human)
lapply(ldf, function(df) {
  layer <- df$Layer.marker.in.human[1]
  ggDF <- merge(df[c("Gene.symbol", "cluster")], seuratO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF$cluster <- as.factor(ggDF$cluster)
  Heatmap_By_Cluster(ggDF, seuratO = seuratO, mapDF = mapDF, lowerLimit = 0
    , upperLimit = 2
    , ggTitle = paste0(graphCodeTitle
      , "\n"
      , "\nGenes differentially expressed in a cluster"
      , "\nHuman layer specific - He et al. 2017"
      , "\nLayer: ", layer
      , "\nSeurat normalized expression"
      , "\n")
  )
  plotHeight <- (5 + 0.15*nrow(df))
  if (plotHeight > 49) {plotHeight <- 49}
  ggsave(paste0(outGraph, "AllGenes_HumanLayer", layer, "_Heatmap_Normalized.png")
    , width = 12, height = plotHeight)
})
################################################################################

### Human specific expression of TFs

# TFs, chromatin remodelers, co-factors
tccDF <- rbind(tfDF, crDF, cfDF)
tccDF$TYPE <- c(rep("TF", nrow(tfDF))
  , rep("Chromatin remodeler", nrow(crDF))
  , rep("Co-factor", nrow(cfDF)))

ssDeDF <- deDF[deDF$ENSEMBL %in% tccDF$V1, ]
ssDeDF[ssDeDF$ENSEMBL %in% aheCpDF$Identifier, ]
ssDeDF[ssDeDF$ENSEMBL %in% aheGzDF$Identifier, ]

ssDeDF[ssDeDF$GENE %in% hsDF$Gene[hsDF$Set == "Human-specific"], ]
ssDeDF[ssDeDF$GENE %in% hsDF$Gene[hsDF$Set == "Primate-specific"], ]
 ################################################################################

### Expression ranking in cell classes of interest

# RG - CARHSP1, ZFHX4
# MALAT1 as check
v1 <- rowMeans(noCentExM[ ,centSO@ident %in% c(7, 9)])
df <- data.frame(v1, rank(-v1))
df[row.names(df) %in% c("CARHSP1", "ZFHX4", "MALAT1"), ]
# ZFHX4   0.7648820       289
# CARHSP1 0.9960629       209
# MALAT1  6.0503535         1

# Excitatory - not deep layer
# MALAT1 as check
v1 <- rowMeans(noCentExM[ ,centSO@ident %in% c(0, 1, 4, 13)])
df <- data.frame(v1, rank(-v1))
df[row.names(df) %in% c("CSRP2", "MALAT1"), ]
# CSRP2  0.9643488       217
# MALAT1 6.2747075         1

# Excitatory - deep layer
# MALAT1 as check
v1 <- rowMeans(noCentExM[ ,centSO@ident %in% c(3, 14)])
df <- data.frame(v1, rank(-v1))
df[row.names(df) %in% c("ST18", "KAT6B", "MALAT1"), ]
# ST18   0.6913903       386
# KAT6B  1.0823028       195
# MALAT1 6.0458288         1

# Interneuron
# MALAT1 as check
v1 <- rowMeans(noCentExM[ ,centSO@ident %in% c(5, 6)])
df <- data.frame(v1, rank(-v1))
df[row.names(df) %in% c("CITED2", "MALAT1"), ]
# CITED2 0.8660202       248
# MALAT1 6.3206253         1
################################################################################

### Feature plot of TFs of interest

genes <- c("CARHSP1", "ZFHX4", "CSRP2", "ST18", "KAT6B", "CITED2")
  ggL <- FeaturePlot_CentScale(
    genes = genes
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO, limLow = -2, limHigh = 2
    )
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\ntSNE colored by normalized centered scaled expression"))
ggsave(paste0(outGraph, "TFsOfInterest_FeaturePlot.png"), width = 13, height = 10)
################################################################################

### Human specific expression of TFs with cell type specific expression

## Check genes of interest

# Human specific enhancers from Luis
genesDF <- bmDF[
  bmDF$hgnc_symbol %in% c("CARHSP1", "ZFHX4", "CSRP2", "ST18", "KAT6B", "CITED2"), ]
genesDF[genesDF$ensembl_gene_id %in% aheCpDF$Identifier, ]
# No intersection with CP list, KAT6B with GZ
genesDF[genesDF$ensembl_gene_id %in% aheGzDF$Identifier, ]
# X ensembl_gene_id hgnc_symbol entrezgene   gene_biotype
# 10240 10240 ENSG00000156650       KAT6B      23522 protein_coding

# Species specific genes from Bakken et al. Allen Macaque
hsDF[hsDF$Gene %in% genesDF$hgnc_symbol, ]
# Gene       Set Rhesus.Human Rhesus.Rat Human.Rat Human.Human.2
# 258  CARHSP1 Conserved    0.9877877  0.9949344 0.9978014            NA
# 435    CSRP2 Conserved    0.9864306  0.9825439 0.9973885     0.9996054
# 1775    ST18 Conserved    0.9631285  0.9890364 0.9228890     0.9648314
# 2039   ZFHX4 Conserved    0.7089200  0.7191158 0.9847448     0.9020230
################################################################################

### Human specific TFs expression

# Data frame of TFs, co-factors, chromatin remodelers
tccDF <- rbind(tfDF, crDF, cfDF)
# 1970
nrow(tccDF)
# Add gene symbol
tccDF$GENE <- bmDF$hgnc_symbol[match(tccDF$V1, bmDF$ensembl_gene_id)]

# Luis CP human specific enhancers intersect with TFs, co-factors, chromatin remodelers
# FALSE  TRUE
# 1928    42
table(tccDF$V1 %in% aheCpDF$Identifier)

# Luis GZ human specific enhancers intersect with TFs, co-factors, chromatin remodelers
# FALSE  TRUE
# 1936    34
table(tccDF$V1 %in% aheGzDF$Identifier)

# Species specific genes from Bakken et al. Allen Macaque intersect with TFs, co-factors, chromatin remodelers
# FALSE  TRUE
# 1890   172
table(hsDF$Gene %in% tccDF$GENE)
hsDF[hsDF$Gene %in% tccDF$GENE, ]
################################################################################

### Check if genes are TFs, chromatin remodelers, or co-factors

genes <- c("CARHSP1", "ZFHX4", "CSRP2", "ST18", "KAT6B", "CITED2")
genesDF <- bmDF[
  bmDF$hgnc_symbol %in% genes, ]

# TFs
genesDF[genesDF$ensembl_gene_id %in% tfDF$V1, ]
# X ensembl_gene_id hgnc_symbol entrezgene   gene_biotype
# 1992 1992 ENSG00000091656       ZFHX4      79776 protein_coding
# 9198 9198 ENSG00000147488        ST18       9705 protein_coding
# 9836 9836 ENSG00000153048     CARHSP1      23589 protein_coding

# Chromatin remodelers
genesDF[genesDF$ensembl_gene_id %in% crDF$V1, ]
#  X ensembl_gene_id hgnc_symbol entrezgene   gene_biotype
# 10240 10240 ENSG00000156650       KAT6B      23522 protein_coding

# Human co-factors
genesDF[genesDF$ensembl_gene_id %in% cfDF$V1, ]
# X ensembl_gene_id hgnc_symbol entrezgene   gene_biotype
# 11613 11613 ENSG00000164442      CITED2      10370 protein_coding
# 14144 14144 ENSG00000175183       CSRP2       1466 protein_coding
################################################################################

### Check if genes have binding site enrichment in DE genes

# None of the genes are in TRANSFAC or JASPAR...

jasparDF$alias[jasparDF$alias %in% genes]
tbejDF$ALIAS[tbejDF$ALIAS %in% genes]

transfacDF$alias[transfacDF$alias %in% genes]
tbetDF$ALIAS[tbetDF$ALIAS %in% genes]
################################################################################

### Overlap genes of interest with appropiate marker genes for ISH

## Number of cells each marker expressed in

genes <- c("CARHSP1", "ZFHX4", "CSRP2", "ST18", "KAT6B", "CITED2"
  , "BCL11B", "SOX5", "SATB2", "LHX6", "SOX2", "PAX6", "HOPX", "CRYAB"
  , "EOMES", "STMN2")
m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
l1 <- apply(m1, 2, function(col) row.names(m1)[col])
# number of cells each marker expressed in
df5 <- data.frame(Number = rowSums(m1))
df5$Gene <- factor(row.names(df5), levels = row.names(df5))

# Plot number of cells each marker expressed in
ggplot(df5, aes(x = Gene, y = Number)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cells expressing genes"
    , "\n(> 0.5 normalized expression)"))
ggsave(paste0(outGraph, "TFsOfInterest_And_Markers_Number_Barplot.png")
  , width = 5, height = 5)


## tSNE colored by intersection and heatmap of numbers of intersections

# RG
genes <- c("CARHSP1", "ZFHX4", "SOX2", "PAX6", "HOPX", "CRYAB")

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"
    , "\nRG"))
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Marker_Intersection_tSNE_RG.png")
  , width = 13, height = 13)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\nRG")
  )
ggsave(paste0(outGraph, "TFsOfInterest_And_Markers_NumberIntersect_Heatmap_RG.png")
  , width = 7, height = 7)


# Excitatory
genes <- c("CSRP2", "SATB2", "LHX2")

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"
    , "\nExcitatory"))
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Marker_Intersection_tSNE_Excitatory.png")
  , width = 13, height = 10)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\nExcitatory")
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Markers_NumberIntersect_Heatmap_Excitatory.png")
  , width = 7, height = 7)


# Deep layer excitatory
genes <- c("ST18", "KAT6B", "BCL11B", "SOX5")

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"
    , "\nExcitatory deep layer"))
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Marker_Intersection_tSNE_ExcDeepLayer.png")
  , width = 13, height = 10)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\nExcitatory deep layer")
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Markers_NumberIntersect_Heatmap_ExcDeepLayer.png")
  , width = 7, height = 7)


# Interneuron
genes <- c("CITED2", "DLX1", "DLX2", "DLX5", "DLX6")

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"
    , "\nInterneuron"))
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Marker_Intersection_tSNE_Interneuron.png")
  , width = 13, height = 17)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\nInterneuron")
)
ggsave(paste0(outGraph, "TFsOfInterest_And_Markers_NumberIntersect_Heatmap_Interneuron.png")
  , width = 7, height = 7)
################################################################################