# Damon Polioudakis
# 2017-02-21
# Plot known marker lists as heatmaps, violin plots, feature plots

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
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# TF enrichment results
tf_enrichment_DF <- read.csv("TF_Enrichment/EnrichmentTFs_ATAConly/EnrichmentTFs_ATACandHiC_logistic_conserved_combined.csv", header = TRUE)

# Cluster DE table
cluster_DE_DF <- read.table(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

## Variables
graphCodeTitle <- "TF_Enrichment_Expression.R"
outGraph <- "../analysis/graphs/TF_Enrichment_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/TF_Enrichment_Expression_"
# outGraph <- "../analysis/graphs/TF_Enrichment_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/TF_Enrichment_Expression_"
# outGraph <- "../analysis/graphs/TF_Enrichment_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/TF_Enrichment_Expression_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################



# Initialize matrix
gene_binary_M <- matrix(NA, nrow(noCentExM), 0)
rownames(gene_binary_M) <- rownames(noCentExM)

Add_Gene_List_To_Binary_Matrix <- function(
  gene_binary_M, gene_list, gene_list_name){
  print("Add_Gene_List_To_Binary_Matrix")
  binary_gene_list <- as.numeric(rownames(gene_binary_M) %in% gene_list)
  gene_binary_M <- cbind(gene_binary_M, binary_gene_list)
  colnames(gene_binary_M)[dim(gene_binary_M)[2]] <- gene_list_name
  return(gene_binary_M)
}

Add_Dropseq_Cluster_Enriched_Genes_To_Binary_Matrix <- function(){
  cluster_genes_LL <- split(cluster_DE_DF$Gene, f = cluster_DE_DF$Cluster)
  for(i in 1:length(cluster_genes_LL)){
    gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
      gene_binary_M = gene_binary_M
      , gene_list = cluster_genes_LL[[i]]
      , gene_list_name = paste0("Dropseq_Cluster", names(cluster_genes_LL)[i])
    )
  }
  return(gene_binary_M)
}

Add_Gene_Lists_To_Binary_Matrix <- function(gene_binary_M, gene_lists_LL){
  for(i in 1:length(gene_lists_LL)){
    gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
      gene_binary_M = gene_binary_M
      , gene_list = gene_lists_LL[[i]]
      , gene_list_name = names(gene_lists_LL)[i]
    )
  }
  return(gene_binary_M)
}

# Initialize matrix
gene_binary_M <- matrix(NA, nrow(noCentExM), 0)
rownames(gene_binary_M) <- rownames(noCentExM)
# Add dropseq cluster enriched genes
cluster_genes_LL <- split(cluster_DE_DF$Gene, f = cluster_DE_DF$Cluster)
names(cluster_genes_LL) <- paste0("Dropseq_Cluster_", names(cluster_genes_LL))
gene_binary_M <- Add_Gene_Lists_To_Binary_Matrix(
  gene_binary_M, gene_lists_LL = cluster_genes_LL
)
# Add TF RE cluster enriched
subset_tf_enrichment_DF <- tf_enrichment_DF[tf_enrichment_DF$pval < 0.05, ]
cluster_genes_LL <- split(
  subset_tf_enrichment_DF$TFname
  , f = subset_tf_enrichment_DF$Cluster
)
names(cluster_genes_LL) <- paste0("TF_Cluster_", names(cluster_genes_LL))
gene_binary_M <- Add_Gene_Lists_To_Binary_Matrix(
  gene_binary_M, gene_lists_LL = cluster_genes_LL
)

log2OR_DF <- expand.grid(colnames(gene_binary_M), colnames(gene_binary_M))
log2OR_DF$Log2_Odds_Ratio <- NA
log2OR_DF$Pvalue <- NA
n <- 0
for(i in 1:ncol(gene_binary_M)){
  for(j in 1:ncol(gene_binary_M)){
    n <- n+1
    print(paste(colnames(gene_binary_M)[i],"vs",colnames(gene_binary_M)[j]))
    # Calculate enrichment with GLM
    dat1 <- as.numeric(gene_binary_M[,i])
    dat2 <- as.numeric(gene_binary_M[,j])
    glm.out <- glm(dat1~dat2, family = binomial)
    enrichP <- summary(glm.out)$coefficients[2,4]
    enrichOR <- summary(glm.out)$coefficients[2,1]
    # log2 odds ratio
    enrich_log2_OR <- log2(exp(enrichOR))
    log2OR_DF$Log2_Odds_Ratio[n] <- enrich_log2_OR
    log2OR_DF$Pvalue[n] <- enrichP
  }
}

# TF RE cluster x cluster enriched genes enrichment heatmap
gg_DF <- log2OR_DF[
  grepl("TF", log2OR_DF$Var1) &
  grepl("Dropseq", log2OR_DF$Var2)
  , ]
gg_DF$Var1 <- gsub("TF_Cluster_", "", gg_DF$Var1)
gg_DF$Var2 <- gsub("Dropseq_Cluster_", "", gg_DF$Var2)
gg_DF$Var1 <- factor(gg_DF$Var1
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
gg_DF$Var2 <- factor(gg_DF$Var2
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
ggplot(gg_DF, aes(x = Var1, y = Var2, fill = -log(Pvalue, 10))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "-log10 p-value") +
  geom_text(aes(x = Var1, y = Var2, label = signif(Log2_Odds_Ratio, 2))
    , color = "black", size = 3) +
  xlab("TF RE enrichment") +
  ylab("Drop-seq cluster enrichment") +
  ggplot_set_theme_publication
ggsave(paste0(outGraph, "_OddsRatio_Heatmap.png")
  , width = 9, height = 7)

# TF RE cluster x TF RE cluster enrichment heatmap
gg_DF <- log2OR_DF[
  grepl("TF", log2OR_DF$Var1) &
  grepl("TF", log2OR_DF$Var2)
  , ]
gg_DF$Var1 <- gsub("TF_Cluster_", "", gg_DF$Var1)
gg_DF$Var2 <- gsub("TF_Cluster_", "", gg_DF$Var2)
gg_DF$Var1 <- factor(gg_DF$Var1
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
gg_DF$Var2 <- factor(gg_DF$Var2
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
gg_DF$NegLog10Pvalue <- -log(gg_DF$Pvalue, 10)
ggplot(gg_DF, aes(x = Var1, y = Var2, fill = -log(Pvalue, 10))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "-log10 p-value") +
  geom_text(aes(x = Var1, y = Var2, label = signif(Log2_Odds_Ratio, 2))
    , color = "black", size = 3) +
  xlab("TF RE enrichment") +
  ylab("TF RE enrichment") +
  ggplot_set_theme_publication
ggsave(paste0(outGraph, "_TF_RE_OddsRatio_Heatmap.png")
  , width = 9, height = 7)


### Plots

## Feature plot mean expression

# Genes
gene_group_DF <- tf_enrichment_DF[tf_enrichment_DF$pval < 0.05, ]
gene_group_DF <- gene_group_DF[c("TFname", "Cluster")]
colnames(gene_group_DF) <- c("Gene", "Grouping")
# Collect tSNE values for ggplot
tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
# Normalized centered scaled
ggL <- FeaturePlot(
  genes = gene_group_DF$Gene
  , tsneDF = tsneDF
  , seuratO = centSO
  , exM = centSO@scale.data
  , limLow = -1.5
  , limHigh = 1.5
  , geneGrouping = gene_group_DF$Grouping
  , centScale = TRUE
)
ggL <- lapply(ggL, function(gg){
  gg <- gg + ggplot_set_theme_publication
  return(gg)
})
ggL[[1]] <- ggL[[1]] + theme(legend.position = "none")
Plot_Grid(
  ggPlotsL = ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of TFs with binding sites enriched in cluster enriched genes"
    , "\nNormalized centered scaled expression"
    , "\n")
)
ggsave(paste0(outGraph
    , "FeaturePlot_NormalizedCenteredScaled_paper.png"
  )
  , width = 12, height = 20, limitsize = FALSE
)

## Heatmap

# Genes
gene_group_DF <- tf_enrichment_DF[tf_enrichment_DF$pval < 0.05, ]
# Order by cluster and enrichment
gene_group_DF$Cluster <- factor(gene_group_DF$Cluster,
   levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
gene_group_DF <- gene_group_DF[
  order(gene_group_DF$Cluster, gene_group_DF$pval), ]
gene_group_DF <- gene_group_DF[c("TFname", "Cluster")]
colnames(gene_group_DF) <- c("Gene", "Group")
gene_group_DF <- gene_group_DF[! duplicated(gene_group_DF), ]
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = gene_group_DF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of TFs with binding sites enriched in cluster enriched genes"
    , "\nSorted by TF binding site enrichment"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "EnrichmentSorted_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 10, height = 60, limitsize = FALSE
)

# Subset to TFs enriched in any cluster
gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% cluster_DE_DF$Gene, ]
gg <- Plot_Marker_Genes_Heatmap_SetColWidths(
  geneGroupDF = gene_group_DF
  , exprM = centSO@scale.data
  , cellID_clusterID <- centSO@ident
)
gg + ggtitle(paste0(
  graphCodeTitle
    , "\n\nExpression of TFs with binding sites enriched in cluster enriched genes and enriched in a cluster"
    , "\nSorted by TF binding site enrichment"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
)
ggsave(paste0(outGraph
    , "EnrichmentSubset_HeatmapSetColWid_NormCenterScale_paper.png"
  )
  , width = 10, height = 60, limitsize = FALSE
)
################################################################################
