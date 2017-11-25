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
################################################################################

### ASD combined de novo proband LOF iossifov + de Rubeis

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

## iossifov + de Rubeis
asdDF = merge(rubDF, iosDF, by = intersect("gene", "gene"))
asdDF$combined_dn_prob_LOF = asdDF$dnv_LGDs_prb + asdDF$Rubeis_dn.LoF
################################################################################

### ASD combined de novo proband LOF iossifov + de Rubeis

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- asdDF$gene[asdDF$combined_dn_prob_LOF > 0]

# Genes DE >0.4 in cluster 4 or cluster 14
df1 <- rbind(data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 4]
  , Cluster = 4)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 14]
    , Cluster = 14)
)

# Subset to ASD
df1 <- df1[df1$Gene %in% genes, ]

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
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.3, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of ASD combined de novo proband LOF iossifov + de Rubeis genes differentially expressed in cluster 4 or 14"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
)
ggsave(paste0(
  outGraph, "DE_IossifovRubeis_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 10, limitsize = FALSE)
################################################################################

### Human specific ASD combined de novo proband LOF iossifov + de Rubeis

# ASD combined de novo proband LOF iossifov + de Rubeis
genes <- asdDF$gene[asdDF$combined_dn_prob_LOF > 0]

# Allen human specific
genes <- genes[genes %in% hsDF$Gene[hsDF$Set == "Human-specific"]]

# Genes DE >0.4 in cluster 4 or cluster 14
df1 <- rbind(data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 4]
  , Cluster = 4)
  , data.frame(Gene = deDF$GENE[deDF$LOG_FC > 0.4 & deDF$CLUSTER == 14]
    , Cluster = 14)
)

# Subset to human specific + ASD
df1 <- df1[df1$Gene %in% genes, ]

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
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.5, align = 'h', axis = 'b'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of human specific ASD combined de novo proband LOF iossifov + de Rubeis genes differentially expressed in cluster 4 or 14"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\nDE filters:"
    , "\n> 0.4 log fold change cluster vs all other cells"
    , "\n")
)
ggsave(paste0(
  outGraph, "DE_HumanSpecific_IossifovRubeis_Heatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 6, limitsize = FALSE)
################################################################################