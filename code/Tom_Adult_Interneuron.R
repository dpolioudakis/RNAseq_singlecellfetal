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
require(viridis)
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# Tom cluster mean expression tables
adult_mn_expr <- read.table(
  "../analysis/analyzed_data/Tom_Adult_Interneuron/Adult.Subtypes.tsv"
  , header = TRUE, row.names = 1)
calb2_mn_expr <- read.table(
  "../analysis/analyzed_data/Tom_Adult_Interneuron/Fetal.CALB2_non_unique.tsv", header = TRUE, row.names = 1)
sst_mn_expr <- read.table(
  "../analysis/analyzed_data/Tom_Adult_Interneuron/Fetal.SST_non_unique.tsv", header = TRUE, row.names = 1)

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

## Variables
graphCodeTitle <- "Tom_Adult_Interneuron.R"
outGraph <- "../analysis/graphs/Tom_Adult_Interneuron/Tom_Adult_Interneuron_"
outTable <- "../analysis/tables/Tom_Adult_Interneuron/Tom_Adult_Interneuron_"

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

################################################################################

### Format

## Adult interneuron enrichment
adult_enrich_M <- matrix(NA, nrow(adult_mn_expr), ncol(adult_mn_expr))
rownames(adult_enrich_M) <- rownames(adult_mn_expr)
colnames(adult_enrich_M) <- colnames(adult_mn_expr)
for (i in 1:ncol(adult_mn_expr)){
  enrich <- adult_mn_expr[ ,i] / rowMeans(adult_mn_expr[,-i])
  adult_enrich_M[,i] <- enrich
}

## TADA Sanders_TADA0.1_exomedel
tadaDF <- tadaDF[,1:21]
tada <- unique(tadaDF[tadaDF$tadaFdrAscSscExomeSscAgpSmallDel<0.1, "RefSeqGeneName"])

## iHART
ihartDF$HGNC.gene.symbol <- gsub("\"", "", ihartDF$HGNC.gene.symbol)

## Subset to high confidence epilepsy genes
epilepsy_high_conf <- epilepsy_DF$Gene[
  epilepsy_DF$Classification == "High-confidence"]
################################################################################

### Data processing

## Adult interneuron enrichment
adult_enrich_M <- matrix(NA, nrow(adult_mn_expr), ncol(adult_mn_expr))
rownames(adult_enrich_M) <- rownames(adult_mn_expr)
colnames(adult_enrich_M) <- colnames(adult_mn_expr)
for (i in 1:ncol(adult_mn_expr)){
  enrich <- adult_mn_expr[ ,i] - rowMeans(adult_mn_expr[,-i])
  adult_enrich_M[,i] <- enrich
}

## Unique disease genes
# ASD
tada_unique <- tada[
  ! tada %in% epilepsy_high_conf
  & ! tada %in% id_genes]
# Epilepsy
epilepsy_unique <- epilepsy_high_conf[
  ! epilepsy_high_conf %in% tada
  & ! epilepsy_high_conf %in% id_genes]
# ID
id_unique <- id_genes[
  ! id_genes %in% tada
  & ! id_genes %in% epilepsy_high_conf]
################################################################################

### ASD TADA heatmaps

# Adult
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% tada, ]
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of ASD risk genes"
    , "\n(Sanders et al)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "ASDtada_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 9)

# Adult unique
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% tada_unique, ]
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of unique ASD risk genes"
    , "\n(Sanders et al)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "ASDtadaUnique_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 9)

# Fetal
fetal_mn_expr <- cbind(calb2_mn_expr, sst_mn_expr)
gg_DF <- fetal_mn_expr[rownames(fetal_mn_expr) %in% tada, ]
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
  # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
  #   , palette = 5, direction = -1, limits = c(lowerLimit, upperLimit)
  #   , na.value = "grey90") +
  scale_fill_distiller(name = "Normalized\nexpression", type = "div"
    , palette = "RdYlBu", direction = -1, na.value = "grey90", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  # theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
ggsave(paste0(outGraph, "ASDtada_Fetal_Interneuron_Heatmap.png"), width = 5, height = 8)
################################################################################

### Epilepsy heatmaps

# Adult
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% epilepsy_high_conf
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of Epilepsy risk genes"
    , "\n(Elizabeth Ruzzo compiled)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "Epilepsy_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 14)

# Adult
# Unique epilepsy genes
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% epilepsy_unique
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of unique Epilepsy risk genes"
    , "\n(Elizabeth Ruzzo compiled)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "EpilepsyUnique_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 14)
################################################################################

### ID heatmaps

# Adult
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% id_genes, ]
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of ID risk genes"
    , "\n(Luis de la Torre-Ubieta compiled)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "ID_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 7)

# Adult
# Unique ID genes
gg_DF <- adult_mn_expr[rownames(adult_mn_expr) %in% id_unique, ]
gg_DF$Gene <- rownames(gg_DF)
gg_DF <- melt(gg_DF)
gg_DF$value[gg_DF$value > 2] <- 2
gg_DF$value[gg_DF$value < 0] <- 0
ggplot(gg_DF, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  # scale_fill_distiller(name = "Normalized\nexpression", type = "seq"
  #   , palette = "YlGnBu", direction = -1, na.value = "grey90"
  #   , limits = c(0,2)) +
  scale_fill_viridis(name = "Expression", limits = c(0,2)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(strip.background = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Clusters") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nExpression of unique ID risk genes"
    , "\n(Luis de la Torre-Ubieta compiled)"
    , "\nAdult interneuron clusters"))
ggsave(paste0(outGraph, "IDunique_Adult_Interneuron_Heatmap.png")
  , width = 4, height = 7)
################################################################################

### Enrichment (log2 odds ratio)

Add_Gene_List_To_Binary_Matrix <- function(
  gene_binary_M, gene_list, gene_list_name){
  print("Add_Gene_List_To_Binary_Matrix")
  binary_gene_list <- as.numeric(rownames(gene_binary_M) %in% gene_list)
  gene_binary_M <- cbind(gene_binary_M, binary_gene_list)
  colnames(gene_binary_M)[dim(gene_binary_M)[2]] <- gene_list_name
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
gene_binary_M <- matrix(NA, nrow(adult_enrich_M), 0)
rownames(gene_binary_M) <- rownames(adult_enrich_M)

# Add cluster enriched
cluster_genes_LL <- lapply(colnames(adult_enrich_M), function(cluster){
  cluster_enrich <- adult_enrich_M[ ,colnames(adult_enrich_M) %in% cluster]
  # names(cluster_enrich)[cluster_enrich > 0.2]
  cluster_enrich <- sort(cluster_enrich, decreasing = TRUE)
  names(cluster_enrich)[1:250]
})
names(cluster_genes_LL) <- paste0("Enriched ", colnames(adult_enrich_M))
gene_binary_M <- Add_Gene_Lists_To_Binary_Matrix(
  gene_binary_M, gene_lists_LL = cluster_genes_LL
)
# Add top expressed for each cluster
cluster_genes_LL <- lapply(colnames(adult_mn_expr), function(cluster){
  cluster_enrich <- as.matrix(adult_mn_expr)[
    ,colnames(adult_mn_expr) %in% cluster]
  print(head(cluster_enrich))
  # names(cluster_enrich)[cluster_enrich > 0.2]
  cluster_enrich <- sort(cluster_enrich, decreasing = TRUE)
  names(cluster_enrich)[1:250]
})
names(cluster_genes_LL) <- paste0("Top expressed ", colnames(adult_mn_expr))
gene_binary_M <- Add_Gene_Lists_To_Binary_Matrix(
  gene_binary_M, gene_lists_LL = cluster_genes_LL
)

# Add disease gene lists
# ASD
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M, gene_list = tada, gene_list_name = "ASD")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M
  , gene_list = epilepsy_high_conf, gene_list_name = "Epilepsy")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M, gene_list = id_genes, gene_list_name = "ID")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M, gene_list = tada_unique
  , gene_list_name = "ASD_unique")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M
  , gene_list = epilepsy_unique, gene_list_name = "Epilepsy_unique")
gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
  gene_binary_M = gene_binary_M, gene_list = id_unique
  , gene_list_name = "ID_unique")

## Calculate odds ratio
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
    log2OR_DF$Number[n] <- sum(dat1 == 1 & dat2 == 1)
  }
}

# Plot
gg_DF <- log2OR_DF[log2OR_DF$Var1 %in% c("ASD", "Epilepsy", "ID"), ]
gg_DF <- gg_DF[! gg_DF$Var2 %in% c(
  "ASD_unique", "Epilepsy_unique", "ID_unique", "ASD", "Epilepsy", "ID"), ]
gg_DF <- gg_DF[gg_DF$Var1 != gg_DF$Var2, ]
gg_DF$FDR <- p.adjust(gg_DF$Pvalue, method="BH")
gg_DF$Method <- ifelse(
  grepl("Enriched", gg_DF$Var2) == TRUE, "Enriched", "Top expressed")
gg_DF$Cluster <- gsub("Enriched ", "", gg_DF$Var2)
gg_DF$Cluster <- gsub("Top expressed ", "", gg_DF$Cluster)
ggplot(gg_DF, aes(x = Cluster, y = -log(FDR, 10), fill = Cluster)) +
  facet_wrap(~Var1+Method, scales = "free", ncol = 2) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(FDR, 10)
    , label = paste0(round(Log2_Odds_Ratio, 1), "\n(", Number, ")"))
    , position = position_dodge(width = 1), color = "black"
    # , angle = 90, hjust = 0.3
    , size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Clusters") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nEnrichment of disease risk genes"
    , "\nNumber = log2 odds ratio; () = size of intersection"
    , "\nAdult interneuron clusters"))
# ggsave("../analysis/graphs/MetaMat/MetaMat_ID_BarPlot_ReorderClusters.png"
#   , width = 10, height = 5)
ggsave(paste0(outGraph, "Enrichment_Disease_Barplot.pdf")
  , width = 6, height = 9)

# Plot unique disease genes
gg_DF <- log2OR_DF[log2OR_DF$Var1 %in% c(
  "ASD_unique", "Epilepsy_unique", "ID_unique"), ]
gg_DF <- gg_DF[! gg_DF$Var2 %in% c(
  "ASD_unique", "Epilepsy_unique", "ID_unique", "ASD", "Epilepsy", "ID"), ]
gg_DF <- gg_DF[gg_DF$Var1 != gg_DF$Var2, ]
gg_DF$FDR <- p.adjust(gg_DF$Pvalue, method="BH")
gg_DF$Method <- ifelse(
  grepl("Enriched", gg_DF$Var2) == TRUE, "Enriched", "Top expressed")
gg_DF$Cluster <- gsub("Enriched ", "", gg_DF$Var2)
gg_DF$Cluster <- gsub("Top expressed ", "", gg_DF$Cluster)
ggplot(gg_DF, aes(x = Cluster, y = -log(FDR, 10), fill = Cluster)) +
  facet_wrap(~Var1+Method, scales = "free", ncol = 2) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(FDR, 10)
    , label = paste0(round(Log2_Odds_Ratio, 1), "\n(", Number, ")"))
    , position = position_dodge(width = 1), color = "black"
    # , angle = 90, hjust = 0.3
    , size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Clusters") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nEnrichment of unique disease risk genes"
    , "\nNumber = log2 odds ratio; () = size of intersection"
    , "\nAdult interneuron clusters"))
# ggsave("../analysis/graphs/MetaMat/MetaMat_ID_BarPlot_ReorderClusters.png"
#   , width = 10, height = 5)
ggsave(paste0(outGraph, "Enrichment_DiseaseGeneUnique_Barplot.pdf")
  , width = 6, height = 9)
################################################################################
