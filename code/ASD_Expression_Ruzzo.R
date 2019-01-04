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
require(tidyverse)
require(viridis)
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
lake_ex_df <- read.table(
  "../lake_2017/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", header = TRUE
)

# Luis metaMat results
mmapDF <- read.csv("../source/Gene_Lists/Overlapped-Genes.csv", header = TRUE)

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

# ASD TADA Sanders 2015 = from Luis metaMat
tadaDF <- read.csv(
  "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv")

# ASD Ruzzo
ihartDF <- read.csv("../source/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv")

# ASD Laura
asd_laura <- read.table(
  "../source/Gene_Lists/TADAms2_10_novel_ASD-risk_genes_FDRlt0.1.txt")[ ,1]

# Known marker Luis table
km_df <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv")

## Variables
script_name <- "ASD_Expression_Ruzzo.R"
out_graph <- "../analysis/graphs/ASD_Expression_Ruzzo/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/ASD_Expression_Ruzzo_"
out_table <- "../analysis/tables/ASD_Expression_Ruzzo/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/ASD_Expression_Ruzzo_"

## Output Directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

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


plot_mean_expression_by_cluster_heatmap <- function(
  expr_m
  , genes
  , cell_cluster_key
  , cluster_annot_key = NA
  , expr_lower_limit = -1
  , expr_upper_limit = 1
  , color_viridis = FALSE
  , keep_gene_order = TRUE){
  gg_tb <- expr_m %>%
    as_tibble(rownames = "gene") %>%
    # calculate mean expression by cluster of genes of interest
    filter(gene %in% genes) %>%
    gather(., key = "cell_id", value = "expression", -gene) %>%
    left_join(.
      , data.frame(cluster = cell_cluster_key
        , cell_id = names(cell_cluster_key))
    ) %>%
    left_join(.
      , data.frame(cluster = names(cluster_annot_key)
        , cluster_annot = cluster_annot_key)
    ) %>%
    group_by(gene, cluster_annot) %>%
    summarise(mean_expression = mean(expression)) %>% ungroup() %>%
    # order clusters
    mutate(cluster_annot = fct_relevel(cluster_annot, cluster_annot_key)) %>%
    # order genes
    mutate(gene = fct_relevel(gene, genes)) %>%
    # set expression limits
    mutate(mean_expression = replace(
      mean_expression, mean_expression > expr_upper_limit, expr_upper_limit)
    ) %>%
    mutate(mean_expression = replace(
      mean_expression, mean_expression < expr_lower_limit, expr_lower_limit))
  # plot
  ggplot(gg_tb, aes(x = cluster_annot, y = fct_rev(gene), fill = mean_expression)) +
      geom_tile() +
      {if(color_viridis == TRUE) {
        scale_fill_viridis(name = "Normalized expression"
        , limits = c(expr_lower_limit, expr_upper_limit))
      } else {
        scale_fill_distiller(name = "Normalized\nexpression\nzscore"
          , type = "div", palette = 5, direction = -1
          , limits = c(expr_lower_limit, expr_upper_limit)
          , na.value = "grey90")
      }} +
      theme_bw() +
      theme(
        axis.ticks = element_blank()
        , text = element_text(size = 12)
        , axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Gene") +
      xlab("Cluster")
}

order_genes_by_group_enrichment <- function(
  genes
  , de_df
  , enrichment_col
  , col_with_groups_to_order_by
  , col_with_gene
  , group_order) {
    # browser()
  print("order_genes_by_cluster_enrichment")
  genes <- as.character(genes)
  enrichment_col <- as.name(enrichment_col)
  col_with_groups_to_order_by <- as.name(col_with_groups_to_order_by)
  col_with_gene <- as.name(col_with_gene)
  ordered_genes <- de_df %>%
    as_tibble %>%
    mutate_if(is.factor, as.character) %>%
    # Filter group not in order list
    filter(UQ(col_with_groups_to_order_by) %in% group_order) %>%
    # Subset to genes
    filter(UQ(col_with_gene) %in% genes) %>%
    # Take group with highest enrichment for each gene
    group_by(UQ(col_with_gene)) %>%
    filter(UQ(enrichment_col) == max(UQ(enrichment_col))) %>%
    # Order by group
    mutate(group = factor(
      UQ(col_with_groups_to_order_by), levels = group_order)) %>%
    arrange(group) %>%
    pull(UQ(col_with_gene))
  # Add genes that were not enriched
  ordered_genes <- c(ordered_genes, genes[! genes %in% ordered_genes])
  return(ordered_genes)
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

### Plot pie charts

plot_pie_charts <- function(){
  print("plot_pie_charts")
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
    , file = paste0(out_table, "iHART69_CellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART 69 classified by cell type expression in human fetal brain"
      , "\nDrop-seq human fetal brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART69_CellTypes_Piechart.pdf"))


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
    , file = paste0(out_table, "iHART16novel_CellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART novel 16 classified by cell type expression in human fetal brain"
      , "\nDrop-seq human fetal brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART16novel_CellTypes_Piechart.pdf"))


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
  ggsave(paste0(out_graph, "iHART69_NowkowskiClusters_JitterPlot.png")
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
    , file = paste0(out_table, "iHART69_NowkowskiCellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART 69 classified by cell type expression in human fetal brain"
      , "\nNowkowski human fetal brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART69_NowkowskiCellTypes_Piechart.pdf"))

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
    , file = paste0(out_table, "iHART16novel_NowkowskiCellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART 16 novel classified by cell type expression in human fetal brain"
      , "\nNowkowski human fetal brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART16novel_NowkowskiCellTypes_Piechart.pdf"))


  ## iHART 69 pie chart with Lake adult dataset
  # iHART genes
  genes <- ihartDF$HGNC.gene.symbol[ihartDF$"iHART.69" == 1]
  # Check expression levels of iHART genes
  subset_ex_M <- lake_ex_df[
    rownames(lake_ex_df) %in% genes
    , ]
  gg_DF <- subset_ex_M
  gg_DF$Gene <- rownames(gg_DF)
  gg_DF <- melt(gg_DF)
  gg_DF <- gg_DF[! is.na(gg_DF$Gene), ]
  ggplot(gg_DF, aes(x = Gene, y = value)) +
    geom_jitter(size = 0.01) +
  ggsave(paste0(out_graph, "iHART69_Lake_JitterPlot.png")
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
    exM = lake_ex_df, de_DF = subset_DE_DF, genes = genes)
  # Number of genes per cell class
  number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
    subset_DE_DF, exM = lake_ex_df, ex_cutoff = 1, genes = genes)
  # Percent
  percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
  # Output table of percents
  write.csv(percent_classes_DF
    , file = paste0(out_table, "iHART69_LakeCellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART 69 classified by cell type expression in human adult brain"
      , "\nLake human adult brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART69_LakeCellTypes_Piechart.pdf"))



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
    exM = lake_ex_df, de_DF = subset_DE_DF, genes = genes)
  # Number of genes per cell class
  number_classes_DF <- Calculate_Number_Of_Genes_Per_Cell_Class(
    subset_DE_DF, exM = lake_ex_df, ex_cutoff = 1, genes = genes)
  # Percent
  percent_classes_DF <- Percent_Of_Genes_Per_Cell_Class(number_classes_DF)
  # Output table of percents
  write.csv(percent_classes_DF
    , file = paste0(out_table, "iHART16novel_LakeCellTypes_Piechart.csv")
  )
  # Plot
  ggplot(percent_classes_DF, aes(x = "", y = Percent, fill = Class)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggplot_set_theme_publication_nolabels +
    ggtitle(paste0(script_name
      , "\n\niHART 16 novel classified by cell type expression in human adult brain"
      , "\nLake human adult brain dataset"
      , "\nGlutamatergic enrichment:"
      , "\n\t\tLog2 odds ratio: ", signif(enrich_log2_OR_pval_L[[1]], 2)
      , "\n\t\tP-value: ", signif(enrich_log2_OR_pval_L[[2]], 2)))
  ggsave(paste0(out_graph, "iHART16novel_LakeCellTypes_Piechart.pdf"))
}
################################################################################

### Plot heatmaps

plot_expression_by_cluster_heatmaps <- function(){
  print("plot_expression_by_cluster_heatmaps")

  ## Heatmap of ASD Laura genes with fixed column widths

  # Drop-seq
  genes <- order_genes_by_group_enrichment(
    genes = asd_laura
    , de_df = deDF
    , enrichment_col = "Log2_Fold_Change"
    , col_with_groups_to_order_by = "Cluster"
    , col_with_gene = "Gene"
    , group_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  )
  geneGroupDF <- data.frame(
    Gene = rev(genes)
    , Group = ""
  )
  # Annotate cluster names
  cluster_annot_key <- c(
    "9" = "vRG"
    , "7" = "oRG"
    , "8" = "PgS"
    , "10" = "PgG2M"
    , "2" = "IP"
    , "0" = "ExN"
    , "1" = "ExM"
    , "4" = "ExCal"
    , "3" = "ExDp1"
    , "13" = "ExDp2"
    , "5" = "InSST"
    , "6" = "InCALB2"
    , "11" = "OPC"
    , "12" = "End"
    , "14" = "Per"
    , "15" = "Mic"
    , "16" = "NA"
  )
  cellID_clusterID_DF <- data.frame(centSO@ident)
  idx <- match(cellID_clusterID_DF[ ,1], names(cluster_annot_key))
  cellID_clusterID_DF$cluster <- cluster_annot_key[idx]
  cellID_clusterID <- cellID_clusterID_DF$cluster
  names(cellID_clusterID) <- row.names(cellID_clusterID_DF)
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID <- cellID_clusterID
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in fetal sorted by cluster enrichment"
        , "\nPolioudakis fetal dropseq dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression, mean centered, variance scaled"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ASDlaura_EnrichmentSorted_HeatmapSetColWid_NormCenterScale.png")
    , width = 10, height = 5
  )
  # Normalized
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID <- cellID_clusterID
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
    , lowerLimit = -0.5
    , upperLimit = 1.5
    , color_viridis = TRUE
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in fetal sorted by cluster enrichment"
        , "\nPolioudakis fetal dropseq dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ASDlaura_EnrichmentSorted_HeatmapSetColWid_Norm.png")
    , width = 10, height = 5
  )

  # Lake
  cellID_clusterID <- gsub("_.*", "", colnames(lake_ex_df))
  names(cellID_clusterID) <- colnames(lake_ex_df)
  clusterOrder <- c("Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
    , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
    , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = t(scale(t(lake_ex_df)))
    , cellID_clusterID = cellID_clusterID
    , clusterOrder = clusterOrder
    , clusters = clusterOrder
    , lowerLimit = -0.5
    , upperLimit = 1.5
    , color_viridis = TRUE
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
        , "\nLake et al. 2017 dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression, mean centered, variance scaled"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ASDlaura_Lake_EnrichmentSorted_HeatmapSetColWid_NormCenterScale.png")
    , width = 10, height = 5
  )
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = t(scale(t(lake_ex_df)))
    , cellID_clusterID = cellID_clusterID
    , clusterOrder = clusterOrder
    , clusters = clusterOrder
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
        , "\nLake et al. 2017 dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ASDlaura_Lake_EnrichmentSorted_HeatmapSetColWid_Norm.png")
    , width = 10, height = 5
  )

  ## Mean expression by cluster heatmaps
  # Dropseq
  # Normalized centered scaled
  genes <- order_genes_by_group_enrichment(
    genes = asd_laura
    , de_df = deDF
    , enrichment_col = "Log2_Fold_Change"
    , col_with_groups_to_order_by = "Cluster"
    , col_with_gene = "Gene"
    , group_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  )
  plot_mean_expression_by_cluster_heatmap(
    expr_m = centSO@scale.data
    , genes = genes
    , cell_cluster_key = centSO@ident
    , cluster_annot_key = cluster_annot_key
    , expr_lower_limit = -1
    , expr_upper_limit = 1) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in fetal sorted by cluster enrichment"
        , "\nPolioudakis fetal dropseq dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nMean of normalized centered variance scaled expression"
        , "\n")
    )
    ggsave(paste0(
        out_graph, "ASDlaura_EnrichmentSorted_HeatmapSetColWid_MeanNormCenterScale.png")
      , width = 10, height = 5
    )
  # Normalized
  plot_mean_expression_by_cluster_heatmap(
    expr_m = noCentExM
    , genes = genes
    , cell_cluster_key = centSO@ident
    , cluster_annot_key = cluster_annot_key
    , expr_lower_limit = -0.5
    , expr_upper_limit = 1.5
    , color_viridis = TRUE) +
    ggtitle(paste0(
      script_name
      , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
      , "\nLake et al. 2017 dataset"
      , "\nx-axis: Genes"
      , "\ny-axis: Clusters"
      , "\nMean of normalized mean expression"
      , "\n")
    )
    ggsave(paste0(
        out_graph
        , "ASDlaura_EnrichmentSorted_HeatmapSetColWid_MeanNorm.png")
      , width = 10, height = 5
    )
  # Lake
  # Normalized centered scaled
  cell_cluster_key <- gsub("_.*", "", colnames(lake_ex_df))
  names(cell_cluster_key) <- colnames(lake_ex_df)
  cluster_annot_key <- c("Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
    , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
    , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
  names(cluster_annot_key) <- cluster_annot_key
  plot_mean_expression_by_cluster_heatmap(
    expr_m = t(scale(t(lake_ex_df)))
    , genes = genes
    , cell_cluster_key = cell_cluster_key
    , cluster_annot_key = cluster_annot_key
    , expr_lower_limit = -1
    , expr_upper_limit = 1) +
    ggtitle(paste0(
      script_name
      , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
      , "\nLake et al. 2017 dataset"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nMean of normalized centered variance scaled expression"
      , "\n")
    )
    ggsave(paste0(
        out_graph, "ASDlaura_Lake_EnrichmentSorted_HeatmapSetColWid_MeanNormCenterScale.png")
      , width = 10, height = 5
    )
  # Normalized
  cell_cluster_key <- gsub("_.*", "", colnames(lake_ex_df))
  names(cell_cluster_key) <- colnames(lake_ex_df)
  cluster_annot_key <- c("Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
    , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
    , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
  names(cluster_annot_key) <- cluster_annot_key
  plot_mean_expression_by_cluster_heatmap(
    expr_m = t(scale(t(lake_ex_df)))
    , genes = genes
    , cell_cluster_key = cell_cluster_key
    , cluster_annot_key = cluster_annot_key
    , expr_lower_limit = -0.5
    , expr_upper_limit = 1.5
    , color_viridis = TRUE) +
    ggtitle(paste0(
      script_name
      , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
      , "\nLake et al. 2017 dataset"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nMean of normalized expression"
      , "\n")
    )
    ggsave(paste0(
        out_graph, "ASDlaura_Lake_EnrichmentSorted_HeatmapSetColWid_MeanNorm.png")
      , width = 10, height = 5
    )
}
################################################################################

### Plot heatmaps of ZNF559

plot_expression_by_cluster_heatmaps_znf559 <- function(){
  print("plot_expression_by_cluster_heatmaps_znf559")

  # Dropseq
  geneGroupDF <- data.frame(
    Gene = "ZNF559"
    , Group = ""
  )
  # Annotate cluster names
  cluster_annot_key <- c(
    "9" = "vRG"
    , "7" = "oRG"
    , "8" = "PgS"
    , "10" = "PgG2M"
    , "2" = "IP"
    , "0" = "ExN"
    , "1" = "ExM"
    , "4" = "ExCal"
    , "3" = "ExDp1"
    , "13" = "ExDp2"
    , "5" = "InSST"
    , "6" = "InCALB2"
    , "11" = "OPC"
    , "12" = "End"
    , "14" = "Per"
    , "15" = "Mic"
    , "16" = "NA"
  )
  cellID_clusterID_DF <- data.frame(centSO@ident)
  idx <- match(cellID_clusterID_DF[ ,1], names(cluster_annot_key))
  cellID_clusterID_DF$cluster <- cluster_annot_key[idx]
  cellID_clusterID <- cellID_clusterID_DF$cluster
  names(cellID_clusterID) <- row.names(cellID_clusterID_DF)
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID <- cellID_clusterID
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in fetal sorted by cluster enrichment"
        , "\nPolioudakis fetal dropseq dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression, mean centered, variance scaled"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ZNF559_EnrichmentSorted_HeatmapSetColWid_NormCenterScale.png")
    , width = 10, height = 3
  )

  # Lake
  cellID_clusterID <- gsub("_.*", "", colnames(lake_ex_df))
  names(cellID_clusterID) <- colnames(lake_ex_df)
  clusterOrder <- c("Ex1", "Ex2", "Ex3e", "Ex4", "Ex5b", "Ex6a","Ex6b"
    , "Ex8","In1a", "In1b", "In1c", "In3", "In4a", "In4b", "In6a", "In6b"
    , "In7", "In8", "OPC", "End", "Oli", "Per", "Mic", "Ast")
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = t(scale(t(lake_ex_df)))
    , cellID_clusterID = cellID_clusterID
    , clusterOrder = clusterOrder
    , clusters = clusterOrder
    , lowerLimit = -1.5
    , upperLimit = 1.5
  ) +
    ggtitle(paste0(
      script_name
        , "\n\nExpression of ASD genes from Laura in adult sorted by fetal Drop-seq cluster enrichment"
        , "\nLake et al. 2017 dataset"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression, mean centered, variance scaled"
        , "\n")
    )
  ggsave(paste0(
      out_graph
      , "ZNF559_Lake_EnrichmentSorted_HeatmapSetColWid_NormCenterScale.png")
    , width = 10, height = 3
  )
}
################################################################################

### Run

main_function <- function(){
  # plot_pie_charts()
  plot_expression_by_cluster_heatmaps()
  plot_expression_by_cluster_heatmaps_znf559()
}
main_function()
################################################################################
