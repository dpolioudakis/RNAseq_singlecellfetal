# Damon Polioudakis
# 2017-11-01
# Separate cycling vRG and oRG

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(reshape2)
require(tidyverse)
require(cowplot)
require(ggbeeswarm)
require(ggpubr)
source("Function_Library.R")
source("GGplot_Theme.R")

options(stringsAsFactors = FALSE)

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

# Cluster DE genes
cluster_de_df <- read.table(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Cell type specific genes
celltype_specific_df <- read.csv("../analysis/tables/Seurat_ClassMarkers/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClassMarkers_Markers_Class.csv")

# ATAC regulatory elements
atac_hic_re_df <- read.csv("../source/Gene_Lists/TorreUbieta_2018_TS2.csv"
  , header = TRUE)
# atac_nohic_re_df <- read.csv("../source/All_ATACcor_toTSS.csv", header = TRUE)
atac_nohic_re_df <- read.csv(
  "../source/All_ATACcor_toTSSwDiffAcc-SignifCor.csv"
  , header = TRUE)


# Nowakowski cluster DE table
nowakowski_de_df <- read.csv(
  "../nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)

# Lake 2017 cluster DE tables
lake_de_df <- read.csv(
  "../source/Lake_2018_TS3_Cluster_DE_Genes.csv"
  , skip = 4, header = TRUE, fill = TRUE
)

# biomaRt gene info
bm_df <- read.csv(
  "../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# house keeping genes
hk_genes_tb <- read_tsv("../source/Gene_Lists/hk_gene_list.txt")

## Variables
script_name <- "ATAC_REs.R"
out_graph <- "../analysis/graphs/ATAC_REs/20181031/ATAC_REs_"
out_table <- "../analysis/tables/ATAC_REs/20181031/ATAC_REs_"
out_data <- "../analysis/analyzed_data/20181031/ATAC_REs/ATAC_REs_"

## Output Directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_data), recursive = TRUE)
################################################################################

### Functions

clean_variable_names <- function(data){
  cleaned <- data %>%
    rename_all(
      funs(
        gsub("* ", "_", .) %>%
        gsub("\\.", "_", .) %>%
        gsub("\\(", "_", .) %>%
        gsub("\\)", "_", .) %>%
        gsub("\\+", "and", .) %>%
        gsub("#", "_number", .) %>%
        gsub("_$", "", .) %>%
        gsub("__", "_", .) %>%
        tolower
      )
    )
  return(cleaned)
}

classify_cells_by_type <- function(subset_de_df, class_cluster_idx){
  print("classify_cells_by_type")
  # Example class_cluster_idx input
  # $`Glia or support cells`
  # [1] "End"     "Per"     "Ast"     "Ast_Cer" "Oli"     "OPC"     "OPC_Cer"
  # [8] "Mic"
  # $`NA`
  # [1] "Gran"
  class_cluster_idx <- melt(class_cluster_idx)
  subset_de_df$Class <- class_cluster_idx$L1[
    match(subset_de_df$Cluster, class_cluster_idx$value)]
  return(subset_de_df)
}
################################################################################

### Clean and format data

hk_genes_tb <- clean_variable_names(hk_genes_tb) %>%
  add_column(housekeeping = TRUE)
################################################################################

### Cell type enriched REs

## Drop-seq

# Remove cluster 16 (low quality cells)
cluster_de_df <- cluster_de_df[cluster_de_df$Cluster != 16, ]

# Merge ATAC+HiC enhancers table and cell type enriched genes table by ensid
atac_hic_re_cluster_df <- merge(
  atac_hic_re_df, cluster_de_df
  , by.x = "ENSGID", by.y = "Ensembl"
)
# Merge ATAC without HiC enhancers table and cell type enriched genes table by ensid
atac_only_re_cluster_df <- merge(
  atac_nohic_re_df, cluster_de_df
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atac_hic_re_cluster_df
  , file = paste0(out_table, "ATACandHiC_CellTypeEnrichment.csv")
  , quote = FALSE
)
write.csv(atac_only_re_cluster_df
  , file = paste0(out_table, "ATAConly_CellTypeEnrichment.csv")
  , quote = FALSE
)
# Format for paper
paper_atac_DF <- atac_only_re_cluster_df
paper_atac_DF <- paper_atac_DF[c("ENSGID", "hgnc_symbol", "enhancerchr"
  , "enhancerstart", "enhancerend", "cor", "p.value", "p.fdr"
  , "gene_biotype", "promoterchr", "promoterstart", "promoterend", "Cluster"
  , "Log2_Fold_Change", "Pvalue", "FDR")]
names(paper_atac_DF)[14:16] <- c(
  "Cluster_enrichment_log2_fold_change"
  , "Cluster_enrichment_pvalue"
  , "Cluster_enrichment_FDR")
paper_atac_DF$Cluster_number <- paper_atac_DF$Cluster
# Cluster annotations
cluster_annot <- c(
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
  # "9" = "vRG"
  # , "7" = "oRG"
  # , "8" = "Cycling progenitor S phase"
  # , "10" = "Cycling progenitor G2/M phase"
  # , "2" = "IPC"
  # , "0" = "Excitatory neuron new born migrating"
  # , "1" = "Excitatory neuron"
  # , "4" = "Excitatory neuron (collosal)"
  # , "3" = "Deep layer excitatory neuron 1"
  # , "13" = "Deep layer excitatory neuron 2"
  # , "5" = "Interneuron (SST)"
  # , "6" = "Interneuron (CALB2)"
  # , "11" = "Oligodendrocyte precursor"
  # , "12" = "Endothelial"
  # , "14" = "Pericyte"
  # , "15" = "Microglia"
)
idx <- match(paper_atac_DF$Cluster_number, names(cluster_annot))
paper_atac_DF$Cluster <- cluster_annot[idx]
paper_atac_DF$Cluster_number <- factor(paper_atac_DF$Cluster_number
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
paper_atac_DF <- paper_atac_DF[order(paper_atac_DF$Cluster_number), ]
paper_atac_DF <- paper_atac_DF[, colnames(paper_atac_DF) != "Cluster_number"]
write.csv(paper_atac_DF
  , file = paste0(out_table, "ATAConly_CellTypeEnrichment_paper.csv")
  , quote = FALSE, row.names = FALSE
)

# Merge ATAC+HiC enhancers table and major cell type specific genes table by ensid
atacRE_celltypeDE_DF <- merge(
  atac_hic_re_df, celltype_specific_df
  , by.x = "ENSGID", by.y = "Ensembl"
)
# Merge ATAC without HiC enhancers table and major cell type specific genes table by ensid
atacRE_noHiC_celltypeDE_DF <- merge(
  atac_nohic_re_df, celltype_specific_df
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atacRE_celltypeDE_DF
  , file = paste0(out_table, "ATACandHiC_CellTypeSpecific.csv")
  , quote = FALSE
)
write.csv(atacRE_noHiC_celltypeDE_DF
  , file = paste0(out_table, "ATAConly_CellTypeSpecific.csv")
  , quote = FALSE
)

## Nowakowski
# Subset to higher enrichment
subset_de_df <- nowakowski_de_df[
  nowakowski_de_df$avg_diff > 0.75, c("gene", "avg_diff", "cluster")]
colnames(subset_de_df) <- c("Gene", "Enrichment", "Cluster")
# Classify cells by type
class_cluster_idx <- list(
  "Astrocyte" = "Astrocyte"
  , "Choroid" = "Choroid"
  , "Cycling progenitors" = c("IPC-div1", "RG-div1", "RG-div2")
  , "Endothelial" = "Endothelial"
  , "IP" = c("IPC-nEN1", "IPC-nEN2", "IPC-nEN3", "IPC-div2")
  , "Microglia" = "Microglia"
  , "OPC" = "OPC"
  , "oRG" = "oRG"
  , "tRG" = "tRG"
  , "vRG" = "vRG"
  , "NA" = c("MGE-div", "MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "MGE-RG1", "MGE-RG2"
    , "nIN1", "nIN2", "nIN3", "nIN4", "nIN5", "Glyc", "U1", "U2", "U3", "U4"
    , "Mural"
    )
  , "Excitatory" = c(
    "EN-PFC2", "EN-PFC3", "EN-V1-2", "EN-PFC1", "EN-V1-1", "EN-V1-3")
  , "Excitatory early" = c("nEN-early2", "nEN-early1")
  , "Excitatory late" = c("nEN-late")
  , "Interneuron" = c("IN-CTX-CGE", "IN-CTX-CGE2", "IN-CTX-MGE2", "IN-CTX-MGE1"
    , "IN-STR"
    )
)
subset_de_df <- classify_cells_by_type(subset_de_df, class_cluster_idx)
subset_de_df <- subset_de_df[subset_de_df$Class != "NA", ]
# Add ensembl ID
idx <- match(subset_de_df$Gene, bm_df$hgnc_symbol)
subset_de_df$Ensembl <- bm_df$ensembl_gene_id[idx]
# Merge with ATAC REs
atacRE_nowakowskiDE_DF <- merge(
  atac_nohic_re_df, subset_de_df
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atacRE_nowakowskiDE_DF
  , file = paste0(
    out_table, "ATAConly_Nowakowski_CellTypeEnrichment.csv")
  , quote = FALSE
)

# ## Lake
# # Subset to higher enrichment
# subset_de_df <- lake_de_df[
#   lake_de_df$"Average.Difference..log.fold.change." > 0.25
#   , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]
# colnames(subset_de_df) <- c("Gene", "Enrichment", "Cluster")
# # Classify cells by type
# class_cluster_idx <- list(
#   "Endothelial" = "End"
#   , "Pericyte" = "Per"
#   , "Astrocyte" = c("Ast", "Ast_Cer")
#   , "Oligodendrocyte" = "Oli"
#   , "OPC" = c("OPC", "OPC_Cer")
#   , "Microglia" = "Mic"
#   , "Granule" = "Gran"
#   , "Excitatory" = c("Ex1", "Ex2", "Ex3a", "Ex3b", "Ex3c", "Ex3d", "Ex3e"
#     , "Ex4", "Ex5a", "Ex5b", "Ex6a", "Ex6b", "Ex8"
#     )
#   , "Interneuron" = c("In1a", "In1b", "In1c", "In2", "In3", "In4a", "In4b"
#     , "In6a", "In6b", "In7", "In8"
#   )
#   , "Purkinje" = c("Purk1", "Purk2")
# )
# subset_de_df <- classify_cells_by_type(subset_de_df, class_cluster_idx)
# subset_de_df <- subset_de_df[subset_de_df$Class != "NA", ]
# # Add ensembl ID
# idx <- match(subset_de_df$Gene, bm_df$hgnc_symbol)
# subset_de_df$Ensembl <- bm_df$ensembl_gene_id[idx]
# # Merge with ATAC REs
# atacRE_lakeDE_DF <- merge(
#   atac_hic_re_df, subset_de_df
#   , by.x = "ENSGID", by.y = "Ensembl"
# )
# write.csv(atacRE_lakeDE_DF
#   , file = paste0(
#     out_table, "ATAConly_Lake_CellTypeEnrichment.csv")
#   , quote = FALSE
# )
################################################################################

### REs for cluster enriched versus house keeping genes

# Merge ATAC+HiC enhancers table and cell type enriched genes table by ensid
atac_hic_re_cluster_df <- full_join(
  x = atac_hic_re_df, y = cluster_de_df, by = c("ENSGID" = "Ensembl"))
# Merge ATAC without HiC enhancers table and cell type enriched genes table by ensid
atac_only_re_cluster_df <- full_join(
  x = atac_nohic_re_df, y = cluster_de_df, by = c("ENSGID" = "Ensembl"))

# merge atac only and atac + hic cluster enriched tables
atac_only_re_cluster_df <- add_column(
  atac_only_re_cluster_df, re_source = "atac only")
atac_hic_re_cluster_df <- add_column(
  atac_hic_re_cluster_df, re_source = "atac and hic")
atac_hic_re_cluster_df <- rename(atac_hic_re_cluster_df
  , p.value = cor.p.value
  , p.fdr = cor.p.fdr
  , padj = padj.diffacc
  , promoterchr = promoterpeakchr
  , promoterstart = promoterpeakstart
  , promoterend = promoterpeakend)
re_cluster_df <- bind_rows(atac_only_re_cluster_df, atac_hic_re_cluster_df)

# add housekeeping genes
re_cluster_df <- re_cluster_df %>%
  full_join(hk_genes_tb, by = c("hgnc_symbol" = "gene_name")) %>% as_tibble %>%
  # filter(housekeeping == TRUE & is.na(Cluster))
  # do not label genes "house keeping" that are cluster enriched
  mutate(housekeeping = if_else(
    housekeeping == TRUE & is.na(Cluster), "TRUE", "FALSE"))

# add gene category column (housekeeping, cluster enriched)
re_cluster_df <- re_cluster_df %>%
  mutate(gene_category = if_else(
    housekeeping == TRUE, "housekeeping", "other genes")) %>%
  mutate(gene_category = if_else(
    ! is.na(Cluster), "cluster enriched", gene_category)) %>%
  mutate(gene_category = if_else(
    is.na(gene_category), "other genes", gene_category))

# add GZ vs CP differential accessibility gene category column
# (housekeeping, cluster enriched, GZ > CP, CP > GZ)
re_cluster_df <- re_cluster_df %>%
  mutate(gene_category_diff_acc = if_else(
    housekeeping == TRUE, "housekeeping", "other genes")) %>%
  mutate(gene_category_diff_acc = if_else(
    ! is.na(Cluster), "cluster enriched", gene_category_diff_acc)) %>%
  mutate(gene_category_diff_acc = if_else(
      ! is.na(Cluster) &
      padj < 0.05 &
      log2FoldChange > 0
    , "GZ > CP & cluster enriched", gene_category_diff_acc)) %>%
  mutate(gene_category_diff_acc = if_else(
      ! is.na(Cluster) &
      padj < 0.05 &
      log2FoldChange < 0
    , "CP > GZ & cluster enriched", gene_category_diff_acc)) %>%
  mutate(gene_category_diff_acc = if_else(
    is.na(gene_category_diff_acc), "other genes", gene_category_diff_acc))

# Remove cluster 16
re_cluster_df <- filter(re_cluster_df, is.na(Cluster) | Cluster != 16)

# Reorder clusters
re_cluster_df <- re_cluster_df %>% mutate(Cluster = factor(re_cluster_df$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15,NA), exclude = NULL))

# calculate or add metrics
# Add biomart gene info
idx <- match(re_cluster_df$ENSGID, bm_df$ensembl_gene_id)
re_cluster_df <- cbind(re_cluster_df
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)
# Distance enhancer end to promoter start
re_cluster_df$distance_enhancer_promoter <- abs(
  re_cluster_df$enhancerend -
  re_cluster_df$promoterstart
)

# distance_enhancer_promoter
means <- re_cluster_df %>% as_tibble %>%
  filter(! is.na(re_source)) %>%
  group_by(gene_category,re_source) %>%
  summarize(
    mean_distance = round(mean(distance_enhancer_promoter, na.rm = TRUE), 1)
    , median_distance = round(median(distance_enhancer_promoter, na.rm = TRUE),1)
  ) %>%
  mutate(label = paste0(
    "mean:\n", mean_distance, "\nmedian:\n", median_distance))
re_cluster_df %>% as_tibble %>%
  filter(! is.na(re_source)) %>%
  ggplot(aes(
    x = gene_category, y = distance_enhancer_promoter, fill = gene_category)) +
    facet_wrap(~re_source) +
    geom_violin(fill = "lightgrey", color = "lightgrey") +
    geom_boxplot(width = 0.05, outlier.shape = NA) +
    geom_text(data = means, aes(
      x = gene_category, y = median_distance + 5e5, label = label)) +
    stat_compare_means(comparisons = list(
      c("cluster enriched", "housekeeping")
      , c("cluster enriched", "other genes")
      , c("housekeeping", "other genes"))) +
    ggplot_set_theme_publication +
    theme(legend.position = "none") +
    ggtitle(paste0(script_name
      , "\n\nATAC enhancers with cell type specific enrichment"
      , "\nDistance of enhancer end to promoter start"))
    ggsave(paste0(out_graph, "atac_hic_housekeeping_distance_enhancer_promoter_boxplot.pdf")
      , width = 9, height = 7)

# distance_enhancer_promoter diff acc
means <- re_cluster_df %>% as_tibble %>%
  filter(! is.na(re_source)) %>%
  group_by(gene_category_diff_acc,re_source) %>%
  summarize(
    mean_distance = round(mean(distance_enhancer_promoter, na.rm = TRUE), 1)
    , median_distance = round(median(distance_enhancer_promoter, na.rm = TRUE),1)
  ) %>%
  mutate(label = paste0(
    "mean:\n", mean_distance, "\nmedian:\n", median_distance))
re_cluster_df %>% as_tibble %>%
  filter(! is.na(re_source)) %>%
  ggplot(aes(x = gene_category_diff_acc, y = distance_enhancer_promoter
    , fill = gene_category_diff_acc)) +
    facet_wrap(~re_source) +
    geom_violin(fill = "lightgrey", color = "lightgrey") +
    geom_boxplot(width = 0.05, outlier.shape = NA) +
    geom_text(data = means, aes(
      x = gene_category_diff_acc, y = Inf, label = label), vjust = 1) +
    ggplot_set_theme_publication +
    theme(legend.position = "none"
      , axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(script_name
      , "\n\nATAC enhancers with cell type specific enrichment"
      , "\nDistance of enhancer end to promoter start"))
    ggsave(paste0(out_graph, "atac_hic_housekeeping_diff_acc_distance_enhancer_promoter_boxplot.pdf")
      , width = 9, height = 7)

plot_enhancers_per_gene_box_violin_plots <- function(
  enhancers_per_gene_df, re_source_filter){
  print("plot_enhancers_per_gene_box_violin_plots")
  # calculate means and medians
  mean_median_df <- enhancers_per_gene_df %>% as_tibble %>%
    filter(re_source == re_source_filter) %>%
    group_by(gene_category, re_source) %>%
    summarize(
      mean_number = round(mean(value, na.rm = TRUE), 1)
      , median_number = round(median(value, na.rm = TRUE), 1)
    ) %>%
    mutate(label = paste0(
      "mean: ", mean_number, "\nmedian: ", median_number))
  # select atac only or atac + hic data and plot
  enhancers_per_gene_df %>% filter(re_source == re_source_filter) %>%
    ggplot(aes(x = gene_category, y = value, fill = gene_category)) +
      geom_violin(fill = "lightgrey", color = "lightgrey") +
      geom_boxplot(width = 0.05, outlier.shape = NA) +
      # coord_cartesian(ylim = c(0, 10)) +
      geom_text(data = mean_median_df, aes(
        x = gene_category, y = Inf, label = label), vjust = 1) +
      ggplot_set_theme_publication +
      theme(legend.position = "none") +
      ylab("Enhancers per gene")
}

# Enhancers per gene
enhancers_per_gene_df <- melt(
  with(re_cluster_df, table(ENSGID, re_source, gene_category)))
enhancers_per_gene_df <- enhancers_per_gene_df[
  enhancers_per_gene_df$value != 0, ]
# atac and hic
atac_hic_gg <- plot_enhancers_per_gene_box_violin_plots(
  enhancers_per_gene_df = enhancers_per_gene_df
  , re_source_filter = "atac and hic")
atac_hic_gg <- atac_hic_gg + stat_compare_means(comparisons = list(
    c("cluster enriched", "housekeeping")
    , c("cluster enriched", "other genes")
    , c("housekeeping", "other genes"))
  , label.y = c(9.5, 8.5, 7.5), tip.length = 0.003) +
  coord_cartesian(ylim = c(0, 20)) +
  ggtitle("atac and hic")
# atac only
atac_only_gg <- plot_enhancers_per_gene_box_violin_plots(
  enhancers_per_gene_df = enhancers_per_gene_df, re_source = "atac only")
atac_only_gg <- atac_only_gg + stat_compare_means(comparisons = list(
    c("cluster enriched", "housekeeping")
    , c("cluster enriched", "other genes")
    , c("housekeeping", "other genes"))
  , label.y = c(55, 50, 45), tip.length = 0.003) +
  coord_cartesian(ylim = c(0, 75)) +
  ggtitle("atac and hic")
# combine plots
Plot_Grid(list(atac_hic_gg, atac_only_gg), rel_height = 0.3, label_size = 10
  , title = paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
  ggsave(paste0(out_graph, "atac_hic_housekeeping_enhancers_per_gene_boxplot.pdf")
    , width = 9, height = 6)

# Enhancers per gene diff accessibility
enhancers_per_gene_df <- melt(
  with(re_cluster_df, table(ENSGID, re_source, gene_category_diff_acc)))
enhancers_per_gene_df <- enhancers_per_gene_df[
  enhancers_per_gene_df$value != 0, ]
enhancers_per_gene_df <- rename(enhancers_per_gene_df
  , gene_category = gene_category_diff_acc)
# atac and hic
atac_hic_gg <- plot_enhancers_per_gene_box_violin_plots(
  enhancers_per_gene_df = enhancers_per_gene_df
  , re_source_filter = "atac and hic")
atac_hic_gg <- atac_hic_gg + coord_cartesian(ylim = c(0, 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("atac and hic")
# atac only
atac_only_gg <- plot_enhancers_per_gene_box_violin_plots(
  enhancers_per_gene_df = enhancers_per_gene_df, re_source = "atac only")
atac_only_gg <- atac_only_gg + coord_cartesian(ylim = c(0, 75)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("atac and hic")
# combine plots
Plot_Grid(list(atac_hic_gg, atac_only_gg), rel_height = 0.3, label_size = 10
  , title = paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(out_graph
    , "atac_hic_housekeeping_diff_acc_enhancers_per_gene_boxplot.pdf")
  , width = 11, height = 6)
################################################################################

### Plots ATAC and HiC

# Add biomart gene info
idx <- match(atac_hic_re_df$ENSGID, bm_df$ensembl_gene_id)
atac_hic_re_df <- cbind(atac_hic_re_df
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)
# Add biomart gene info
idx <- match(atac_hic_re_cluster_df$ENSGID, bm_df$ensembl_gene_id)
atac_hic_re_cluster_df <- cbind(atac_hic_re_cluster_df
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)

# Reorder clusters
atac_hic_re_cluster_df$Cluster <- factor(atac_hic_re_cluster_df$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Remove cluster 16
atac_hic_re_cluster_df <- atac_hic_re_cluster_df[
  ! is.na(atac_hic_re_cluster_df$Cluster), ]
# Distance enhancer end to promoter start
atac_hic_re_cluster_df$distance_enhancer_promoter <- abs(
  atac_hic_re_cluster_df$enhancerend -
  atac_hic_re_cluster_df$promoterpeakstart
)
atac_hic_re_cluster_df$Enhancer_Size <- abs(
  atac_hic_re_cluster_df$enhancerend -
  atac_hic_re_cluster_df$enhancerstart
)

# Plots
# distance_enhancer_promoter
ggplot(atac_hic_re_cluster_df, aes(
  x = Cluster, y = distance_enhancer_promoter, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  geom_quasirandom(size = 0.1) +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAChiC_distance_enhancer_promoter_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atac_hic_re_cluster_df, aes(
  x = gene_biotype, y = distance_enhancer_promoter, fill = gene_biotype)) +
  geom_boxplot(outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  geom_quasirandom(size = 0.1) +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAChiC_distance_enhancer_promoter_By_Biotype_boxplot.pdf")
  , width = 6, height = 5)
ggplot(atac_hic_re_cluster_df, aes(x = distance_enhancer_promoter)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAChiC_distance_enhancer_promoter_density.pdf"))
# Enhancer_Size
ggplot(atac_hic_re_cluster_df, aes(
  x = Cluster, y = Enhancer_Size, fill = Cluster)) +
  # facet_wrap(~gene_biotype) +
  geom_quasirandom(size = 0.1) +
  ggplot_set_theme_publication +
  geom_boxplot(outlier.shape = NA) +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(out_graph, "ATAChiC_Enhancer_Size_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atac_hic_re_cluster_df, aes(x = Enhancer_Size)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(out_graph, "ATAChiC_Enhancer_Size_density.pdf")
  , width = 7, height = 5)


# Enhancers per gene boxplot
ggplot(ggDF, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = 0.1) +
  # coord_cartesian(ylim = c(0, 10)) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ylab("Enhancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(out_graph, "ATAChiC_Enhancers_Per_Gene_boxplot.pdf")
  , width = 7, height = 5)
# Enhancers per gene density plot
ggplot(ggDF, aes(x = value)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(out_graph, "ATAChiC_Enhancer_Per_Gene_density.pdf")
  , width = 7, height = 5)

# Enhancers per gene versus CDS length
ggplot(ggDF, aes(x = value, y = cds_length, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("CDS length") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus CDS length"))
ggsave(paste0(out_graph, "ATAChiC_Enhancers_Per_Gene_vs_CDS_Length.pdf")
  , width = 7, height = 5)

# Enhancers per gene versus GC content
ggplot(ggDF, aes(x = value, y = percentage_gc_content, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("GC content") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus gene GC content"))
ggsave(paste0(out_graph, "ATAChiC_Enhancers_Per_Gene_vs_GC_Content.pdf")
  , width = 7, height = 5)
################################################################################

### Plots ATAC only

# Add biomart gene info
idx <- match(atac_nohic_re_df$ENSGID, bm_df$ensembl_gene_id)
atac_nohic_re_df <- cbind(atac_nohic_re_df
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)
# Add biomart gene info
idx <- match(atac_only_re_cluster_df$ENSGID, bm_df$ensembl_gene_id)
atac_only_re_cluster_df <- cbind(atac_only_re_cluster_df
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)

# Reorder clusters
atac_only_re_cluster_df$Cluster <- factor(atac_only_re_cluster_df$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Remove cluster 16
atac_only_re_cluster_df <- atac_only_re_cluster_df[
  ! is.na(atac_only_re_cluster_df$Cluster), ]
# Distance enhancer end to promoter start
atac_only_re_cluster_df$distance_enhancer_promoter <- abs(
  atac_only_re_cluster_df$enhancerend -
  atac_only_re_cluster_df$promoterstart
)
atac_only_re_cluster_df$Enhancer_Size <- abs(
  atac_only_re_cluster_df$enhancerend -
  atac_only_re_cluster_df$enhancerstart
)

# Remove enhancer entries that are actually promoters
ss_atac_only_re_cluster_df <- atac_only_re_cluster_df[
  atac_only_re_cluster_df$enhancerstart != atac_only_re_cluster_df$promoterstart
  , ]

# Plots
# distance_enhancer_promoter
ggplot(atac_only_re_cluster_df, aes(
  x = Cluster, y = distance_enhancer_promoter, fill = Cluster)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAConly_distance_enhancer_promoter_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atac_only_re_cluster_df, aes(
  x = gene_biotype, y = distance_enhancer_promoter, fill = gene_biotype)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAConly_distance_enhancer_promoter_By_Biotype_boxplot.pdf")
  , width = 6, height = 5)
ggplot(atac_only_re_cluster_df, aes(x = distance_enhancer_promoter)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(out_graph, "ATAConly_distance_enhancer_promoter_density.pdf"))
# Enhancer_Size
ggplot(atac_only_re_cluster_df, aes(
  x = Cluster, y = Enhancer_Size, fill = Cluster)) +
  # facet_wrap(~gene_biotype) +
  coord_cartesian(ylim = c(0, 5000)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(out_graph, "ATAConly_Enhancer_Size_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atac_only_re_cluster_df, aes(x = Enhancer_Size)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(out_graph, "ATAConly_Enhancer_Size_density.pdf")
  , width = 7, height = 5)


# Enhancers per gene
ggDF <- melt(with(atac_only_re_cluster_df, table(ENSGID, Cluster)))
# Add biomart gene info
idx <- match(ggDF$ENSGID, bm_df$ensembl_gene_id)
ggDF <- cbind(ggDF
  , bm_df[idx, c("percentage_gc_content", "cds_length")]
)
ggDF <- ggDF[ggDF$value != 0, ]
ggDF$Cluster <- factor(ggDF$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Enhancers per gene boxplot
ggplot(ggDF, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 50)) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ylab("Enhancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(out_graph, "ATAConly_Enhancers_Per_Gene_boxplot.pdf")
  , width = 7, height = 5)
# Enhancers per gene density plot
ggplot(ggDF, aes(x = value)) +
  geom_histogram(binwidth = 1) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(out_graph, "ATAConly_Enhancer_Per_Gene_histogram.pdf")
  , width = 7, height = 5)

# Enhancers per gene versus CDS length
sapply(split(ggDF, ggDF$Cluster), function(subset_ggDF){
  cor(subset_ggDF$cds_length, subset_ggDF$value, use = "complete")
})
ggplot(ggDF, aes(x = value, y = cds_length, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  coord_cartesian(ylim = c(0, 10000)) +
  theme(legend.position = "none") +
  ylab("CDS length") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus CDS length"))
ggsave(paste0(out_graph, "ATAConly_Enhancers_Per_Gene_vs_CDS_Length.pdf")
  , width = 7, height = 6)

# Enhancers per gene versus GC content
ggplot(ggDF, aes(x = value, y = percentage_gc_content, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("GC content") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(script_name
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus gene GC content"))
ggsave(paste0(out_graph, "ATAConly_Enhancers_Per_Gene_vs_GC_Content.pdf")
  , width = 7, height = 6)
################################################################################
