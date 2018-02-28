# Damon Polioudakis
# 2017-04-05
# Pooled versus bulk
# R 3.4.0
################################################################################

rm(list = ls())
sessionInfo()
set.seed(27)

require(Seurat)
require(biomaRt)
require(gProfileR)
require(reshape2)
require(eulerr)
require(gridExtra)
require(cqn)
require(cowplot)
source("Function_Library.R")

# Processed data
load("../analysis/analyzed_data/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_Processed_Data.RData")

# biomart
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Gene lengths and GC content for Union Exon model to normalize bulk by length
load("../source/ENSEMBLhg19_UnionAnno.rda")

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

# # Markers from Seurat
# mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_Top20.txt"
#   , header = TRUE)

# Brain expressed genes from Vivek
brain_exp_DF <- read.csv("../source/BrainExpressed_GTex.csv", header = FALSE)

## Variables
graphCodeTitle <- "Pooled_Vs_Bulk.R"
outGraph <- "../analysis/graphs/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outGraph), recursive = TRUE)
outTable <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outTable), recursive = TRUE)
outTableGpro <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler/Pooled_Vs_Bulk_gprofiler_"
dir.create("../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler", recursive = TRUE)
outTableGproIds <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler_goIDs/Pooled_Vs_Bulk_gprofiler_goIDs_"
dir.create("../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler_goIDs", recursive = TRUE)
outAnalysis <- "../analysis/analyzed_data/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outAnalysis), recursive = TRUE)

## ggplot2 theme
## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Plots

Venn_Diagram_3 <- function(
  v1, v2, v3, fill = c("#e41a1c", "#377eb8", "#4daf4a"), ...) {
  print(length(intersect(v1, intersect(v2, v3))))
  venn <- euler(c(
    "A" = length(v1)
    , "B" = length(v2)
    , "C" = length(v3)
    , "A&B" = length(intersect(v1, v2))
    , "A&C" = length(intersect(v1, v3))
    , "B&C" = length(intersect(v2, v3))
    , "A&B&C" = length(intersect(v1, intersect(v2, v3)))))
  plot(venn, counts = TRUE, fill = fill, ...)
}

Venn_Diagram_3_Wrapper <- function(
  datasets
  , bias_direction
  , fill = c("#e41a1c", "#377eb8", "#4daf4a")
  , title = paste0(paste(datasets, collapse = ", "), ": ", bias_direction)
  , ...){

  ensembl_gene_id_L <- lapply(datasets, function(dataset){
    print(dataset)
    ensembl_gene_id <- row.names(bias_flag_DF)[
      bias_flag_DF[[dataset]] == bias_direction]
    return(ensembl_gene_id)
  })
  print(head(ensembl_gene_id_L[[1]]))
  Venn_Diagram_3(
    v1 = ensembl_gene_id_L[[1]]
    , v2 = ensembl_gene_id_L[[2]]
    , v3 = ensembl_gene_id_L[[3]]
    , labels = datasets
    , fill = fill
    , main = paste0(graphCodeTitle
      , "\n\nIntersection of single-cell or bulk RNA-seq biased gene sets"
      , "\n", title)
    , ...
  )
}

## Correlation to bulk as pool size increases
gg <- ggplot(pool_corr_DF, aes(
  x = Number_of_cells, y = Spearman, color = Dataset_label)) +
  facet_wrap(~Dataset_label, ncol = 2, scales = "free") +
  geom_line() +
  xlab("Number of cells in pool") +
  ylab("Spearman") +
  ylim(0, 1) +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCorrelation of different size pools to bulk"
    , "\nGZ and CP fetal human brain"
    , "\n"
  ))
ggsave(paste0(outGraph, "Pool_Size_Spearman.pdf"))
gg + ggplot_set_theme_publication
ggsave(paste0(outGraph, "Pool_Size_Spearman_paper.pdf"))

## Number of genes detected as pool size increases
ggDF <- melt(number_gene_detected_DF
  , measure.vars = c("Detection_threshold_0"
  , "Detection_threshold_1", "Detection_threshold_2"
  , "Detection_threshold_10")
)
# Plot as grid
ggplot(ggDF, aes(x = Number_of_cells, y = value, col = variable)) +
  facet_wrap(~Dataset, ncol = 3, scales = "free_x") +
  geom_line() +
  xlab("Number of cells in pool") +
  ylab("Number of genes detected") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nNumber genes detected as pool size increases"
    , "\nGZ and CP fetal human brain"
    , "\n"
  ))
ggsave(paste0(outGraph, "Pool_Size_Genes_Detected.pdf"), height = 7, width = 11)

## Scatter plots of pooled vs bulk with fit and spearman
# Format for ggplot
df1 <- melt(mean_cpm_DF[
      c("Bulk", "Bulk_GZ", "Dropseq", "Fluidigm_HT", "Fluidigm_LT")
      ]
    , measure.vars = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
  )
df2 <- melt(mean_cpm_DF[
    c("Bulk", "Bulk_GZ", "Dropseq_GZ", "Pollen", "Fluidigm_HT_GZ")
  ]
  , measure.vars = c("Dropseq_GZ", "Pollen", "Fluidigm_HT_GZ")
)
df2$Bulk <- df2$Bulk_GZ
ggDF <- rbind(df1, df2)
# Spearman
sprL <- lapply(split(ggDF, ggDF$variable), function(df) {
  cor(df$Bulk, df$value, method = "spearman")})
sprL <- lapply(sprL, round, 2)
ggDF$value <- log(ggDF$value + 1, 2)
ggDF$Bulk <- log(ggDF$Bulk + 1, 2)
labelsDatasetSpr <- c(paste0(unique(ggDF$variable), "\nSpearman: ", sprL))
names(labelsDatasetSpr) <- unique(ggDF$variable)
# ggplot
gg <- ggplot(ggDF, aes(x = value, y = Bulk, color = variable)) +
  facet_wrap("variable", ncol = 2
    , labeller = labeller(variable = labelsDatasetSpr), scales = "free") +
  geom_point(alpha = 0.15, shape = 1, size = 0.25) +
  xlim(c(0, 15)) +
  ylim(c(0, 15)) +
  stat_smooth(col = "black") +
  ylab("Bulk: log2(mean CPM + 1)") +
  xlab("Pooled: log2(mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPooled scRNA-seq vs Bulk RNA-seq"
    , "\n"
    , "\nMean of CPM across samples and pools"))
ggsave(paste0(outGraph, "Scatter.png"), width = 6.5, height = 9)
gg + ggplot_set_theme_publication
ggsave(paste0(outGraph, "Scatter_paper.png"), width = 6.5, height = 9)

## Table of number of genes biased toward bulk or single-cell
dfl1 <- lapply(names(bias_flag_DF), function(name){
  y <- bias_flag_DF[name]
  df1 <- data.frame(table(y))
  row.names(df1) <- df1$y
  df1 <- df1[ -1]
  names(df1) <- name
  return(df1)
})
df1 <- do.call("cbind", dfl1)
write.csv(df, file = paste0(outTable, "_NumberGenesSubset.csv"), quote = FALSE)

## Plot subset scatter plots and density plots
ldf1 <- lapply(c("Dropseq", "Fluidigm_HT", "Fluidigm_LT"), function(name){
  data.frame(
    Mean_Pooled = log2_mean_cpm_DF[[name]]
    , Mean_Bulk = log2_mean_cpm_DF$Bulk
    , Subset = bias_flag_DF[[name]]
    , Dataset = name
  )
})
ldf2 <- lapply(c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen"), function(name){
  data.frame(
    Mean_Pooled = log2_mean_cpm_DF[[name]]
    , Mean_Bulk = log2_mean_cpm_DF$Bulk_GZ
    , Subset = bias_flag_DF[[name]]
    , Dataset = name
  )
})
ggDF <- rbind(do.call("rbind", ldf1)
  , do.call("rbind", ldf2)
)
# scatter
gg <- ggplot(ggDF, aes(x = Mean_Pooled, y = Mean_Bulk)) +
  facet_wrap("Dataset", ncol = 2, scales = "free") +
  geom_point(alpha = 0.15, shape = 16, size = 0.25, aes(color = Subset)) +
  scale_colour_manual(values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  xlim(c(0, 15)) +
  ylim(c(0, 15)) +
  ylab("Bulk: log2(Mean CPM + 1)") +
  xlab("Pooled: log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSubsetting pooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain GZ and CP"
    , "\nMean of CPM across samples and pools"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"))
ggsave(paste0(outGraph, "ScatterSubset.png"), width = 6, height = 9)
gg + ggplot_set_theme_publication
ggsave(paste0(outGraph, "ScatterSubset_paper.png"), width = 6, height = 9)
# density plot
ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  facet_grid(variable ~ Dataset) +
  geom_density(aes(col = Subset)) +
  scale_colour_manual(values = c("grey","#1f78b4", "#e31a1c")) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Density") +
  xlab("log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSubsetting pooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain GZ and CP"
    , "\nMean of CPM across samples and pools"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"))
ggsave(paste0(outGraph, "DensitySubset.pdf"), width = 14, height = 8)

## Venn diagrams of gene intersections high in bulk or high in pooled
pdf(paste0(outGraph, "VennDiagrams.pdf"))
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
    , bias_direction = "Bulk high"
    , fill = c("#e41a1c", "#377eb8", "#984ea3")
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
    , bias_direction = "Pool high"
    , fill = c("#e41a1c", "#377eb8", "#984ea3")
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen")
    , bias_direction = "Bulk high"
    , fill = c("#fb9a99", "#a6cee3", "#b2df8a")
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen")
    , bias_direction = "Pool high"
    , fill = c("#fb9a99", "#a6cee3", "#b2df8a")
  )
dev.off()


Ttest_Table <- function(data_table, facet_col, group_col, value_col){
  print("T-test...")
  pvals_M <- sapply(split(data_table, data_table[[facet_col]]), function(df1) {
    print(head(df1))
    combinations <- combn(unique(df1[[group_col]]), 2)
    pvals_M <- matrix(NA, 1, ncol(combinations))
    for(i in 1:ncol(combinations)){
      results <- t.test(
        df1[[value_col]][df1$value == combinations[1,i]]
        , df1[[value_col]][df1$value == combinations[2,i]]
      )
      pvals_M[1,i] <- results$p.value
    }
    colnames(pvals_M) <- paste(combinations[1, ], "vs", combinations[2, ], "\npvalue")
    pvals_M <- signif(pvals_M, 2)
    return(as.data.frame(pvals_M))
  })
  pvals_M <- data.frame(pvals_M)
  pvals_M <- apply(pvals_M, 2, unlist)
  return(pvals_M)
}

## CDS length of biased gene sets
# Format for ggplot
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
idx <- match(ggDF$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
ggDF$cds_length <- mean_cpm_DF$cds_length[idx]
ggDF <- ggDF[! is.na(ggDF$cds_length), ]
ggDF$value[ggDF$value == "None"] <- "All other genes"
ggDF$value <- factor(ggDF$value
  , levels = c("Bulk high", "All other genes", "Pool high")
)
# Plot
ggplot(ggDF, aes(x = variable, y = cds_length, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 5000)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("CDS length") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "Length_Boxplot.pdf"), width = 8, height = 5)
# Subset for paper
ggDF <- ggDF[ggDF$variable %in% c("GZ_Datasets", "GZCP_Datasets"), ]
ggDF <- droplevels(ggDF)
# T test
pval_DF <- Ttest_Table(
  data_table = ggDF, facet_col = "variable", group_col = "value"
  , value_col = "cds_length"
)
ggt <- tableGrob(pval_DF)
gg <- ggplot(ggDF, aes(x = value, y = cds_length, fill = value)) +
  facet_wrap(~variable, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  # T test for plot
  stat_compare_means(cds_length ~ value, data = ggDF, group.by = "variable"
    , comparisons = list(
      c("Bulk high", "All other genes")
      , c("Bulk high", "Pool high")
      , c("All other genes", "Pool high")
    )
    , method = "t.test", p.adjust.method = "none"
    , label.y = c(8500, 9500, 7500), tip.length = 0.003) +
  coord_cartesian(ylim = c(0, 10000)) +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("CDS length") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
# Combine plot and table of p-values
plot_grid(gg, ggt, ncol = 1, rel_heights = c(1,0.5))
ggsave(paste0(outGraph, "Length_Boxplot_paper.pdf"), width = 6, height = 8)

# CDS length for each dataset
Plot_CDS_Length_Histogram <- function(exM, title, geneIdType){
  cGene <- rowSums(exM)
  ggDF <- data.frame(GENE = names(cGene), COUNTS = cGene)
  ggDF <- merge(ggDF, mean_cpm_DF[c("ensembl_gene_id", "cds_length")]
    , by.x = "GENE", by.y = "ensembl_gene_id")
  ggDF <- data.frame(CDS_LENGTH = rep(ggDF$cds_length, ggDF$COUNTS))
  gg <- ggplot(ggDF, aes(x = CDS_LENGTH)) +
    geom_histogram(binwidth = 500) +
    # coord_cartesian(xlim = c(0, 10000)) +
    # scale_x_continuous(limits = c(0, 10000)) +
    xlab("CDS length") +
    ylab("Counts") +
    ggtitle(title)
  return(gg)
}
ggL <- list(
  Plot_CDS_Length_Histogram(bulk_ex_DF, "Bulk RNA-seq")
  , Plot_CDS_Length_Histogram(ds_ex_DF, "Drop-seq")
  , Plot_CDS_Length_Histogram(gz_df_ex_DF, "Drop-seq GZ")
  , Plot_CDS_Length_Histogram(fldm_HT_ex_DF, "Fluidigm HT")
  , Plot_CDS_Length_Histogram(gz_fldm_HT_ex_DF, "Fluidigm HT GZ")
  , Plot_CDS_Length_Histogram(fldm_LT_ex_DF, "Fluidigm LT")
  , Plot_CDS_Length_Histogram(pollen_ex_DF, "Pollen")
)
Plot_Grid(ggL, ncol = 3, rel_height = 0.1, title = paste0(graphCodeTitle
  , "\n"
  , "\nHistogram of raw counts versus CDS length"
  , "\n")
)
ggsave(paste0(outGraph, "LengthCounts_Histogram.png"), width = 8, height = 20)

## GC content for each dataset
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
idx <- match(ggDF$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
ggDF$percentage_gc_content <- mean_cpm_DF$percentage_gc_content[idx]
ggDF <- ggDF[! is.na(ggDF$percentage_gc_content), ]
ggDF$value[ggDF$value == "None"] <- "All other genes"
ggDF$value <- factor(ggDF$value
  , levels = c("Bulk high", "All other genes", "Pool high")
)
# Plot
ggplot(ggDF, aes(x = variable, y = percentage_gc_content, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 100)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("Percent GC") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercent GC content of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "GC_Boxplot.pdf"), width = 8, height = 5)
# Subset for paper
ggDF <- ggDF[ggDF$variable %in% c("GZ_Datasets", "GZCP_Datasets"), ]
ggDF <- droplevels(ggDF)
# T test
pval_DF <- Ttest_Table(
  data_table = ggDF, facet_col = "variable", group_col = "value"
  , value_col = "percentage_gc_content"
)
ggt <- tableGrob(pval_DF)
gg <- ggplot(ggDF, aes(x = variable, y = percentage_gc_content, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 100)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggplot_set_theme_publication +
  xlab("Dataset") +
  ylab("Percent GC") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercent GC content of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
# Combine plot and table of p-values
plot_grid(gg, ggt, ncol = 1, rel_widths = c(1,0.5))
ggsave(paste0(outGraph, "GC_Boxplot_paper.pdf"), width = 5, height = 8)

## Biotype
# Format for ggplot2 and calculate percent of each biotype in each subset
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
idx <- match(ggDF$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
ggDF$gene_biotype <- mean_cpm_DF$gene_biotype[idx]
# # Remove TR_ and IG_ biotypes
# ggDF <- ggDF[-grep("TR_|IG_", ggDF$gene_biotype), ]
ggDF <- aggregate(gene_biotype~variable+value, ggDF, table)
# Calculate percentage
total <- rowSums(ggDF$gene_biotype)
percent_DF <- (ggDF$gene_biotype/total) * 100
# Format for ggplot
colnames(ggDF) <- c("Dataset", "Subset")
ggDF <- cbind(ggDF[ ,1:2], percent_DF)
ggDF <- melt(ggDF, measure.vars = 3:ncol(ggDF))
# ggplot
ggplot(ggDF, aes(x = variable, y = value, fill = Subset)) +
  facet_wrap(~Dataset) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_bar(aes(fill = SUBSET), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_flip() +
  theme(text = element_text(size = 12)) +
  xlab("Gene biotypes") +
  ylab("Percent of genes") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nFrequency of gene biotypes"
    , "\npooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"
    , "\nPercent of biotype relative to all genes in subset"))
ggsave(paste0(outGraph, "Biotypes.pdf"), width = 9.5, height = 20)
# Plot for paper
ggDF <- ggDF[ggDF$Dataset %in% c("GZ_Datasets", "GZCP_Datasets"), ]
ggplot(ggDF, aes(x = variable, y = value, fill = Subset)) +
  facet_grid(Dataset~Subset) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_bar(aes(fill = SUBSET), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_flip() +
  theme(text = element_text(size = 12)) +
  xlab("Gene biotypes") +
  ylab("Percent of genes") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nFrequency of gene biotypes"
    , "\npooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"
    , "\nPercent of biotype relative to all genes in subset"))
ggsave(paste0(outGraph, "Biotypes_paper.pdf"), width = 9.5, height = 12)

## Percent MT
# Format
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
idx <- match(ggDF$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
ggDF$hgnc_symbol <- mean_cpm_DF$hgnc_symbol[idx]
# Percentage of MT genes
ggDF$MT <- 0
ggDF$MT[grepl("^MT-", ggDF$hgnc_symbol)] <- 1
ggDF <- dcast(ggDF, variable + value ~ MT)
ggDF$Percent_MT <- (ggDF[ ,4] / (ggDF[,3] + ggDF[,4])) * 100
ggplot(ggDF, aes(x = value, y = Percent_MT, fill = value)) +
  facet_wrap(~variable, scale = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Gene subset") +
  ylab("Percent MT genes") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercentage of MT genes in genes biased towards single-cell or bulk RNA-seq"
    , "\n"))
ggsave(paste0(outGraph, "MT_Percent.pdf"), width = 7, height = 9)

## Percentage of brain expressed genes
# Format
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
# Percentage of brain expressed genes
ggDF$Brain_expressed <- 0
ggDF$Brain_expressed[ggDF$ensembl_gene_id %in% brain_exp_DF$V1] <- 1
ggDF <- dcast(ggDF, variable + value ~ Brain_expressed)
ggDF$Percent_brain_expressed <- (ggDF[ ,4] / (ggDF[,3] + ggDF[,4])) * 100
# ggplot
gg <- ggplot(ggDF, aes(x = value, y = Percent_brain_expressed, fill = value)) +
  facet_wrap(~variable, scale = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Gene subset") +
  ylab("Percent brain expressed genes") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPercentage of brain expressed genes in genes biased towards"
    , "\nsingle-cell or bulk RNA-seq"
    , "\n"))
ggsave(paste0(outGraph, "BrainExpressed_Percent.pdf"), width = 7, height = 9)
gg + ggplot_set_theme_publication
ggsave(paste0(outGraph, "BrainExpressed_Percent.pdf"), width = 7, height = 9)

## HGNC in Drop-seq 0
df1 <- bias_flag_DF[ ,"Dropseq_0", drop = FALSE]
df1$ensembl_gene_id <- row.names(df1)
idx <- match(df1$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
df1$hgnc_symbol <- mean_cpm_DF$hgnc_symbol[idx]
df1$hgnc_symbol_TF <- ! df1$hgnc_symbol == ""
df1 <- dcast(df1, Dropseq_0~hgnc_symbol_TF, length)
names(df1)[2:3] <- c("HGNC_Symbol", "No_HGNC_Symbol")
write.csv(df, paste0(outTable, "Dropseq0_Number_HGNCsymbol.csv"))

## GO analysis
# Write out for REVIGO
lapply(names(gpro_results_DFLL), function(name1) {
  print(name1)
  lapply(names(gpro_results_DFLL[[name1]]), function(name2){
    print(name2)
    # Full gprofiler table
    gpro_results_DF <- gpro_results_DFLL[[name1]][[name2]]
    write.csv(gpro_results_DF
      , file = paste0(outTableGpro, "GO_", name1, "_", name2, ".csv")
      , quote = FALSE, row.names = FALSE)
    # Term IDs for inputing into REVIGO
    go_term_id <- gpro_results_DF$term.id[
      gpro_results_DF$domain %in% c("BP", "MF", "CC")]
    write.table(go_term_id
      , paste0(outTableGproIds, "GO_TermID_", name1, "_", name2, ".txt")
      , quote = FALSE, row.names = FALSE)
  })
})

## GO term pvalue barplot
df1 <- gpro_results_DFLL[["GZCP_Datasets"]][["Bulk high"]]
df1$Dataset <- "GZCP_Datasets"
df2 <- gpro_results_DFLL[["GZ_Datasets"]][["Bulk high"]]
df2$Dataset <- "GZ_Datasets"
ggDF <- rbind(df1, df2)
ggDF <- ggDF[c("p.value", "term.name", "Dataset")]
ggDF <- ggDF[
  ggDF$term.name %in% c("cell adhesion"
    , "generation of neurons"
    , "axon development"
    , "Epigenetic regulation of gene expression"
    , "neuron projection development")
  , ]
ggDF$term.name <- factor(ggDF$term.name
  , levels = unique(ggDF[order(-ggDF$p.value), ]$term.name)
)
ggplot(ggDF, aes(x = term.name, y = -log(p.value, 10), fill = Dataset)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Dataset, scale = "free") +
  # scale_fill_manual(values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  geom_hline(yintercept = -log(0.5, 10), color = "red") +
  coord_flip() +
  xlab("GO terms") +
  ylab("-log10(p-value)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSelected GO terms from bulk biased gene lists"
    , "\nGZ CP single-cell datasets (Drop-seq, Fluidigm HT, Fluidigm LT)"
    , "\nGZ single-cell datasets (Drop-seq GZ, Fluidigm HT GZ, Pollen)"
    , "\n"))
ggsave(paste0(outGraph, "GO_paper.pdf"), width = 14, height = 6)

## Marker expression
df1 <- ratio_log2_mean_cpm_DF[ ,"Dropseq", drop = FALSE]
df1$ensembl_gene_id <- row.names(df1)
idx <- match(df1$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
df1$hgnc_symbol <- mean_cpm_DF$hgnc_symbol[idx]
df1$Grouping <- kmDF$Grouping[match(df1$hgnc_symbol, kmDF$Gene.Symbol)]
df1$Grouping <- factor(df1$Grouping
  , levels = c("Non-marker gene", as.character(unique(kmDF$Grouping))))
df1$Grouping[is.na(df1$Grouping)] <- "Non-marker gene"
# Markers combined as boxplot
ggplot(df1, aes(x = Grouping, y = Dropseq)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Z-score (<- Bulk, Pooled ->)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nEnrichment or depletion of cell type markers in pooled versus bulk"
    , "\n"
    , "\nDrop-seq GZ CP"
    , "\nRatio of log2(CPM + 1) pooled versus bulk"))
ggsave(paste0(outGraph, "KnownMarkersCombined.pdf")
  , width = 10, height = 6)
# Markers individually bar graph
ggplot(df1, aes(x = hgnc_symbol, y = Dropseq)) +
  facet_wrap("Grouping", scales = "free_x") +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size = 12)) +
  ylab("Z-score (<- Bulk, Pooled ->)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nEnrichment or depletion of cell type markers in pooled versus bulk"
    , "\n"
    , "\nDrop-seq GZ CP"
    , "\nRatio of log2(CPM + 1) pooled versus bulk"))
ggsave(paste0(outGraph, "KnownMarkers.png"), width = 15, height = 12)
# Subset to marker groups of interest for paper
df1 <- df1[df1$Grouping %in% c(
  "Non-marker gene"
  , "GABAergic interneuron"
  , "RG"
  , "OPC"
  , "Microglia"
  , "oRG"
  , "Glutamatergic Neuron"
  , "Pericyte"
  , "Neuron"
  , "Endothelial Cell"
  , "IP"
  , "vRG")
  , ]
# Markers combined as boxplot
ggplot(df1, aes(x = Grouping, y = Dropseq)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Z-score (<- Bulk, Pooled ->)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nEnrichment or depletion of cell type markers in pooled versus bulk"
    , "\n"
    , "\nDrop-seq GZ CP"
    , "\nRatio of log2(CPM + 1) pooled versus bulk"))
ggsave(paste0(outGraph, "KnownMarkersCombined_paper.pdf")
  , width = 5, height = 6)
################################################################################

# ## Drop-seq union genes versus intersection genes
# # Intersection (only genes in both)
# mnBothDF <- merge(mean_bulk_ex_DF, mnPods_ex_DF, by.x = "row.names"
#   , by.y = "ensembl_gene_id")
# # Convert NAs to 0s
# mnBothDF[is.na(mnBothDF)] <- 0
# # Spearman correlation
# # Union
# sprUn <- round(cor(mnExDF$Mean_Bulk, mnExDF$MEAN_POOLED_DROPSEQ
#   , method = "spearman"), 2)
# sprUn
# # Intersection
# sprIn <- round(cor(mnBothDF$Mean_Bulk, mnBothDF$MEAN_POOLED
#   , method = "spearman"), 2)
# sprIn
# # Format for ggplot2
# df1 <- data.frame(Pooled = log(mnExDF$Mean_Bulk + 1, 2)
#   , Bulk = log(mnExDF$MEAN_POOLED_DROPSEQ + 1, 2))
# df1$DATASET <- "Union genes"
# df2 <- data.frame(Pooled = log(mnBothDF$Mean_Bulk + 1, 2)
#   , Bulk = log(mnBothDF$MEAN_POOLED + 1, 2))
# df2$DATASET <- "Intersection genes"
# ggDF <- rbind(df1, df2)
# # Add spearman to labels
# sprL <- c(sprUn, sprIn)
# labelsDatasetSpr <- c(paste0(unique(ggDF$DATASET), "\nSpearman: ", sprL))
# names(labelsDatasetSpr) <- unique(ggDF$DATASET)
# # ggplot
# ggplot(ggDF, aes(x = Pooled, y = Bulk, col = DATASET)) +
#   facet_wrap("DATASET", ncol = 2
#     , labeller = labeller(DATASET = labelsDatasetSpr)) +
#   geom_point(alpha = 0.5, shape = 1, size = 0.25) +
#   stat_smooth(col = "black") +
#   scale_color_manual(values = c("#e41a1c", "#377eb8")) +
#   theme(legend.position = "none") +
#   ylab("Bulk: log2(Mean CPM + 1)") +
#   xlab("Pooled: log2(Mean CPM + 1)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nPooled Drop-seq vs Bulk RNAseq"
#     , "\n"
#     , "\nIntersection (genes only detected in both) or union (genes in either)"
#     , "\nHuman Fetal Brain CP + GZ"
#     , "\nMean of CPM across samples and pools"))
# ggsave(paste0(outGraph, "ScatterIntersectionUnion.png"), width = 8, height = 5
#   , dpi = 600)
################################################################################



### Write table of subset gene lists

# List of dataframes, each dataframe is a subset
cols <- grep("SUBSET", names(ssl2MnCpmDF), value = TRUE)
names(cols) <- cols
l <- lapply(cols, function(dfcol) {
  l <- list(
    as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Bulk high", 1])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Bulk high", 2])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Pool high", 1])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Pool high", 2])
    )
  # Data frame of lists of different lengths
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  mat <- sapply(l, "[", i = seq.max)
  print(str(mat))
  colnames(mat) <- c("BULK_HIGH_ENSEMBL", "BULK_HIGH_HGNC", "POOL_HIGH_ENSEMBL", "POOL_HIGH_HGNC")
  # Convert dataframe to lists
  l <- as.list(as.data.frame(mat, stringsAsFactors = FALSE))
  return(l)
})
# Format column names
ll <- lapply(names(l), function(datasetName) {
  dataset <- l[[datasetName]]
  colNames <- paste0(datasetName, "_", names(dataset))
  names(dataset) <- colNames
  print(str(dataset))
  return(dataset)
})
# Append NAs to dataframes to make same length
n.obs <- sapply(ll, function(l) sapply(l, length))
seq.max <- seq_len(max(n.obs[1, ]))
lm <- lapply(ll, function(l) sapply(l, "[", i = seq.max))
m <- do.call("cbind", lm)
m[m == "0"] <- NA

# Write out
write.csv(m, file = paste0(outTable, "Subset_Gene_Lists.csv"))
################################################################################






# ## Intersect GO terms
#
# # Drop-seq GZ, Fluidigm HT GZ, Pollen
# # High pool
# go <- intersect(
#   intersect(gpPhiBloDsGzDF$term.id[
#       gpPhiBloDsGzDF$domain == "BP" |
#       gpPhiBloDsGzDF$domain == "MF" |
#       gpPhiBloDsGzDF$domain == "CC"]
#     , gpPhiBloFhGzDF$term.id)
#   , gpPhiBloPnDF$term.id)
# write.table(go, paste0(outTable, "GO_DsGz_FhGz_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(gpPloBhiDsGzDF$term.id[
#       gpPloBhiDsGzDF$domain == "BP" |
#       gpPloBhiDsGzDF$domain == "MF" |
#       gpPloBhiDsGzDF$domain == "CC"]
#     , gpPloBhiFhGzDF$term.id)
#   , gpPloBhiPnDF$term.id)
# write.table(go, paste0(outTable, "GO_DsGz_FhGz_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
#
# # Drop-seq, Fluidigm HT, Fluidigm LT
# # High pool
# go <- intersect(
#   intersect(gpPhiBloDsDF$term.id[
#       gpPhiBloDsDF$domain == "BP" |
#       gpPhiBloDsDF$domain == "MF" |
#       gpPhiBloDsDF$domain == "CC"]
#     , gpPhiBloFhDF$term.id)
#   , gpPhiBloFlDF$term.id)
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(gpPloBhiDsDF$term.id[
#       gpPloBhiDsDF$domain == "BP" |
#       gpPloBhiDsDF$domain == "MF" |
#       gpPloBhiDsDF$domain == "CC"]
#     , gpPloBhiFhDF$term.id)
#   , gpPloBhiFlDF$term.id)
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
#
# # Drop-seq, Fluidigm HT, Fluidigm LT, Pollen
# # High pool
# go <- intersect(
#   intersect(
#     intersect(gpPhiBloDsDF$term.id[
#         gpPhiBloDsDF$domain == "BP" |
#         gpPhiBloDsDF$domain == "MF" |
#         gpPhiBloDsDF$domain == "CC"]
#       , gpPhiBloFhDF$term.id)
#     , gpPhiBloFlDF$term.id)
#   , gpPhiBloPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(
#     intersect(gpPloBhiDsDF$term.id[
#         gpPloBhiDsDF$domain == "BP" |
#         gpPloBhiDsDF$domain == "MF" |
#         gpPloBhiDsDF$domain == "CC"]
#       , gpPloBhiFhDF$term.id)
#     , gpPloBhiFlDF$term.id)
#   , gpPloBhiPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
#
# # Drop-seq, not in Fluidigm HT, Fluidigm LT
# # High pool
# go <- setdiff(
#   setdiff(gpPhiBloDsDF$term.id[
#       gpPhiBloDsDF$domain == "BP" |
#       gpPhiBloDsDF$domain == "MF" |
#       gpPhiBloDsDF$domain == "CC"]
#     , gpPhiBloFhDF$term.id)
#   , gpPhiBloPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_NotInFhFl_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- setdiff(
#   setdiff(gpPloBhiDsDF$term.id[
#       gpPloBhiDsDF$domain == "BP" |
#       gpPloBhiDsDF$domain == "MF" |
#       gpPloBhiDsDF$domain == "CC"]
#     , gpPloBhiFhDF$term.id)
#   , gpPloBhiFlDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_NotInFhFl_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)

## Gene length for each GO ID

dir.create("../analysis/graphs/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_GO_CDS_Length", recursive = TRUE)
# Loop through gprofiler results for each subset
lapply(names(gproLDF), function(name) {
  print(name)
  gproDF <- gproLDF[[name]]
  if (nrow(gproDF) > 0) {
    # Only BP GO terms
    gproDF <- gproDF[gproDF$domain == "BP", ]
    # Genes
    l <- strsplit(gproDF$intersection, ",")
    names(l) <- gproDF$term.name
    # Data frame of genes and term names
    ggDF <- do.call("rbind", lapply(names(l), function(name) {
      data.frame(hgnc_symbol = l[[name]], term.name = name)
    }))
    # Add CDS length
    ggDF <- merge(ggDF, mnExDF[c("hgnc_symbol", "cds_length")], by = "hgnc_symbol")
    # Take longest CDS (some hgnc_symbol occur twice in from matching multiple
    # ensembl IDs)
    ggDF <- aggregate(cds_length ~ hgnc_symbol + term.name, data = ggDF, max)
    # Add all genes
    df <- mnExDF[mnExDF$cds_length > 0, c("hgnc_symbol", "cds_length")]
    df[ ,3] <- df$cds_length
    df[ ,2] <- "All genes"
    names(df)[2:3] <- c("term.name", "cds_length")
    ggDF <- rbind(ggDF, df)
    # Term to color all genes box only
    ggDF$AREA_COLOR <- "Genes GO intersection"
    ggDF$AREA_COLOR[ggDF$term.name == "All genes"] <- "All genes"
    # ggplot
    ggplot(ggDF, aes(y = cds_length, x = term.name, fill = AREA_COLOR)) +
      geom_boxplot() +
      scale_fill_brewer(type = "qual", palette = 6) +
      coord_flip(ylim = c(0, 6000)) +
      theme(text = element_text(size = 12))
    ggsave(paste0(
      "../analysis/graphs/Pooled_Vs_Bulk_GO_CDS_Length/Pooled_Vs_Bulk_GO_CDS_Length_"
      , name, ".pdf"), height = (length(unique(ggDF$term.name))*.1 + 4)
      , width = 12)
  }
  else {print(paste0("No significant GO terms for ", name))}
})

# # "BP"  "CC"  "MF"  "cor" "hp"  "hpa" "keg" "mi"  "omi" "rea" "tf"
#
# pdf(paste0(outGraph, "gProfileR.pdf"), width = 9, height = 7.5)
# # Title page
# plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
# text(5, 8, "GO term enrichment using gProfileR of genes")
# text(5, 7, "higher in pooled versus bulk or higher in bulk versus pooled")
# text(5, 6, "Multiple hypothesis correction using g:SCS")
# text(5, 5, "g:SCS stringency is in between FDR and Bonferroni and was designed for GO term enrichment analysis")
# text(5, 4, "All human genes used as background")
# # Drop-seq
# Graph_GO_Enrichment(gpPloBhiDsDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Drop-seq"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloDsDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Drop-seq"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Fluidigm HT
# Graph_GO_Enrichment(gpPloBhiFhDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm HT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFhDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidigm HT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Fluidigm LT
# Graph_GO_Enrichment(gpPloBhiFlDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm LT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFlDF, "BP", 3
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidgm LT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Drop-seq GZ
# Graph_GO_Enrichment(gpPloBhiDsGzDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Drop-seq"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloDsGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Drop-seq"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# # Fluidigm HT GZ
# Graph_GO_Enrichment(gpPloBhiFhGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm HT"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFhGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidigm HT"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# # Pollen
# Graph_GO_Enrichment(gpPloBhiPnDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Pollen"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloPnDF, "BP", 2
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Pollen"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# dev.off()

# # GO Terms molecular function
# ggDF <- subset(gpPlBhDF, domain == "MF" & relative.depth == 5)
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
#
# # Kegg
# ggDF <- subset(gpPlBhDF, domain == "keg")
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
#
# # Omi
# ggDF <- subset(gpPlBhDF, domain == "omi")
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
################################################################################

### Expression of genes biased towards pool or bulk by cluster

Heatmap_Expression_Clusters_Format <- function (geneListDF, expLimitHigh
  , expLimitLow){
  # Subset to genes, merge to keep duplicate genes
  # e.g. if plotting cluster markers, may be the same marker of multiple clusters
  ggDF <- merge(geneListDF, centSO@scale.data, by.x = 1, by.y = "row.names"
    , all.x = TRUE)
  row.names(ggDF) <- paste0(length(ggDF[ ,1]):1, "_", ggDF[ ,1])
  ggDF <- ggDF[ ,-c(1)]
  # Order by clustering
  idx <- match(colnames(ggDF), names(sort(centSO@ident)))
  ggDF <- ggDF[ ,idx]
  # Format for ggplot2
  ggDF <- as.matrix(ggDF)
  ggDF <- melt(ggDF)
  # Add clusters
  idx <- match(ggDF$Var2, names(centSO@ident))
  ggDF$CLUSTERS <- centSO@ident[idx]
  # Set sample order by clustering
  ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
  # Set expression limits
  ggDF$value[ggDF$value > expLimitHigh] <- expLimitHigh
  ggDF$value[ggDF$value < expLimitLow] <- expLimitLow
  # Return data frame formatted for ggplot
  return(ggDF)
}

Heatmap_Expression_Clusters <- function (ggDF, graphTitle) {
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(~CLUSTERS, space = "free", scales = "free") +
    # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
    scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", graphTitle
      , "\nMean centered and scaled normalized expression"
      , "\n"))
}

geneListDF <- data.frame(GENE = ssl2MnCpmDF$HGNC_SYMBOL[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high"])
geneListDF <- data.frame(GENE = geneListDF[! geneListDF$GENE == 0, ])
ggDF <- Heatmap_Expression_Clusters_Format(geneListDF, 5, -1)
gg <- Heatmap_Expression_Clusters(
  ggDF, "Genes higher in bulk versus pooled Drop-seq")
ggsave(paste0(outGraph, "Heatmap_Dropseq_BulkHigh.png")
  , width = 10, height = 10, dpi = 600)

geneListDF <- data.frame(GENE = ssl2MnCpmDF$HGNC_SYMBOL[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high"])
geneListDF <- data.frame(GENE = geneListDF[! geneListDF$GENE == 0, ])
ggDF <- Heatmap_Expression_Clusters_Format(geneListDF, 5, -1)
gg <- Heatmap_Expression_Clusters(
  ggDF, "Genes higher in pooled Drop-seq versus bulk")
ggsave(paste0(outGraph, "Heatmap_Dropseq_PoolHigh.png")
  , width = 10, height = 10, dpi = 600)
################################################################################

### Known cell type markers from Luis

## Intersection with pooled vs bulk subsets
Intersect_KnownMarks_Subsets <- function(sub, bias) {
  genes <- ssl2MnCpmDF[ssl2MnCpmDF[ ,c(sub)] == bias, ]$HGNC_SYMBOL
  kmDF$Gene.Symbol %in% genes
}
l <- list(Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Pool high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Pool high")
)
nr <- max(sapply(l, length))
nc <- length(l)
m <- matrix(NA, nr, nc)
for(i in 1:length(l)){
  m[ ,i] <- l[[i]]
}
colnames(m) <- c(
  "SUBSET_DROPSEQ_BULK_HIGH"
  , "SUBSET_DROPSEQ_POOL_HIGH"
  , "SUBSET_DROPSEQ_0_BULK_HIGH"
  , "SUBSET_DROPSEQ_0_POOL_HIGH"
  , "SINGLE_CELL_BIAS_GZ_BULK_HIGH"
  , "SINGLE_CELL_BIAS_GZ_POOL_HIGH"
  , "SINGLE_CELL_BIAS_BULK_HIGH"
  , "SINGLE_CELL_BIAS_POOL_HIGH"
)
rnames <- paste0(kmDF$Group, " ", kmDF$Gene.Symbol)
row.names(m) <- rnames
ggDF <- melt(m)
ggDF$Group <- gsub(" .*", "", ggDF$Var1)
ggDF$Var1 <- gsub(".* ", "", ggDF$Var1)
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_wrap("Group", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  ylab("Variant location") +
  xlab("Cell ID") +
  ggtitle(paste0(graphCodeTitle
    , "\nVariants called for each Cell: Fluidigm HT, Human VZ and CP"
    , "\n"))
ggsave(paste0(outGraph, "KnownMarkers_In_Subsets_Heatmap.pdf"), height = 24, width = 16)


## Intersection with pooled vs bulk subsets
Intersect_KnownMarks_Subsets <- function(sub, bias) {
  genes <- ssl2MnCpmDF[ssl2MnCpmDF[ ,c(sub)] == bias, ]$HGNC_SYMBOL
  as.list(kmDF[kmDF$Gene.Symbol %in% genes, 1:3])
}
ll <- list(Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Pool high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Pool high")
  )
l <- unlist(ll, recursive = FALSE)
nr <- max(sapply(ll, length))
nc <- 3*length(ll)
m <- matrix(NA, nr, nc)
for(i in 1:length(ll)){
  m[ ,i] <- as.character(l[[i]][1:nr])
}
m[is.na(m)] <- ""
colnames(m) <- c(
  rep("SUBSET_DROPSEQ_BULK_HIGH", 3)
  , rep("SUBSET_DROPSEQ_POOL_HIGH", 3)
  , rep("SUBSET_DROPSEQ_0_BULK_HIGH", 3)
  , rep("SUBSET_DROPSEQ_0_POOL_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_GZ_BULK_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_GZ_POOL_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_BULK_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_POOL_HIGH", 3)
)
write.csv(m, paste0(outTable, "KnownMarkers_In_Subsets.csv"), quote = FALSE
  , row.names = FALSE)

################################################################################

# ### Cluster markers from Seurat
#
# # Z-score
# sdDs <- sd(l2MnCpmDF$RATIO_DROPSEQ)
# mnDs <- mean(l2MnCpmDF$RATIO_DROPSEQ)
# # Subset to marker genes of interest for Luis' excel file
# df <- merge(l2MnCpmDF, mkDF, by.x = "HGNC_SYMBOL", by.y = "gene")
# # df <- l2MnCpmDF[l2MnCpmDF$HGNC_SYMBOL %in% kmDF$Gene.Symbol, ]
# df$Z_SCORE <- (df$RATIO_DROPSEQ - mnDs) / sdDs
# df$cluster <- as.factor(df$cluster)
#
# # Markers individually bar graph
# ggplot(df, aes(x = HGNC_SYMBOL, y = Z_SCORE)) +
#   facet_wrap("Grouping", scales = "free_x") +
#   geom_bar(stat = "identity") +
#   coord_cartesian(ylim = c(-5, 5)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme(text = element_text(size = 12)) +
#   ylab("Z-score (<- Bulk, Pooled ->)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nEnrichment or depletion of cluster markers in pooled versus bulk"
#     , "\n"
#     , "\nDrop-seq GZ CP"
#     , "\nRatio of log2(CPM + 1) pooled versus bulk"))
# ggsave(paste0(outGraph, "SeuratMarkers.pdf"), width = 12, height = 12)
#
# # Markers combined as boxplot
# ggplot(df, aes(x = cluster, y = Z_SCORE)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ylab("Z-score (<- Bulk, Pooled ->)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nEnrichment or depletion of cluster markers in pooled versus bulk"
#     , "\n"
#     , "\nDrop-seq GZ CP"
#     , "\nRatio of log2(CPM + 1) pooled versus bulk"))
# ggsave(paste0(outGraph, "SeuratMarkersCombined.pdf"), width = 8, height = 6)
# ################################################################################
