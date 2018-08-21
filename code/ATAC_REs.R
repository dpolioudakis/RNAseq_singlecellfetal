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
require(ggplot2)
require(cowplot)
require(ggbeeswarm)
source("Function_Library.R")

options(stringsAsFactors = FALSE)

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

# Cluster DE genes
cluster_DE_DF <- read.table(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Cell type specific genes
celltype_specific_DF <- read.csv("../analysis/tables/Seurat_ClassMarkers/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClassMarkers_Markers_Class.csv")

# ATAC regulatory elements
atac_HiC_RE_DF <- read.csv("../source/Gene_Lists/TorreUbieta_2018_TS2.csv"
  , header = TRUE)
atac_noHiC_RE_DF <- read.csv("../source/All_ATACcor_toTSS.csv", header = TRUE)

# Nowakowski cluster DE table
nowakowski_DE_DF <- read.csv(
  "../nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)

# Lake 2017 cluster DE tables
lake_DE_DF <- read.csv(
  "../source/Lake_2018_TS3_Cluster_DE_Genes.csv"
  , skip = 4, header = TRUE, fill = TRUE
)

# biomaRt gene info
bmDF <- read.csv(
  "../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "ATAC_REs.R"
outGraph <- "../analysis/graphs/ATAC_REs/ATAC_REs_"
outTable <- "../analysis/tables/ATAC_REs/ATAC_REs_"
outData <- "../analysis/analyzed_data/ATAC_REs/ATAC_REs_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
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
################################################################################

### Cell type enriched REs

## Drop-seq

# Remove cluster 16 (low quality cells)
cluster_DE_DF <- cluster_DE_DF[cluster_DE_DF$Cluster != 16, ]

# Merge ATAC+HiC enhancers table and cell type enriched genes table by ensid
atacHiC_RE_clusterDE_DF <- merge(
  atac_HiC_RE_DF, cluster_DE_DF
  , by.x = "ENSGID", by.y = "Ensembl"
)
# Merge ATAC without HiC enhancers table and cell type enriched genes table by ensid
ataconly_RE_clusterDE_DF <- merge(
  atac_noHiC_RE_DF, cluster_DE_DF
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atacHiC_RE_clusterDE_DF
  , file = paste0(outTable, "ATACandHiC_CellTypeEnrichment.csv")
  , quote = FALSE
)
write.csv(ataconly_RE_clusterDE_DF
  , file = paste0(outTable, "ATAConly_CellTypeEnrichment.csv")
  , quote = FALSE
)
# Format for paper
paper_atac_DF <- ataconly_RE_clusterDE_DF
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
  , "8" = "Cycling progenitor S phase"
  , "10" = "Cycling progenitor G2/M phase"
  , "2" = "IPC"
  , "0" = "Excitatory neuron new born migrating"
  , "1" = "Excitatory neuron"
  , "4" = "Excitatory neuron (collosal)"
  , "3" = "Deep layer excitatory neuron 1"
  , "13" = "Deep layer excitatory neuron 2"
  , "5" = "Interneuron (SST)"
  , "6" = "Interneuron (CALB2)"
  , "11" = "Oligodendrocyte precursor"
  , "12" = "Endothelial"
  , "14" = "Pericyte"
  , "15" = "Microglia"
)
idx <- match(paper_atac_DF$Cluster_number, names(cluster_annot))
paper_atac_DF$Cluster <- cluster_annot[idx]
paper_atac_DF$Cluster_number <- factor(paper_atac_DF$Cluster_number
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
paper_atac_DF <- paper_atac_DF[order(paper_atac_DF$Cluster_number), ]
write.csv(paper_atac_DF
  , file = paste0(outTable, "ATAConly_CellTypeEnrichment_paper.csv")
  , quote = FALSE, row.names = FALSE
)

# Merge ATAC+HiC enhancers table and major cell type specific genes table by ensid
atacRE_celltypeDE_DF <- merge(
  atac_HiC_RE_DF, celltype_specific_DF
  , by.x = "ENSGID", by.y = "Ensembl"
)
# Merge ATAC without HiC enhancers table and major cell type specific genes table by ensid
atacRE_noHiC_celltypeDE_DF <- merge(
  atac_noHiC_RE_DF, celltype_specific_DF
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atacRE_celltypeDE_DF
  , file = paste0(outTable, "ATACandHiC_CellTypeSpecific.csv")
  , quote = FALSE
)
write.csv(atacRE_noHiC_celltypeDE_DF
  , file = paste0(outTable, "ATAConly_CellTypeSpecific.csv")
  , quote = FALSE
)

## Nowakowski
# Subset to higher enrichment
subset_DE_DF <- nowakowski_DE_DF[
  nowakowski_DE_DF$avg_diff > 0.75, c("gene", "avg_diff", "cluster")]
colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
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
subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
subset_DE_DF <- subset_DE_DF[subset_DE_DF$Class != "NA", ]
# Add ensembl ID
idx <- match(subset_DE_DF$Gene, bmDF$hgnc_symbol)
subset_DE_DF$Ensembl <- bmDF$ensembl_gene_id[idx]
# Merge with ATAC REs
atacRE_nowakowskiDE_DF <- merge(
  atac_noHiC_RE_DF, subset_DE_DF
  , by.x = "ENSGID", by.y = "Ensembl"
)
write.csv(atacRE_nowakowskiDE_DF
  , file = paste0(
    outTable, "ATAConly_Nowakowski_CellTypeEnrichment.csv")
  , quote = FALSE
)

# ## Lake
# # Subset to higher enrichment
# subset_DE_DF <- lake_DE_DF[
#   lake_DE_DF$"Average.Difference..log.fold.change." > 0.25
#   , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]
# colnames(subset_DE_DF) <- c("Gene", "Enrichment", "Cluster")
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
# subset_DE_DF <- Classify_Cells_By_Type(subset_DE_DF, class_cluster_idx)
# subset_DE_DF <- subset_DE_DF[subset_DE_DF$Class != "NA", ]
# # Add ensembl ID
# idx <- match(subset_DE_DF$Gene, bmDF$hgnc_symbol)
# subset_DE_DF$Ensembl <- bmDF$ensembl_gene_id[idx]
# # Merge with ATAC REs
# atacRE_lakeDE_DF <- merge(
#   atac_HiC_RE_DF, subset_DE_DF
#   , by.x = "ENSGID", by.y = "Ensembl"
# )
# write.csv(atacRE_lakeDE_DF
#   , file = paste0(
#     outTable, "ATAConly_Lake_CellTypeEnrichment.csv")
#   , quote = FALSE
# )
################################################################################

### Plots ATAC + HiC

# Add biomart gene info
idx <- match(atac_HiC_RE_DF$ENSGID, bmDF$ensembl_gene_id)
atac_HiC_RE_DF <- cbind(atac_HiC_RE_DF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
)
# Add biomart gene info
idx <- match(atacHiC_RE_clusterDE_DF$ENSGID, bmDF$ensembl_gene_id)
atacHiC_RE_clusterDE_DF <- cbind(atacHiC_RE_clusterDE_DF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
)

# Reorder clusters
atacHiC_RE_clusterDE_DF$Cluster <- factor(atacHiC_RE_clusterDE_DF$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Remove cluster 16
atacHiC_RE_clusterDE_DF <- atacHiC_RE_clusterDE_DF[
  ! is.na(atacHiC_RE_clusterDE_DF$Cluster), ]
# Distance enhancer end to promoter start
atacHiC_RE_clusterDE_DF$Distance_Enhancer_Promoter <- abs(
  atacHiC_RE_clusterDE_DF$enhancerend -
  atacHiC_RE_clusterDE_DF$promoterpeakstart
)
atacHiC_RE_clusterDE_DF$Enhancer_Size <- abs(
  atacHiC_RE_clusterDE_DF$enhancerend -
  atacHiC_RE_clusterDE_DF$enhancerstart
)

# Plots
# Distance_Enhancer_Promoter
ggplot(atacHiC_RE_clusterDE_DF, aes(
  x = Cluster, y = Distance_Enhancer_Promoter, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  geom_quasirandom(size = 0.1) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAChiC_Distance_Enhancer_Promoter_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atacHiC_RE_clusterDE_DF, aes(
  x = gene_biotype, y = Distance_Enhancer_Promoter, fill = gene_biotype)) +
  geom_boxplot(outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  geom_quasirandom(size = 0.1) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAChiC_Distance_Enhancer_Promoter_By_Biotype_boxplot.pdf")
  , width = 6, height = 5)
ggplot(atacHiC_RE_clusterDE_DF, aes(x = Distance_Enhancer_Promoter)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAChiC_Distance_Enhancer_Promoter_density.pdf"))
# Enhancer_Size
ggplot(atacHiC_RE_clusterDE_DF, aes(
  x = Cluster, y = Enhancer_Size, fill = Cluster)) +
  # facet_wrap(~gene_biotype) +
  geom_quasirandom(size = 0.1) +
  ggplot_set_theme_publication +
  geom_boxplot(outlier.shape = NA) +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(outGraph, "ATAChiC_Enhancer_Size_boxplot.pdf")
  , width = 7, height = 5)
ggplot(atacHiC_RE_clusterDE_DF, aes(x = Enhancer_Size)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(outGraph, "ATAChiC_Enhancer_Size_density.pdf")
  , width = 7, height = 5)


# Enhancers per gene
ggDF <- melt(with(atacHiC_RE_clusterDE_DF, table(ENSGID, Cluster)))
# Add biomart gene info
idx <- match(ggDF$ENSGID, bmDF$ensembl_gene_id)
ggDF <- cbind(ggDF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
)
ggDF <- ggDF[ggDF$value != 0, ]
ggDF$Cluster <- factor(ggDF$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Enhancers per gene boxplot
ggplot(ggDF, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = 0.1) +
  # coord_cartesian(ylim = c(0, 10)) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ylab("Enhancers per gene") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(outGraph, "ATAChiC_Enhancers_Per_Gene_boxplot.pdf")
  , width = 7, height = 5)
# Enhancers per gene density plot
ggplot(ggDF, aes(x = value)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(outGraph, "ATAChiC_Enhancer_Per_Gene_density.pdf")
  , width = 7, height = 5)

# Enhancers per gene versus CDS length
sapply(split(ggDF, ggDF$Cluster), function(subset_ggDF){
  cor(subset_ggDF$cds_length, subset_ggDF$value, use = "complete")
})
ggplot(ggDF, aes(x = value, y = cds_length, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("CDS length") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus CDS length"))
ggsave(paste0(outGraph, "ATAChiC_Enhancers_Per_Gene_vs_CDS_Length.pdf")
  , width = 7, height = 5)

# Enhancers per gene versus GC content
ggplot(ggDF, aes(x = value, y = percentage_gc_content, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("GC content") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC+HiC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus gene GC content"))
ggsave(paste0(outGraph, "ATAChiC_Enhancers_Per_Gene_vs_GC_Content.pdf")
  , width = 7, height = 5)
################################################################################

### Plots ATAC only

# Add biomart gene info
idx <- match(atac_noHiC_RE_DF$ENSGID, bmDF$ensembl_gene_id)
atac_noHiC_RE_DF <- cbind(atac_noHiC_RE_DF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
)
# Add biomart gene info
idx <- match(ataconly_RE_clusterDE_DF$ENSGID, bmDF$ensembl_gene_id)
ataconly_RE_clusterDE_DF <- cbind(ataconly_RE_clusterDE_DF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
)

# Reorder clusters
ataconly_RE_clusterDE_DF$Cluster <- factor(ataconly_RE_clusterDE_DF$Cluster
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
# Remove cluster 16
ataconly_RE_clusterDE_DF <- ataconly_RE_clusterDE_DF[
  ! is.na(ataconly_RE_clusterDE_DF$Cluster), ]
# Distance enhancer end to promoter start
ataconly_RE_clusterDE_DF$Distance_Enhancer_Promoter <- abs(
  ataconly_RE_clusterDE_DF$enhancerend -
  ataconly_RE_clusterDE_DF$promoterstart
)
ataconly_RE_clusterDE_DF$Enhancer_Size <- abs(
  ataconly_RE_clusterDE_DF$enhancerend -
  ataconly_RE_clusterDE_DF$enhancerstart
)

# Remove enhancer entries that are actually promoters
ss_ataconly_RE_clusterDE_DF <- ataconly_RE_clusterDE_DF[
  ataconly_RE_clusterDE_DF$enhancerstart != ataconly_RE_clusterDE_DF$promoterstart
  , ]

# Plots
# Distance_Enhancer_Promoter
ggplot(ataconly_RE_clusterDE_DF, aes(
  x = Cluster, y = Distance_Enhancer_Promoter, fill = Cluster)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAConly_Distance_Enhancer_Promoter_boxplot.pdf")
  , width = 7, height = 5)
ggplot(ataconly_RE_clusterDE_DF, aes(
  x = gene_biotype, y = Distance_Enhancer_Promoter, fill = gene_biotype)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAConly_Distance_Enhancer_Promoter_By_Biotype_boxplot.pdf")
  , width = 6, height = 5)
ggplot(ataconly_RE_clusterDE_DF, aes(x = Distance_Enhancer_Promoter)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nDistance of enhancer end to promoter start"))
ggsave(paste0(outGraph, "ATAConly_Distance_Enhancer_Promoter_density.pdf"))
# Enhancer_Size
ggplot(ataconly_RE_clusterDE_DF, aes(
  x = Cluster, y = Enhancer_Size, fill = Cluster)) +
  # facet_wrap(~gene_biotype) +
  coord_cartesian(ylim = c(0, 5000)) +
  geom_violin(fill = "lightgrey", color = "lightgrey") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(outGraph, "ATAConly_Enhancer_Size_boxplot.pdf")
  , width = 7, height = 5)
ggplot(ataconly_RE_clusterDE_DF, aes(x = Enhancer_Size)) +
  geom_density() +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancer size"))
ggsave(paste0(outGraph, "ATAConly_Enhancer_Size_density.pdf")
  , width = 7, height = 5)


# Enhancers per gene
ggDF <- melt(with(ataconly_RE_clusterDE_DF, table(ENSGID, Cluster)))
# Add biomart gene info
idx <- match(ggDF$ENSGID, bmDF$ensembl_gene_id)
ggDF <- cbind(ggDF
  , bmDF[idx, c("percentage_gc_content", "cds_length")]
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
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(outGraph, "ATAConly_Enhancers_Per_Gene_boxplot.pdf")
  , width = 7, height = 5)
# Enhancers per gene density plot
ggplot(ggDF, aes(x = value)) +
  geom_histogram(binwidth = 1) +
  ggplot_set_theme_publication +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene"))
ggsave(paste0(outGraph, "ATAConly_Enhancer_Per_Gene_histogram.pdf")
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
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus CDS length"))
ggsave(paste0(outGraph, "ATAConly_Enhancers_Per_Gene_vs_CDS_Length.pdf")
  , width = 7, height = 6)

# Enhancers per gene versus GC content
ggplot(ggDF, aes(x = value, y = percentage_gc_content, color = Cluster)) +
  facet_wrap(~Cluster) +
  geom_point(size = 0.1) +
  theme(legend.position = "none") +
  ylab("GC content") +
  xlab("Enchancers per gene") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nATAC enhancers with cell type specific enrichment"
    , "\nEnhancers per gene versus gene GC content"))
ggsave(paste0(outGraph, "ATAConly_Enhancers_Per_Gene_vs_GC_Content.pdf")
  , width = 7, height = 6)
################################################################################
