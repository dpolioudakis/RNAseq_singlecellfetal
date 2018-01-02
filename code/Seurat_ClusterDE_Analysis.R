# Damon Polioudakis
# 2017-08-28
# DE genes for each cluster

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
# require(irlba)
require(gridExtra)
require(fdrtool)
require(cowplot)
require("ggdendro")
source("Function_Library.R")

## Inputs

# Seurat
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Biomart to add ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Kang
# Non regressed data
load("../neurogenesis/orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
kangExM <- t(datExpr)
rm(datExpr)
# sampleKey is ProcessedKangMetaData.csv
kangMtDF <- sampleKey
kangAnnotRaw <- read.csv("../neurogenesis/orig.data/InVivoData/annot.csv", row.names = 1)
kangAnnotnet2 <- kangAnnotRaw[which(rownames(kangAnnotRaw) %in% rownames(kangExM)), ]

# Miller
load("../neurogenesis/orig.data/LCMDE/AllenLCM.Rdata")
millerExDF = AllenLCM$datExpr
millerMtDF = AllenLCM$datTraits
millerZonesDF = read.csv("../neurogenesis/orig.data/LCMDE/LCM_Zones_CPio.csv")
millerAnnotRawDF = read.csv("../neurogenesis/orig.data/LCMDE/annot.csv", row.names = 1)

# DE of clusters
inTables <- list.files("../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/"
  , full.names = TRUE)
inTables <- inTables[grep("Cluster\\d", inTables, perl = TRUE)]

## Variables
graphCodeTitle <- "Seurat_ClusterDE_Analysis.R"
outGraph <- "../analysis/graphs/Seurat_ClusterDE_Analysis_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"
outTable <- "../analysis/tables/Seurat_ClusterDE_Analysis_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"
outData <- "../analysis/Seurat_ClusterDE_Analysis_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outData)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 14)))
theme_update(plot.title = element_text(size = 14))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

################################################################################

### Format and compile

## DE tables into one table

# Loop through DE text files and compile into one table
ldf <- lapply(inTables, function(inDE) {
  df <- read.table(inDE, header = TRUE)
  df <- df[order(-df$LOG_FC), ]
})
clusterDeDF <- do.call("rbind", ldf)
head(clusterDeDF)
# Add ensembl IDs
clusterDeDF$ENSEMBL <- bmDF$ensembl_gene_id[
  match(clusterDeDF$GENE, bmDF$hgnc_symbol)]
# Fill in IDs from GENE column that were kept as ensembl ID b/c no gene sym
clusterDeDF$ENSEMBL[is.na(clusterDeDF$ENSEMBL)] <-
  clusterDeDF$GENE[is.na(clusterDeDF$ENSEMBL)]
# Check
length(grep("ENSG", clusterDeDF$ENSEMBL))
nrow(clusterDeDF)
# # Filter FDR < 0.05
# deDF <- deDF[deDF$FDR < 0.05, ]

# Write out as tab delimited
write.table(x = clusterDeDF
  , file = paste0(outTable, "ClusterX_Vs_All_Clusters.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)


## Kang convert probe ids to ensembl
etzId <- kangAnnotRaw$ENTREZ_ID[match(row.names(kangExM), rownames(kangAnnotRaw))]
ensId <- bmDF$ensembl_gene_id[match(etzId, bmDF$entrezgene)]
row.names(kangExM) <- ensId


## Miller

# Subset to samples from neocortix layers
toMatch <- paste(millerZonesDF$Label, collapse = "|")
wellIDs <- millerMtDF$well_id[grep(toMatch, millerMtDF$structure_acronym)]
millerMtDF <- millerMtDF[millerMtDF$well_id %in% wellIDs, ]
millerExDF <- millerExDF[ ,colnames(millerExDF) %in% wellIDs]

# Subset to entrez annotated probes
millerExDF <- millerExDF[! is.na(millerAnnotRawDF$ENTREZ_ID), ]
millerAnnot <- millerAnnotRawDF[! is.na(millerAnnotRawDF$ENTREZ_ID), ]

# Convert Miller expression matrix to entrez
# Miller rows match Miller annotation rows
row.names(millerExDF) <- millerAnnot$ENTREZ_ID

# Get maximum expression probe
genes = unique(millerAnnot$ENTREZ_ID)
keepind = matrix(nrow = 0, ncol = 0);
for (ii in 1:length(genes)) {
  genematchind = which(millerAnnot$ENTREZ_ID == genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(millerExDF[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}
millerExDF = millerExDF[keepind, ];
rownames(millerExDF) = genes;
millerAnnot = millerAnnot[keepind, ];

# Convert entrez to ensembl
# Subset to genes that have ensembl
millerExDF <- millerExDF[row.names(millerExDF) %in% bmDF$entrezgene, ]
millerExDF$Ensembl <- bmDF$ensembl_gene_id[
  match(row.names(millerExDF), bmDF$entrezgene)]
# Take average expression of genes that have same ensembl
# millerExDF <- aggregate(.~Ensembl, millerExDF, mean)
##### Averaging not working correctly, removed duplicates as place holder
millerExDF <- millerExDF[! duplicated(millerExDF$Ensembl), ]
row.names(millerExDF) <- millerExDF$Ensembl
millerExDF <- millerExDF[ ,-1]


# row.names(millerExDF) <- millerAnnot$SYMBOL


# Add Zone column combining laminar zones from different regions
millerMtDF$Zone <- millerMtDF$structure_acronym
millerMtDF$Zone <- gsub("[f|t|p|o]SZ.*o",   "SZo", millerMtDF$Zone)
millerMtDF$Zone <- gsub("[f|t|p|o]SZ.*i",   "SZi", millerMtDF$Zone)
for (i in 1:nrow(millerZonesDF)) {
  millerMtDF$Zone <- gsub(
    as.character(millerZonesDF[i,"Label"])
    , as.character(millerZonesDF[i,"Zone"])
    , millerMtDF$Zone)
}
millerMtDF$Zone <- gsub("CPori", "CPi", millerMtDF$Zone)
millerMtDF$Zone <- gsub("SZori", "SZi", millerMtDF$Zone)

# # Subset to frontal cortex
# millerMtDF <- millerMtDF[grep("^f.*", millerMtDF$structure_acronym), ]
# millerExDF <- millerExDF[ ,grep("^f.*", millerMtDF$structure_acronym)]
################################################################################

### DE Heatmaps, top 20 DE table, and save seurat object

# Top 20 DE genes
clusterDe20DF <- data.frame(
  clusterDeDF %>% group_by(CLUSTER) %>% top_n(20, LOG_FC))
# Write out as tab delimited
write.table(x = clusterDe20DF
  , file = paste0(outTable, "ClusterX_Vs_All_Clusters_Top20.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)

# Save seurat object with cluster DE data frame
# save(centSO, noCentExM, metDF, clusterDeDF, file = paste0(outData, "seuratO.Robj"))

# Violin plots of fold changes by cluster
ggDF <- clusterDeDF
ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
ggplot(ggDF, aes(x = CLUSTER, y = LOG_FC)) +
  geom_violin(aes(fill = CLUSTER)) +
  geom_jitter(size = 0.05) +
  ylab("Log fold change") +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nFold changes for DE genes by cluster"
    , "\n"))
ggsave(paste0(outGraph, "DE_Violin.png"), width = 12, height = 8)

# Histograms of fold changes by cluster
ggDF <- clusterDeDF
ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
ggplot(ggDF, aes(x = LOG_FC)) +
  geom_histogram() +
  facet_wrap(~CLUSTER) +
  ylab("Counts") +
  xlab("Log fold change") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nHistogram of fold changes for DE genes by cluster"
    , "\n"))
ggsave(paste0(outGraph, "DE_Histogram.png"), width = 12, height = 8)

# MA plots for each cluster
ggL <- lapply(unique(clusterDeDF$CLUSTER), function(clusterID) {
  print(clusterID)
  # Subset expression matrix to cluster
  cdf <- noCentExM[ ,centSO@ident == clusterID, drop = FALSE]
  # Subset expression matrix to all other cells
  ndf <- noCentExM[ , ! centSO@ident == clusterID]
  # Fold change
  lfc <- rowMeans(cdf) - rowMeans(ndf)
  # Mean expression
  mn <- rowMeans(noCentExM)
  # Combine in DF for ggplot
  ggDF <- data.frame(LOG_FOLD_CHANGE = lfc, MEAN_EXPRESSION = mn
    , GENE = names(lfc))
  ggDF$DE_GENE <- row.names(ggDF) %in% clusterDeDF$GENE[clusterDeDF$CLUSTER == clusterID]
  # ggplot
  gg <- ggplot(ggDF, aes(x = MEAN_EXPRESSION, y = LOG_FOLD_CHANGE)) +
    geom_point(size = 0.1, alpha = 0.25, aes(color = DE_GENE)) +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "none") +
    xlab("Mean normalized expression") +
    ylab("Log fold change") +
    ggtitle(clusterID)
  return(gg)
})
# plot_grid combine cluster heatmaps
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMA plots: Cells in cluster Vs all other cells"
  , "\nColor indicates > 0.2 log fold change and FDR > 0.05"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "MAplot.png"), width = 12, height = 4+length(ggL))

# Number of DE genes per cluster barplots
ggDF1 <- data.frame(table(clusterDeDF$CLUSTER))
ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4]))
ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7]))
p1 <- ggplot(ggDF1, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.2 log fold change")
p2 <- ggplot(ggDF2, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.4 log fold change")
p3 <- ggplot(ggDF3, aes(x = Var1, y = Freq)) +
  geom_col() +
  ylab("Number of DE genes") +
  xlab("Clusters") +
  ggtitle("DE > 0.7 log fold change")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_Number_Barplot.pdf"), width = 13, height = 6)

# Number of DE genes per cluster versus cluster size
ggDF1 <- data.frame(table(clusterDeDF$CLUSTER)
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7])
  , NUMBER_CELLS = as.vector(table(centSO@ident)))
p1 <- ggplot(ggDF1, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.2 log fold change")
p2 <- ggplot(ggDF2, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.4 log fold change")
p3 <- ggplot(ggDF3, aes(x = NUMBER_CELLS, y = Freq)) +
  geom_point() +
  xlab("Number of cells in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.7 log fold change")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster versus number of cells per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_NumberVsCells_ScatterPlot.pdf")
  , width = 13, height = 6)

# Number of DE genes per cluster versus mean genes detected per cluster
ggDF1 <- data.frame(table(clusterDeDF$CLUSTER)
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7])
  , nGene = tapply(centSO@meta.data$nGene, centSO@ident, mean))
p1 <- ggplot(ggDF1, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.2 log fold change")
p2 <- ggplot(ggDF2, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.4 log fold change")
p3 <- ggplot(ggDF3, aes(x = nGene, y = Freq)) +
  geom_point() +
  xlab("Mean genes detected in cluster") +
  ylab("Number of DE genes in cluster") +
  ggtitle("DE > 0.7 log fold change")
# plot_grid combine cluster heatmaps
pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nDE genes per cluster versus number of cells per cluster"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "DE_NumberVsnGene_ScatterPlot.pdf")
  , width = 13, height = 6)
################################################################################

### Hierarchical cluster by Seurat cluster mean expression of DE genes

# Top 10 DE genes
clusterDe10DF <- data.frame(
  clusterDeDF %>% group_by(CLUSTER) %>% top_n(10, LOG_FC))

# Subset expression matrix
df <- data.frame(t(centSO@scale.data[row.names(centSO@scale.data) %in% clusterDeDF$GENE, ]))
df$ClusterID <- centSO@ident
df <- aggregate(.~ClusterID, df, mean)
row.names(df) <- df$ClusterID
df <- df[ ,colnames(df) != "ClusterID"]

# Obtain the dendrogram
dend <- as.dendrogram(hclust(
  d = dist(df, method = "euclidean")
  , method = "ward.D2"))
# Extract dend data
dend_data <- dendro_data(dend)
# Exract the hcluster order
hclust_order <- as.numeric(as.character(dend_data$labels$label))

# Plot dendrogram
# ggplot
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
    hjust = 1, angle = 90, size = 3)
ggsave(paste0(outGraph, "hclust_dend.pdf"))

# Heatmap plot
# Use the dendrogram label data to position the labels
# Order genes by dendro order
geneGroupDF <- data.frame(GENE = clusterDe10DF$GENE, GROUP = clusterDe10DF$CLUSTER)
geneGroupDF$GROUP <- factor(geneGroupDF$GROUP, levels = hclust_order)
geneGroupDF <- geneGroupDF[order(geneGroupDF$GROUP), ]
# Plot heatmap
ggL <- lapply(hclust_order, function(cluster){
  # tryCatch(
    Heatmap_By_Cluster(
      geneGroupDF = geneGroupDF
      , exprM = as.matrix(centSO@scale.data)
      , seuratO = centSO
      , clusters = cluster
      , lowerLimit = -1.5
      , upperLimit = 1.5
      , geneOrder = TRUE
    )
    # , error = function(e) NULL)
})
# Remove nulls from ggplot list
ggL <- ggL[! sapply(ggL, is.null)]
# Extract legend
legend <- get_legend(ggL[[1]])
# Duplicate one plot to save labels
ggLabels <- ggL[[1]]
# Remove axis labels
ggL[1:length(ggL)] <- lapply(ggL[1:length(ggL)], function(gg) {
  # gg <- gg + facet_grid(~SEURAT_CLUSTERS, space = "free", scales = "free")
  gg + theme(
    strip.text.y = element_blank()
    , axis.text.y = element_blank()
    , legend.position = "none"
    , axis.title.y = element_blank()
    , axis.ticks.y = element_blank()
    , axis.title.x = element_blank()
    # margin: top, right, bottom, and left
    , plot.margin = unit(c(1, 0.02, 1, 0.02), "cm")
  )
})
# Add genes as labels by plotting and reducing margin to cutoff heatmap
ggLabels <- lbak
ggLabels <- ggLabels + theme(
  strip.text.y = element_blank()
  # , axis.text.y = element_blank()
  , legend.position = "none"
  , axis.title.y = element_blank()
  , axis.ticks.y = element_blank()
  , axis.title.x = element_blank()
  # margin: top, right, bottom, and left
  , plot.margin = unit(c(1, -50, 1, 1), "cm")
)

# Combine individual heatmaps and dendrogram
# Determine relative widths by scaling to log 10 cluster size
df <- data.frame(table(centSO@ident))
df$rel_widths <- as.vector(log((df$Freq + 1), 10)) + 1
df <- df[match(hclust_order, df$Var1), ]
rel_widths <- c(10, df$rel_widths, 1)
# Combine
pg <- plot_grid(plotlist = append(list(ggLabels), ggL), ncol = 19
  , rel_widths = rel_widths, align = 'h', axis = 't')
# Add legend
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.1))
# Add title
title <- paste0(graphCodeTitle
  , "\n\nTop 10 DE genes for each cluster"
  , "\nMean centered, variance scaled normalized expression"
  , "\nCells sorted by hierarchical clustering of Seurat clusters (columns)")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))

ggsave(paste0(outGraph, "hclust_heatmap.png"), width = 12, height = 28)
################################################################################

### Heatmaps

# Top 10 markers for each cluster
clusterDeDF %>% group_by(CLUSTER) %>% top_n(10, LOG_FC) -> top10
geneGroupDF <- data.frame(GENE = top10$GENE, GROUP = top10$CLUSTER)

Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, clusters = c(1:17), lowerLimit = -1.5, upperLimit = 1.5
  , geneOrder = NULL
  , clusterOrder = c(17,16,13,15,8,10,11,7,9,3,14,2,5,6,4,1,0,12))
ggsave(paste0(outGraph, "Top10DE_ExprHeatmap_CentScale.png")
  , width = 13, height = 24)
# 
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
#   , seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
# )
# # ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')

# Split by cluster
ldf <- split(top10, top10$CLUSTER)

# Mean centered variance scaled
# Heatmap for each cluster
pgL <- lapply(names(ldf), function(cl) {
  deDF <- ldf[[cl]]
  # Centered scaled
  geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
  # Heatmaps
  ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
    , seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
  )
  # ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
  # plot_grid combine
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # now add the title
  title <- ggdraw() + draw_label(paste0("Cluster: ", cl))
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  return(pg)
})
# plot_grid combine cluster heatmaps
pg <- plot_grid(plotlist = pgL, ncol = 2, align = 'v', axis = 'r')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nTop 10 DE genes for each cluster"
  , "\nMean centered, variance scaled normalized expression"
  , "\nCells sorted by cluster (columns)"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "ExprHeatmap_CentScale.png")
  , width = 18, height = length(pgL)*2.5)

# No mean centered variance scaled
# Heatmap for each cluster
pgL <- lapply(names(ldf), function(cl) {
  deDF <- ldf[[cl]]
  # Centered scaled
  geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
  # Heatmaps
  ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
    , seuratO = centSO, lowerLimit = 0, upperLimit = 3
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
  )
  # ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
  # plot_grid combine
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # now add the title
  title <- ggdraw() + draw_label(paste0("Cluster: ", clusterID))
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  return(pg)
})
# plot_grid combine cluster heatmaps
pg <- plot_grid(plotlist = pgL, ncol = 2, align = 'v', axis = 'l')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nTop 10 DE genes for each cluster"
  , "\nNormalized expression"
  , "\nCells sorted by cluster (columns)"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "ExprHeatmap_NoCentScale_TEST.png")
  , width = 18, height = length(pgL)*2.5)

# Heatmaps - mean expression
# Mean centered variance scaled
# Top 20 markers for each cluster
clusterDeDF %>% group_by(CLUSTER) %>% top_n(20, LOG_FC) -> top20
# Split by cluster
ldf <- split(top20, top20$CLUSTER)
# Heatmap for each cluster
ggL <- lapply(names(ldf), function(cl) {
  print(cl)
  deDF <- ldf[[cl]]
  gg <- DE_Mean_Heatmap(clusterDeDF = deDF
    , exDF = centSO@scale.data
    , clusterIDs = centSO@ident
    , upLim = 1.5
    , lowLim = -1.5
    , ggtitle = paste0("Cluster: ", cl)
  )
  gg <- gg + theme(axis.text.y = element_text(size = 10)) +
    return(gg)
})
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine cluster heatmaps
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nTop 20 DE genes for each cluster"
  , "\nMean centered, variance scaled normalized expression"
  , "\nMean expression"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "ExprMeanHeatmap_CentScale.png")
  , width = 16, height = length(ggL)*1.5)
  

# # Violin plots of top 10 markers
# pdf(paste0(outGraph, "ViolinPlot_Top10Markers.pdf"), width = 10)
# top10L <- split(top10, top10$cluster)
# lapply(top10L, function(top10cluster) {
#   VlnPlot(centSO, top10cluster$gene, size.use = 0.5)
# })
# dev.off()
# 
# # Feature plot of top 10 markers
# pdf(paste0(outGraph, "FeaturePlot_Top10Markers.pdf"), width = 10)
# top10L <- split(top10, top10$cluster)
# lapply(top10L, function(top10cluster) {
#   FeaturePlot(centSO, top10cluster$gene, cols.use = c("grey","blue")
#     , pt.size = 0.7)
# })
# dev.off()

## To redo expression heatmaps individually for each cluster

ldf <- split(clusterDeDF, clusterDeDF$CLUSTER)

lapply(names(ldf)[1], function(clusterID) {
  
  print(clusterID)
  deDF <- ldf[[clusterID]]
  
  # Expression heatmap of DE genes
  # Centered scaled
  p1 <- DE_Heatmap(clusterDeDF = deDF
    , exDF = centSO@scale.data
    , clusterIDs = centSO@ident
    , upLim = 1.5
    , lowLim = -1.5
    , ggtitle = paste0(
      "\nMean centered, variance scaled, normalized expression"
      , "\n")
  )
  # Not centered scaled
  p2 <- DE_Heatmap(clusterDeDF = deDF
    , exDF = noCentExM
    , clusterIDs = centSO@ident
    , upLim = 3
    , lowLim = 0
    , ggtitle = paste0(
      "\nNormalized expression"
      , "\n")
  )
  p1 <- p1 + theme(axis.text.y = element_blank())
  p2 <- p2 + theme(axis.text.y = element_blank())
  # plot_grid
  pg <- plot_grid(p1, p2, ncol = 2)
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nSignificant (p-value < 0.05) DE genes for cluster ", clusterID
    , "\nCells sorted by cluster (columns)"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  ggsave(paste0(outGraph, "ExprHeatmap_Cluster", clusterID, ".png")
    , width = 12, height = 8)
  
})

# Heatmap - mean expression
lapply(names(ldf), function(clusterID) {
  
  print(clusterID)
  deDF <- ldf[[clusterID]]
  
  deDF <- deDF[1:40, ]
  
  # Expression heatmap of DE genes
  # Centered scaled
  p1 <- DE_Mean_Heatmap(clusterDeDF = deDF
    , exDF = centSO@scale.data
    , clusterIDs = centSO@ident
    , upLim = 1.5
    , lowLim = -1.5
    , ggtitle = paste0(
      "\nMean centered, variance scaled, normalized expression"
      , "\n")
  )
  # Not centered scaled
  p2 <- DE_Mean_Heatmap(clusterDeDF = deDF
    , exDF = noCentExM
    , clusterIDs = centSO@ident
    , upLim = 3
    , lowLim = 0
    , ggtitle = paste0(
      "\nNormalized expression"
      , "\n")
  )
  # plot_grid
  pg <- plot_grid(p1, p2, ncol = 2)
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nSignificant (p-value < 0.05) DE genes for cluster ", clusterID
    , "\nTop 40 DE genes"
    , "\nMean expression"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  ggsave(paste0(outGraph, "ExprMeanHeatmap_Cluster", clusterID, ".png")
    , width = 12, height = 8)
})


################################################################################

### DE genes expression across Kang et al. cortex stages 1-8

## Plot DE genes expression across Kang et al. cortex stages 1-8

ggL <- lapply(sort(unique(clusterDeDF$CLUSTER)), function(cluster){
  # Subset DE data frame to cluster
  specificClusterDeDF <- clusterDeDF[clusterDeDF$CLUSTER == cluster, ]
  # Subset to genes > 0.4 log fold change
  specificClusterDeDF <- specificClusterDeDF[specificClusterDeDF$LOG_FC > 0.4, ]
  # Subset Kang expression matrix to DE genes
  exM <- kangExM[row.names(kangExM) %in% specificClusterDeDF$ENSEMBL, ]
  # Format for ggplot2 and add stage
  df <- melt(exM)
  df$Stage <- kangMtDF$Stage[match(df$Var2, kangMtDF$X)]
  # Plot
  gg <- ggplot(df, aes(x = Stage, y = value)) +
    geom_jitter(size = 0.01, alpha = 0.2) +
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal, 
      fun.args = list(conf.int = 0.95), color = "red", fatten = 0.25) +
    # stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point") +
    # geom_smooth(method = "loess") +
    ylab("Normalized expression") +
    ggtitle(paste0("Cluster ", cluster))
  return(gg)
})
Plot_Grid(ggPlotsL = ggL, ncol = 4, rel_height = 0.15
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster DE genes expression in Kang cortex stage 1-8"
    , "\nDE > 0.4 log fold change"
    , "\nRed points = mean"
    , "\nRed bar = 95% confidence intervals"
    , "\nBlack points = expression of each gene in each sample"))
ggsave(paste0(outGraph, "KangExpr.png"), width = 13, height = 18)


## Correlation to Kang stages

# Convert hgnc symbols to ensembl and leave ensembl IDs unchanged
# (Gene IDs are a mix of hgnc symbols and ensembl IDs for those genes that have
# no hgnc symbol)


Convert_Mixed_GeneSym_EnsID_To_EnsID <- function(ids){
  idx <- match(ids, bmDF$hgnc_symbol)
  ens <- bmDF$ensembl_gene_id[idx]
  ids[! is.na(ens)] <- as.character(ens[! is.na(ens)])
  return(ids)
}

Cluster_Correlation_To_Kang_Stages <- function (topGenes) {
  ll <- lapply(sort(unique(clusterDeDF$CLUSTER)), function(cluster){
    # Subset expression matrix to cells in cluster
    exM <- noCentExM[ ,centSO@ident %in% cluster]
    # Convert hgnc symbols to ensembl and leave ensembl IDs unchanged
    # (Gene IDs are a mix of hgnc symbols and ensembl IDs for those genes that have
    # no hgnc symbol)
    row.names(exM) <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(exM))
    # Subset to genes above mean expression threshold
    rMns <- rowMeans(exM)
    rMns <- sort(rMns, decreasing = TRUE)
    ids <- names(rMns)[1:topGenes]
    ssKangExM <- kangExM[row.names(kangExM) %in% ids, ]
    df <- melt(ssKangExM)
    df$Stage <- kangMtDF$Stage[match(df$Var2, kangMtDF$X)]
    df <- aggregate(value~Stage+Var1, df, mean)
    df <- dcast(df, Var1~Stage)
    df$Cluster_Mean_Expression <- rMns[match(df$Var1, names(rMns))]
    scor <- apply(df[ ,-1], 2, function(col){
      cor(df$Cluster_Mean_Expression, col, method = "spearman")
    })
    return(scor)
  })
  df <- do.call("cbind", ll)
  df <- df[-9, ]
  colnames(df) <- sort(unique(clusterDeDF$CLUSTER))
  df <- melt(df)
  return(df)
}

Plot_Cluster_Correlation_To_Kang_Stages <- function (ggDF, topGenes) {
  ggplot(ggDF, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", space = "Lab"
      , name = "Spearman") +
    geom_text(aes(label = round(value, 2))) +
    scale_x_continuous(breaks = unique(ggDF$Var1)) +
    scale_y_continuous(breaks = unique(ggDF$Var2)) +
    xlab("Stage") +
    ylab("Cluster") +
    ggtitle(paste0("Used top ", topGenes, " expressed genes in cluster"))
}

df <- Cluster_Correlation_To_Kang_Stages(topGenes = 2500)
gg1 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 2500)

df <- Cluster_Correlation_To_Kang_Stages(topGenes = 5000)
gg2 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 5000)

df <- Cluster_Correlation_To_Kang_Stages(topGenes = 10000)
gg3 <- Plot_Cluster_Correlation_To_Kang_Stages(ggDF = df, topGenes = 10000)

Plot_Grid(ggPlotsL = list(gg1, gg2, gg3), ncol = 3, rel_height = 0.15
  , title = paste0(graphCodeTitle
    , "\n\nCorrelation of Seurat cluster mean expression profile to"
    , "\nKang stage expression")
)
ggsave(paste0(outGraph, "KangCorrelation.png"), width = 16, height = 9)
################################################################################

### DE genes across Miller zones

## Expression across Miller zones

ggL <- lapply(sort(unique(clusterDeDF$CLUSTER)), function(cluster){
  # Subset DE data frame to cluster
  specificClusterDeDF <- clusterDeDF[clusterDeDF$CLUSTER == cluster, ]
  # Subset to genes > 0.4 log fold change
  specificClusterDeDF <- specificClusterDeDF[specificClusterDeDF$LOG_FC > 0.4, ]
  # Subset Kang expression matrix to DE genes
  exM <- millerExDF[row.names(millerExDF) %in% specificClusterDeDF$ENSEMBL, ]
  # Format for ggplot2 and add stage
  df <- melt(exM)
  df$Zone <- millerMtDF$Zone[match(df$variable, millerMtDF$well_id)]
  df$Zone <- factor(df$Zone, levels = c("VZ", "SZi", "SZo", "IZ", "SP", "CPi", "CPo", "MZ", "SG"))
  df <- df[! is.na(df$Zone), ]
  # Plot
  gg <- ggplot(df, aes(x = Zone, y = value)) +
    geom_jitter(size = 0.01, alpha = 0.2) +
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal, 
      fun.args = list(conf.int = 0.95), color = "lightcoral", alpha = 0.75) +
    stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point") +
    # geom_smooth(method = "loess") +
    ylab("Normalized expression") +
    ggtitle(paste0("Cluster ", cluster))
  return(gg)
})
Plot_Grid(ggPlotsL = ggL, ncol = 4, rel_height = 0.15
  , title = paste0(graphCodeTitle
    , "\n\nSeurat cluster DE genes expression in Miller zones"
    , "\nDE > 0.4 log fold change"
    , "\nRed points = mean"
    , "\nBlue bar = 95% confidence intervals"
    , "\nBlack points = expression of each gene in each sample"))
ggsave(paste0(outGraph, "MillerExpr.png"), width = 13, height = 18)
################################################################################




