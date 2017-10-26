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
source("Function_Library.R")

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
clusterID <- (as.numeric(args[1]) - 1)
print(paste0("Cluster ID: ", clusterID))

## Inputs
# Seurat
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO_TEST.Robj")
# Biomart to add ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterDE.R"
outGraph <- "../analysis/graphs/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"
outTable <- "../analysis/tables/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"
outData <- "../analysis/Seurat_ClusterDE_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_ClusterDE_DS2-11_"

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

# # Filter expression matrix by:
# # percent of cells in cluster gene is expressed in
# # fold change of gene in cluster versus all other cells
# DE_Filters_ExpMatrix <- function(minPercent, foldChange, clusterID) {
#   # Expressed > 0 counts in > X% of cells in cluster
#   # Subset expression matrix to cluster
#   cdf <- as.matrix(centSO@data)[ ,centSO@ident == clusterID]
#   # Expressed > 0 counts in > X% of cells in cluster
#   idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (minPercent / 100)
#   print(paste0("Genes expressed in > ", minPercent, "% of cells in cluster"))
#   print(table(idxp))
#   
#   # Fold change > Y of gene in cluster versus all other cells
#   # Subset expression matrix to cluster
#   cdf <- noCentExM[ ,centSO@ident == clusterID]
#   # Subset expression matrix to all other cells
#   ndf <- noCentExM[ , ! centSO@ident == clusterID]
#   # Fold change
#   v1 <- rowMeans(cdf) - rowMeans(ndf)
#   idxf <- v1 > foldChange
#   print(paste0("Genes > ", foldChange, " fold change in cluster versus all other cells"))
#   print(table(idxf))
#   
#   # Filter exDF
#   exDF <- as.matrix(centSO@data[idxp & idxf, ])
#   return(exDF)
# }
# 
# ## Function: DE Linear model
# # termsDF:
# # ExpCondition RIN.y     Seq.PC1      Seq.PC2
# # 1            CP   8.4  0.04792498 -0.090448567
# # 2            CP   8.1  0.53502697 -0.287629654
# # 3            CP   8.0  0.18824922 -0.155651102
# # 4            VZ   8.4  0.02529722 -0.100858264
# # 5            VZ   8.7  0.45139297  0.856908177
# # 6            VZ   9.1  0.27861748 -0.248868277
# # mod: "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2"
# DE_Linear_Model <- function (exDatDF, termsDF, mod) {
#   
#   lmmod <- apply(as.matrix(exDatDF), 1
#     , function(y) {
#       mod <- as.formula(mod)
#       lm(mod, data = termsDF)})
#   
#   coefmat <- matrix(NA, nrow = nrow(exDatDF)
#     , ncol = length(coef(lmmod[[1]])))
#   pvalmat <- matrix(NA, nrow = nrow(exDatDF)
#     , ncol = length(summary(lmmod[[1]])[[4]][ ,4]))
#   colnames(coefmat) <- names(coef(lmmod[[1]]))
#   rownames(coefmat) <- rownames(exDatDF)
#   colnames(pvalmat) <- names(summary(lmmod[[1]])[[4]][ ,4])
#   rownames(pvalmat) <- rownames(exDatDF)
#   for (i in 1:nrow(exDatDF)) {
#     if (i%%100 == 0) {cat(".")}
#     coefmat[i, ] <- coef(lmmod[[i]])
#     pvalmat[i, ] <- summary(lmmod[[i]])[[4]][ ,4]
#   }
#   deCoefPvalLM <- list(coefmat = coefmat, pvalmat = pvalmat)
#   return(deCoefPvalLM)
# }
# 
# # Format output of linear model into data frame
# Format_DE <- function (deLM, clusterID) {
#   # Combine log2 fold changes, p-values
#   deDF <- data.frame(GENE = row.names(deLM$coefmat)
#     , LOG_FC = deLM$coefmat[ ,2]
#     , PVALUE = deLM$pvalmat[ ,2])
#   # Order by pvalue
#   deDF <- deDF[order(deDF$PVALUE), ]
#   # Add cluster ID
#   deDF$CLUSTER <- clusterID
#   # Percent of cells in cluster expressing gene > 0 counts
#   cdf <- as.matrix(centSO@data)[row.names(centSO@data) %in% deDF$GENE, centSO@ident == clusterID]
#   deDF$PERCENT_CLUSTER <- (rowSums(cdf > 0) / ncol(cdf)) * 100
#   # Percent of all cells expressing gene > 0 counts
#   deDF$PERCENT_ALL <- (rowSums(as.matrix(centSO@data)[row.names(centSO@data) %in% deDF$GENE, ] > 0)
#     / ncol(centSO@data[row.names(centSO@data) %in% deDF$GENE, ])) * 100
#   # Order by log fold change
#   deDF <- deDF[order(-deDF$LOG_FC), ]
#   return(deDF)
# }

# # Heatmap of expression
# DE_Heatmap <- function(clusterDeDF, exDF, clusterIDs, ggtitle, upLim, lowLim) {
#   # Subset to genes, merge to keep duplicate genes (from more than one cluster)
#   ggDF <- merge(clusterDeDF, exDF, by.x = "GENE", by.y = "row.names"
#     , all.x = TRUE)
#   # For gene ordering by cluster
#   idx <- match(rev(clusterDeDF$GENE), ggDF$GENE)
#   ggDF <- ggDF[idx, ]
#   row.names(ggDF) <- paste0(length(ggDF$GENE):1, "_", ggDF$GENE)
#   ggDF <- ggDF[ ,-c(1:ncol(clusterDeDF))]
#   # Order by clustering
#   idx <- match(colnames(ggDF), names(sort(clusterIDs)))
#   ggDF <- ggDF[ ,idx]
#   # Format for ggplot2
#   ggDF <- as.matrix(ggDF)
#   ggDF <- melt(ggDF)
#   # Add clusters
#   idx <- match(ggDF$Var2, names(clusterIDs))
#   ggDF$CLUSTERS <- clusterIDs[idx]
#   # Set sample order by clustering
#   ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(clusterIDs)))
#   # Set expression limits
#   ggDF$value[ggDF$value > upLim] <- upLim
#   ggDF$value[ggDF$value < lowLim] <- lowLim
#   # ggplot
#   ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
#     geom_tile() +
#     facet_grid(~CLUSTERS, scales = "free") +
#     # facet_grid(~CLUSTERS, space = "free", scales = "free") +
#     # scale_fill_gradient2(name = "Normalized\nexpression"
#     #   , high = "#d7191c", low = "white") +
#     scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
#     scale_fill_distiller(name = "Normalized\nexpression", type = "div"
#       , palette = 5, direction = -1) +
#     theme_bw() +
#     theme(axis.text.x = element_blank()) +
#     theme(axis.ticks = element_blank()) +
#     theme(text = element_text(size = 12)) +
#     # theme(axis.text.y = element_blank()) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     ylab("Genes") +
#     xlab("Cells") +
#     ggtitle(ggtitle)
# }

# Heatmap of mean expression per cluster
DE_Mean_Heatmap <- function(
  clusterDeDF, exDF, clusterIDs, ggtitle, upLim, lowLim) {
  
  # Subset to DE genes
  ggDF <- exDF[match(clusterDeDF$GENE, row.names(exDF)), ]
  
  # Change column names from Cell IDs to Cluster ID
  colnames(ggDF) <- clusterIDs
  # Melt
  ggDF <- melt(ggDF)
  # Mean expression per cluster
  ggDF <- aggregate(value~Var2+Var1, mean, data = ggDF)
  # Order genes by fold change
  ggDF$Var1 <- factor(ggDF$Var1, rev(levels(ggDF$Var1)))
  # Set expression limits
  ggDF$value[ggDF$value > upLim] <- upLim
  ggDF$value[ggDF$value < lowLim] <- lowLim
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    # scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(lowLim, upLim)) +
    scale_x_continuous(breaks = unique(ggDF$Var2), labels = unique(ggDF$Var2)) +
    theme_bw() +
    # theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    # theme(axis.text.y = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(ggtitle)
}
################################################################################

### Finding differentially expressed genes (cluster biomarkers)

print(paste0("### Finding differentially expressed genes for cluster ", clusterID))

# Filter cells
df <- DE_Filters_ExpMatrix(centSO
  , minPercent = 10, foldChange = 0.2, clusterID = clusterID)

# DE Linear model
termsDF <- centSO@meta.data[c("nUMI", "librarylab", "individual")]
# Add term TRUE/FALSE cell is in cluster
termsDF$cluster <- FALSE
termsDF$cluster[centSO@ident == clusterID] <- TRUE
mod <- "y ~ cluster+nUMI+librarylab+individual"
deLM <- DE_Linear_Model(exDatDF = df, termsDF = termsDF, mod = mod)

# Format LM output into data frame
deDF <- Format_DE(deLM, centSO, clusterID)
print(head(deDF))

# FDR correct
# NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
corrected <- fdrtool(deDF$PVALUE, statistic = "pvalue", plot = FALSE)
deDF$FDR <- corrected$lfdr
# Check
table(deDF$PVALUE < 0.05)
table(deDF$FDR < 0.05)
print(head(deDF))
# Plot
# ggplot histogram p-values
p1 <- ggplot(deDF, aes(x = PVALUE)) +
  geom_histogram() +
  xlab("P-value") +
  ylab("Count") +
  ggtitle("P-value") +
  theme(plot.title = element_text(hjust = 0.5))
# ggplot histogram FDR
p2 <- ggplot(deDF, aes(x = FDR)) +
  geom_histogram() +
  xlab("FDR corrected p-value") +
  ylab("Count") +
  ggtitle("Benjamini Hochberg corrected p-value") +
  theme(plot.title = element_text(hjust = 0.5))
# plot_grid
pg <- plot_grid(p1, p2, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nHistograms of p-values before and after Benjamini Hochberg correction"
  , "\nCluster: ", clusterID))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "Pvalue_Histogram_Cluster", clusterID, ".pdf")
  , width = 8, height = 7)

# v1 <- c("ETV1", "FOXP1", "TBR1", "BCL11B", "SOX5", "FEZF2")

# Write to tab delimited table
print("Writing DE table to text file")
write.table(x = deDF
  , file = paste0(outTable, "Cluster", clusterID, "_Vs_All_Clusters.txt")
  , sep = "\t", quote = FALSE, row.names = FALSE)

## Expression heatmaps of DE genes

# Centered scaled
geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
# Heatmaps
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
)
# Remove y-axis labels
ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
  , "\nMean centered, variance scaled, normalized expression"
  , "\nCells sorted by cluster (columns)"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale_Cluster", clusterID, ".png")
  , width = 12, height = 8)

# Not centered scaled
geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
# Heatmaps
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = noCentExM
  , seuratO = centSO, lowerLimit = 0, upperLimit = 3
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
)
# Remove y-axis labels
ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
  , "\nNormalized expression"
  , "\nCells sorted by cluster (columns)"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "ExprHeatmap_Norm_Cluster", clusterID, ".png")
  , width = 12, height = 8)
################################################################################

### DE Heatmaps, top 20 DE table, and save seurat object

# If cluster DE is done for all clusters compile DE tables and make DE heatmaps
if (clusterID == max(as.numeric(as.character(unique(centSO@ident))))) {
  
  print("### DE calculated for last cluster, compiling DE table, plotting metrics and heatmaps...")
  
  # Loop through DE text files and compile into one table
  ldf <- lapply(sort(unique(centSO@ident)), function(clid) {
    df <- read.table(
      paste0(outTable, "Cluster", clid, "_Vs_All_Clusters.txt")
      , header = TRUE)
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
  
  # Top 20 DE genes
  clusterDe20DF <- data.frame(
    clusterDeDF %>% group_by(CLUSTER) %>% top_n(20, LOG_FC))
  # Write out as tab delimited
  write.table(x = clusterDe20DF
    , file = paste0(outTable, "ClusterX_Vs_All_Clusters_Top20.txt")
    , sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save seurat object with cluster DE data frame
  save(centSO, noCentExM, metDF, clusterDeDF, file = paste0(outData, "seuratO.Robj"))
  
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
    # Subset expression matrix to cluster
    cdf <- noCentExM[ ,centSO@ident == clusterID]
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
    print(head(ggDF))
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
  
  ## Heatmaps
  
  # Top 10 markers for each cluster
  clusterDeDF %>% group_by(CLUSTER) %>% top_n(10, LOG_FC) -> top10
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
    title <- ggdraw() + draw_label(paste0("Cluster: ", clusterID))
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
  
}

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