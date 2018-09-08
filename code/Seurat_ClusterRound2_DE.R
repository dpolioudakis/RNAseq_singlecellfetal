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

# Seurat clustering round 2 object
load(paste0("../analysis/analyzed_data/Seurat_ClusterRound2/DS2-11"
  , "/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40"
  , "/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_Cluster"
  , clusterID, "_seuratO.Robj"))

# Biomart to add ensembl IDs
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterRound2_DE.R"
out_sub_path <- paste0(
  "Seurat_ClusterRound2"
  , "/DE"
  , "/DS2-11"
  , "/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40"
  , "/VarGenes"
  , "/RegNumiLibBrain"
  , "/PC1-40"
)
outGraph <- paste0(
  "../analysis/graphs/", out_sub_path, "/Seurat_ClusterRound2_")
outTable <- paste0(
  "../analysis/tables/", out_sub_path, "/Seurat_ClusterRound2_")
outData <- paste0(
  "../analysis/analyzed_data/", out_sub_path, "/Seurat_ClusterRound2_")

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

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

### Differentially expressed genes for each cluster vs other cells in cluster

print(paste0(
  "### Finding differentially expressed genes for cluster ", clusterID)
)

Run_DE <- function(sub_cluster_ID){

  print("Run_DE")
  print(paste0("Calculating DE for subcluster: ", sub_cluster_ID))

  # Filter cells
  exDF <- DE_Filters_ExpMatrix(
    so = so, minPercent = 10, clusterID = sub_cluster_ID
  )

  # DE Linear model
  termsDF <- so@meta.data[c("nUMI", "librarylab", "individual")]
  termsDF <- termsDF[row.names(termsDF) %in% colnames(exDF), ]

  # Add term TRUE/FALSE cell is in cluster
  termsDF$cluster <- FALSE
  termsDF$cluster[so@ident == sub_cluster_ID] <- TRUE
  mod <- "y ~ cluster+nUMI+librarylab+individual"
  deLM <- DE_Linear_Model(exDatDF = exDF, termsDF = termsDF, mod = mod)

  # Format LM output into data frame
  deDF <- Format_DE(deLM, so, sub_cluster_ID)

  # Add ensembl
  deDF$Ensembl <- Convert_Mixed_GeneSym_EnsID_To_EnsID(as.character(deDF$Gene))

  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method = "BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)

  # Format
  # Cluster and subcluster columns
  names(deDF)[names(deDF) == "Cluster"] <- "Subcluster"
  deDF$Cluster <- clusterID
  # Order columns
  deDF <- deDF[ ,c("Cluster", "Subcluster", "Gene", "Ensembl"
    , "Log2_Fold_Change", "Pvalue", "FDR", "Percent_Cluster", "Percent_All")]

  # Write to tab delimited table
  print("Writing DE table as csv")
  write.csv(x = deDF
    , file = paste0(outTable, "Cluster", clusterID, "_Subcluster", sub_cluster_ID, "_Vs_All_Subcluster.csv")
    , quote = FALSE, row.names = FALSE)
}

# Run DE
lapply(unique(so@ident), function(sub_cluster_ID){
  Run_DE(sub_cluster_ID)
})
################################################################################

# ### Finding differentially expressed genes (cluster biomarkers)
#
# print("### Finding differentially expressed genes for each cluster")
#
#
# ldf <- lapply(lso, function(so) {
#
#   # so <- lso[[6]]
#   # DE for each cluster
#   ldf <- lapply(sort(unique(so@ident)), function(cluster) {
#     # Filter cells
#     df <- DE_Filters_ExpMatrix(
#       so = so, minPercent = 10, clusterID = cluster)
#
#     # DE Linear model
#     termsDF <- so@meta.data[c("nUMI", "librarylab", "individual")]
#     # Add term TRUE/FALSE cell is in cluster
#     termsDF$cluster <- FALSE
#     termsDF$cluster[so@ident == cluster] <- TRUE
#     print(head(termsDF))
#     mod <- "y ~ cluster+nUMI+librarylab+individual"
#     deLM <- DE_Linear_Model(exDatDF = df, termsDF = termsDF, mod = mod)
#
#     # Format LM output into data frame
#     deDF <- Format_DE(deLM, so, cluster)
#     print(head(deDF))
#
#     # FDR correct
#     # NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
#     corrected <- fdrtool(deDF$PVALUE, statistic = "pvalue", plot = FALSE)
#     deDF$FDR <- corrected$lfdr
#     # Check
#     print(table(deDF$PVALUE < 0.05))
#     print(table(deDF$FDR < 0.05))
#     print(head(deDF))
#
#     return(deDF)
#   })
#   df <- do.call("rbind", ldf)
#   return(df)
# })
#
# save.image(paste0(outData, "Workspace.RData"))


#
#
#
# # Write to tab delimited table
# print("Writing DE table to text file")
# write.table(x = deDF
#   , file = paste0(outTable, "Cluster", clusterID, "_Vs_All_Clusters.txt")
#   , sep = "\t", quote = FALSE, row.names = FALSE)
#
# ## Expression heatmaps of DE genes
#
# # Centered scaled
# geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
# # Heatmaps
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = so@scale.data
#   , seuratO = so, lowerLimit = -1.5, upperLimit = 1.5
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
# )
# # Remove y-axis labels
# ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
#   , "\nMean centered, variance scaled, normalized expression"
#   , "\nCells sorted by cluster (columns)"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# ggsave(paste0(outGraph, "ExprHeatmap_NormCentScale_Cluster", clusterID, ".png")
#   , width = 12, height = 8)
#
# # Not centered scaled
# geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
# # Heatmaps
# ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = noCentExM
#   , seuratO = so, lowerLimit = 0, upperLimit = 3
#   , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
# )
# # Remove y-axis labels
# ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
# # plot_grid combine
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nSignificant (FDR < 0.05) DE genes for cluster ", clusterID
#   , "\nNormalized expression"
#   , "\nCells sorted by cluster (columns)"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# ggsave(paste0(outGraph, "ExprHeatmap_Norm_Cluster", clusterID, ".png")
#   , width = 12, height = 8)
# ################################################################################
#
# ### DE Heatmaps, top 20 DE table, and save seurat object
#
# # If cluster DE is done for all clusters compile DE tables and make DE heatmaps
# if (clusterID == max(as.numeric(as.character(unique(so@ident))))) {
#
#   print("### DE calculated for last cluster, compiling DE table, plotting metrics and heatmaps...")
#
#   # Loop through DE text files and compile into one table
#   ldf <- lapply(sort(unique(so@ident)), function(clid) {
#     df <- read.table(
#       paste0(outTable, "Cluster", clid, "_Vs_All_Clusters.txt")
#       , header = TRUE)
#     df <- df[order(-df$LOG_FC), ]
#   })
#   clusterDeDF <- do.call("rbind", ldf)
#   head(clusterDeDF)
#   # Add ensembl IDs
#   clusterDeDF$ENSEMBL <- bmDF$ensembl_gene_id[
#     match(clusterDeDF$GENE, bmDF$hgnc_symbol)]
#   # Fill in IDs from GENE column that were kept as ensembl ID b/c no gene sym
#   clusterDeDF$ENSEMBL[is.na(clusterDeDF$ENSEMBL)] <-
#     clusterDeDF$GENE[is.na(clusterDeDF$ENSEMBL)]
#   # Check
#   length(grep("ENSG", clusterDeDF$ENSEMBL))
#   nrow(clusterDeDF)
#   # # Filter FDR < 0.05
#   # deDF <- deDF[deDF$FDR < 0.05, ]
#
#   # Write out as tab delimited
#   write.table(x = clusterDeDF
#     , file = paste0(outTable, "ClusterX_Vs_All_Clusters.txt")
#     , sep = "\t", quote = FALSE, row.names = FALSE)
#
#   # Top 20 DE genes
#   clusterDe20DF <- data.frame(
#     clusterDeDF %>% group_by(CLUSTER) %>% top_n(20, LOG_FC))
#   # Write out as tab delimited
#   write.table(x = clusterDe20DF
#     , file = paste0(outTable, "ClusterX_Vs_All_Clusters_Top20.txt")
#     , sep = "\t", quote = FALSE, row.names = FALSE)
#
#   # Save seurat object with cluster DE data frame
#   save(so, noCentExM, metDF, clusterDeDF, file = paste0(outData, "seuratO.Robj"))
#
#   # Violin plots of fold changes by cluster
#   ggDF <- clusterDeDF
#   ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
#   ggplot(ggDF, aes(x = CLUSTER, y = LOG_FC)) +
#     geom_violin(aes(fill = CLUSTER)) +
#     geom_jitter(size = 0.05) +
#     ylab("Log fold change") +
#     xlab("Cluster") +
#     ggtitle(paste0(graphCodeTitle
#       , "\n\nFold changes for DE genes by cluster"
#       , "\n"))
#   ggsave(paste0(outGraph, "DE_Violin.png"), width = 12, height = 8)
#
#   # Histograms of fold changes by cluster
#   ggDF <- clusterDeDF
#   ggDF$CLUSTER <- as.factor(ggDF$CLUSTER)
#   ggplot(ggDF, aes(x = LOG_FC)) +
#     geom_histogram() +
#     facet_wrap(~CLUSTER) +
#     ylab("Counts") +
#     xlab("Log fold change") +
#     ggtitle(paste0(graphCodeTitle
#       , "\n\nHistogram of fold changes for DE genes by cluster"
#       , "\n"))
#   ggsave(paste0(outGraph, "DE_Histogram.png"), width = 12, height = 8)
#
#   # MA plots for each cluster
#   ggL <- lapply(unique(clusterDeDF$CLUSTER), function(clusterID) {
#     # Subset expression matrix to cluster
#     cdf <- noCentExM[ ,so@ident == clusterID]
#     # Subset expression matrix to all other cells
#     ndf <- noCentExM[ , ! so@ident == clusterID]
#     # Fold change
#     lfc <- rowMeans(cdf) - rowMeans(ndf)
#     # Mean expression
#     mn <- rowMeans(noCentExM)
#     # Combine in DF for ggplot
#     ggDF <- data.frame(LOG_FOLD_CHANGE = lfc, MEAN_EXPRESSION = mn
#       , GENE = names(lfc))
#     ggDF$DE_GENE <- row.names(ggDF) %in% clusterDeDF$GENE[clusterDeDF$CLUSTER == clusterID]
#     print(head(ggDF))
#     # ggplot
#     gg <- ggplot(ggDF, aes(x = MEAN_EXPRESSION, y = LOG_FOLD_CHANGE)) +
#       geom_point(size = 0.1, alpha = 0.25, aes(color = DE_GENE)) +
#       scale_color_manual(values = c("black", "red")) +
#       theme(legend.position = "none") +
#       xlab("Mean normalized expression") +
#       ylab("Log fold change") +
#       ggtitle(clusterID)
#     return(gg)
#   })
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nMA plots: Cells in cluster Vs all other cells"
#     , "\nColor indicates > 0.2 log fold change and FDR > 0.05"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
#   ggsave(paste0(outGraph, "MAplot.png"), width = 12, height = 4+length(ggL))
#
#   # Number of DE genes per cluster barplots
#   ggDF1 <- data.frame(table(clusterDeDF$CLUSTER))
#   ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4]))
#   ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7]))
#   p1 <- ggplot(ggDF1, aes(x = Var1, y = Freq)) +
#     geom_col() +
#     ylab("Number of DE genes") +
#     xlab("Clusters") +
#     ggtitle("DE > 0.2 log fold change")
#   p2 <- ggplot(ggDF2, aes(x = Var1, y = Freq)) +
#     geom_col() +
#     ylab("Number of DE genes") +
#     xlab("Clusters") +
#     ggtitle("DE > 0.4 log fold change")
#   p3 <- ggplot(ggDF3, aes(x = Var1, y = Freq)) +
#     geom_col() +
#     ylab("Number of DE genes") +
#     xlab("Clusters") +
#     ggtitle("DE > 0.7 log fold change")
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nDE genes per cluster"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
#   ggsave(paste0(outGraph, "DE_Number_Barplot.pdf"), width = 13, height = 6)
#
#   # Number of DE genes per cluster versus cluster size
#   ggDF1 <- data.frame(table(clusterDeDF$CLUSTER)
#     , NUMBER_CELLS = as.vector(table(so@ident)))
#   ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4])
#     , NUMBER_CELLS = as.vector(table(so@ident)))
#   ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7])
#     , NUMBER_CELLS = as.vector(table(so@ident)))
#   p1 <- ggplot(ggDF1, aes(x = NUMBER_CELLS, y = Freq)) +
#     geom_point() +
#     xlab("Number of cells in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.2 log fold change")
#   p2 <- ggplot(ggDF2, aes(x = NUMBER_CELLS, y = Freq)) +
#     geom_point() +
#     xlab("Number of cells in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.4 log fold change")
#   p3 <- ggplot(ggDF3, aes(x = NUMBER_CELLS, y = Freq)) +
#     geom_point() +
#     xlab("Number of cells in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.7 log fold change")
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nDE genes per cluster versus number of cells per cluster"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
#   ggsave(paste0(outGraph, "DE_NumberVsCells_ScatterPlot.pdf")
#     , width = 13, height = 6)
#
#   # Number of DE genes per cluster versus mean genes detected per cluster
#   ggDF1 <- data.frame(table(clusterDeDF$CLUSTER)
#     , nGene = tapply(so@meta.data$nGene, so@ident, mean))
#   ggDF2 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.4])
#     , nGene = tapply(so@meta.data$nGene, so@ident, mean))
#   ggDF3 <- data.frame(table(clusterDeDF$CLUSTER[clusterDeDF$LOG_FC > 0.7])
#     , nGene = tapply(so@meta.data$nGene, so@ident, mean))
#   p1 <- ggplot(ggDF1, aes(x = nGene, y = Freq)) +
#     geom_point() +
#     xlab("Mean genes detected in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.2 log fold change")
#   p2 <- ggplot(ggDF2, aes(x = nGene, y = Freq)) +
#     geom_point() +
#     xlab("Mean genes detected in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.4 log fold change")
#   p3 <- ggplot(ggDF3, aes(x = nGene, y = Freq)) +
#     geom_point() +
#     xlab("Mean genes detected in cluster") +
#     ylab("Number of DE genes in cluster") +
#     ggtitle("DE > 0.7 log fold change")
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(p1, p2, p3, ncol = 3, align = 'v', axis = 'r')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nDE genes per cluster versus number of cells per cluster"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
#   ggsave(paste0(outGraph, "DE_NumberVsnGene_ScatterPlot.pdf")
#     , width = 13, height = 6)
#
#   ## Heatmaps
#
#   # Top 10 markers for each cluster
#   clusterDeDF %>% group_by(CLUSTER) %>% top_n(10, LOG_FC) -> top10
#   # Split by cluster
#   ldf <- split(top10, top10$CLUSTER)
#
#   # Mean centered variance scaled
#   # Heatmap for each cluster
#   pgL <- lapply(names(ldf), function(cl) {
#     deDF <- ldf[[cl]]
#     # Centered scaled
#     geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
#     # Heatmaps
#     ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = so@scale.data
#       , seuratO = so, lowerLimit = -1.5, upperLimit = 1.5
#       , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
#     )
#     # ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
#     # plot_grid combine
#     pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
#     # now add the title
#     title <- ggdraw() + draw_label(paste0("Cluster: ", clusterID))
#     # rel_heights values control title margins
#     pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
#     return(pg)
#   })
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = pgL, ncol = 2, align = 'v', axis = 'r')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nTop 10 DE genes for each cluster"
#     , "\nMean centered, variance scaled normalized expression"
#     , "\nCells sorted by cluster (columns)"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
#   ggsave(paste0(outGraph, "ExprHeatmap_CentScale.png")
#     , width = 18, height = length(pgL)*2.5)
#
#   # No mean centered variance scaled
#   # Heatmap for each cluster
#   pgL <- lapply(names(ldf), function(cl) {
#     deDF <- ldf[[cl]]
#     # Centered scaled
#     geneGroupDF <- data.frame(GENE = deDF$GENE, GROUP = "")
#     # Heatmaps
#     ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = so@scale.data
#       , seuratO = so, lowerLimit = 0, upperLimit = 3
#       , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
#     )
#     # ggL <- lapply(ggL, function(gg) {gg + theme(axis.text.y = element_blank())})
#     # plot_grid combine
#     pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
#     # now add the title
#     title <- ggdraw() + draw_label(paste0("Cluster: ", clusterID))
#     # rel_heights values control title margins
#     pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
#     return(pg)
#   })
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = pgL, ncol = 2, align = 'v', axis = 'l')
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nTop 10 DE genes for each cluster"
#     , "\nNormalized expression"
#     , "\nCells sorted by cluster (columns)"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
#   ggsave(paste0(outGraph, "ExprHeatmap_NoCentScale_TEST.png")
#     , width = 18, height = length(pgL)*2.5)
#
#   # Heatmaps - mean expression
#   # Mean centered variance scaled
#   # Top 20 markers for each cluster
#   clusterDeDF %>% group_by(CLUSTER) %>% top_n(20, LOG_FC) -> top20
#   # Split by cluster
#   ldf <- split(top20, top20$CLUSTER)
#   # Heatmap for each cluster
#   ggL <- lapply(names(ldf), function(cl) {
#     print(cl)
#     deDF <- ldf[[cl]]
#     gg <- DE_Mean_Heatmap(clusterDeDF = deDF
#       , exDF = so@scale.data
#       , clusterIDs = so@ident
#       , upLim = 1.5
#       , lowLim = -1.5
#       , ggtitle = paste0("Cluster: ", cl)
#     )
#     gg <- gg + theme(axis.text.y = element_text(size = 10)) +
#       return(gg)
#   })
#   # extract the legend from one of the plots
#   legend <- get_legend(ggL[[1]])
#   # Remove legends from plots
#   ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
#   # plot_grid combine cluster heatmaps
#   pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
#   # add the legend to the row we made earlier. Give it one-third of the width
#   # of one plot (via rel_widths).
#   pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nTop 20 DE genes for each cluster"
#     , "\nMean centered, variance scaled normalized expression"
#     , "\nMean expression"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
#   ggsave(paste0(outGraph, "ExprMeanHeatmap_CentScale.png")
#     , width = 16, height = length(ggL)*1.5)
#
# }
#
# # # Violin plots of top 10 markers
# # pdf(paste0(outGraph, "ViolinPlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   VlnPlot(so, top10cluster$gene, size.use = 0.5)
# # })
# # dev.off()
# #
# # # Feature plot of top 10 markers
# # pdf(paste0(outGraph, "FeaturePlot_Top10Markers.pdf"), width = 10)
# # top10L <- split(top10, top10$cluster)
# # lapply(top10L, function(top10cluster) {
# #   FeaturePlot(so, top10cluster$gene, cols.use = c("grey","blue")
# #     , pt.size = 0.7)
# # })
# # dev.off()
#
# ## To redo expression heatmaps individually for each cluster
#
# ldf <- split(clusterDeDF, clusterDeDF$CLUSTER)
#
# lapply(names(ldf)[1], function(clusterID) {
#
#   print(clusterID)
#   deDF <- ldf[[clusterID]]
#
#   # Expression heatmap of DE genes
#   # Centered scaled
#   p1 <- DE_Heatmap(clusterDeDF = deDF
#     , exDF = so@scale.data
#     , clusterIDs = so@ident
#     , upLim = 1.5
#     , lowLim = -1.5
#     , ggtitle = paste0(
#       "\nMean centered, variance scaled, normalized expression"
#       , "\n")
#   )
#   # Not centered scaled
#   p2 <- DE_Heatmap(clusterDeDF = deDF
#     , exDF = noCentExM
#     , clusterIDs = so@ident
#     , upLim = 3
#     , lowLim = 0
#     , ggtitle = paste0(
#       "\nNormalized expression"
#       , "\n")
#   )
#   p1 <- p1 + theme(axis.text.y = element_blank())
#   p2 <- p2 + theme(axis.text.y = element_blank())
#   # plot_grid
#   pg <- plot_grid(p1, p2, ncol = 2)
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nSignificant (p-value < 0.05) DE genes for cluster ", clusterID
#     , "\nCells sorted by cluster (columns)"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
#   ggsave(paste0(outGraph, "ExprHeatmap_Cluster", clusterID, ".png")
#     , width = 12, height = 8)
#
# })
#
# # Heatmap - mean expression
# lapply(names(ldf), function(clusterID) {
#
#   print(clusterID)
#   deDF <- ldf[[clusterID]]
#
#   deDF <- deDF[1:40, ]
#
#   # Expression heatmap of DE genes
#   # Centered scaled
#   p1 <- DE_Mean_Heatmap(clusterDeDF = deDF
#     , exDF = so@scale.data
#     , clusterIDs = so@ident
#     , upLim = 1.5
#     , lowLim = -1.5
#     , ggtitle = paste0(
#       "\nMean centered, variance scaled, normalized expression"
#       , "\n")
#   )
#   # Not centered scaled
#   p2 <- DE_Mean_Heatmap(clusterDeDF = deDF
#     , exDF = noCentExM
#     , clusterIDs = so@ident
#     , upLim = 3
#     , lowLim = 0
#     , ggtitle = paste0(
#       "\nNormalized expression"
#       , "\n")
#   )
#   # plot_grid
#   pg <- plot_grid(p1, p2, ncol = 2)
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nSignificant (p-value < 0.05) DE genes for cluster ", clusterID
#     , "\nTop 40 DE genes"
#     , "\nMean expression"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
#   ggsave(paste0(outGraph, "ExprMeanHeatmap_Cluster", clusterID, ".png")
#     , width = 12, height = 8)
# })
# ################################################################################
