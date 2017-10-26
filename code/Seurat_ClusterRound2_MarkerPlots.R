# Damon Polioudakis
# 2017-05-28
# 2nd iteration of Seurat clustering

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
require(irlba)
require(gridExtra)
require(cowplot)
source("Function_Library.R")

## Inputs

# Seurat clustering object
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")
# Clear parts of object to save memory
centSO@raw.data <- NULL
centSO@data <- NULL

# Seurat clustering round 2 object
# load("../analysis/Seurat_ClusterRound2_DS2-11/Seurat_ClusterRound2_DS2-11_VarGenes_PC1-30_seuratO.Robj")
# load("../analysis/Seurat_ClusterRound2_DS2-11/Seurat_ClusterRound2_DS2-11_VarGenes_RegNumiLibBrain_PC1-30_seuratO.Robj")
load("../analysis/Seurat_ClusterRound2_DS2-11/Seurat_ClusterRound2_DS2-11_AllGenes_RegNumiLibBrain_PC1-30_seuratO.Robj")

# Cell cycle markers from Macosko 2015 Table S2 to remove from variable gene
# list used for clustering
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# Molyneaux markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

# Markers from Lake 2016
lmDF <- read.csv("../source/Lake_2016_ST5.csv", skip = 4, header = TRUE, fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_ClusterRound2_MarkerPlots.R"
# outGraph <- "../analysis/graphs/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_VarGenes_RegNumiLibBrain_PC1-30_"
# outTable <- "../analysis/tables/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_VarGenes_RegNumiLibBrain_PC1-30_"
# outData <- "../analysis/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_VarGenes_RegNumiLibBrain_PC1-30_"
outGraph <- "../analysis/graphs/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_AllGenes_RegNumiLibBrain_PC1-30_"
outTable <- "../analysis/tables/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_AllGenes_RegNumiLibBrain_PC1-30_"
outData <- "../analysis/Seurat_ClusterRound2_MarkerPlots_DS2-11/Seurat_ClusterRound2_MarkerPlots_DS2-11_AllGenes_RegNumiLibBrain_PC1-30_"

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

print("### Feature plot of covariates")

# Loop through clusters
pgL <- lapply(names(lso), function(cl) {
  
  so <- lso[[cl]]
  
  ## Plot tSNE graph colored by lanes, samples, or GZ CP
  # Collect tSNE values
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  # Add cluster identity
  ggDF$CLUSTER <- so@ident
  # Add metadata
  ggDF <- data.frame(ggDF, so@meta.data[match(row.names(ggDF), so@meta.data$CELL), ])
  ggDF$BRAIN <- as.factor(ggDF$BRAIN)
  
  # ggplot Donor
  gg1 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = BRAIN)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggtitle(paste0("Donor"
      , "\nCluster: ", cl))
  
  # ggplot GZ/CP
  gg2 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = REGION)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggtitle(paste0("GZ/CP"
      , "\nCluster: ", cl))
  
  # ggplot Lab library
  gg3 <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = LIBRARY)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggtitle(paste0("Lab library"
      , "\nCluster: ", cl))
  
  # tSNE colored by clustering
  ggTsne <- TSNE_Plot(so) + theme(legend.position = "none")
  
  # plot grid
  pg <- plot_grid(ggTsne, gg1, gg2, gg3, align = "v", axis = "l", ncol = 4)
  return(pg)
})
pg <- plot_grid(plotlist = pgL, ncol = 1)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nSeurat clustering and covariates"
  , "\n"
  , "\nRemove genes > 0 counts in < 3 cells"
  , "\nRemove cells < 200 genes detected"
  , "\nRemove cells > 3192 (3 SD) genes detected"
  , "\nRemove genes detected in < 3 cells"
  , "\nNormalize expression"
  , "\nRegress out covariates"
  , "\ntSNE + Seurat cluster PC 1-40, cluster round 2 tSNE PCA 1-30"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
# save
ggsave(paste0(outGraph, "tSNE_covariates.png")
  , width = 20, height = 80, limitsize = FALSE)
################################################################################

print("### Feature plot of top PC scores")

# Loop through clusters
lapply(names(lso), function(cl) {
  
  so <- lso[[cl]]
  
  # Feature plots of top PC scores
  ggL <- lapply(c(1:8), function(i) {
    # Collect tSNE values
    ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
    # Add PC score
    ggDF$PC <- so@dr$pca@cell.embeddings[ ,i]
    gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PC)) +
      geom_point(size = 0.05, alpha = 0.5) +
      # guides(colour = guide_legend(override.aes = list(size = 7))) +
      scale_color_distiller(name = "PC score", type = "div"
        , palette = 5, direction = -1) +
      ggtitle(paste0("PC: ", i, "\n")) +
      theme(text = element_text(size = 14))
    return(gg)
  })
  # Add tSNE colored by clustering for reference
  ggTsne <- TSNE_Plot(so) + theme(legend.position = "none")
  ggL <- append(list(ggTsne), ggL)
  # plot grid
  pg <- plot_grid(plotlist = ggL, align = "v", axis = "r", ncol = 3)
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nSeurat cluster and tSNE colored by PC scores for cluster: ", cl
    , "\n"
    , "\nRemove genes > 0 counts in < 3 cells"
    , "\nRemove cells < 200 genes detected"
    , "\nRemove cells > 3192 (3 SD) genes detected"
    , "\nRemove genes detected in < 3 cells"
    , "\nNormalize expression"
    , "\nRegress out covariates"
    , "\ntSNE + Seurat cluster PC 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
  # save
  ggsave(paste0(outGraph, "PCscore_FeaturePlot_Cluster", cl, ".png")
    , width = 15, height = 15)
})
################################################################################

# ### Feature plot genes with high PC loadings
# 
# # Genes with highest PC loadings
# ldf <- lapply(1:4, function(pc) {
#   data.frame(GENES = names(sort(abs(so@dr$pca@gene.loadings[ ,pc])
#     , decreasing = TRUE)[1:5]), PC = paste0("PC", pc))
# })
# genesDF <- do.call("rbind", ldf)
# 
# # Loop through clusters
# lapply(names(lso), function(cl) {
#   
#   so <- lso[[cl]]
#   
#   # Collect tSNE values for ggplot
#   ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
#   
#   ## Feature plot - normalized, mean centered scaled on full dataset
#   # Loop through and plot each group of genes
#   ggL <- apply(genesDF, 1, function(v1) {
#     gene <- v1[["GENES"]]
#     grouping <- v1[["PC"]]
#     print(gene)
#     ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
#     ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
#     ggFp <- Feature_Plot(ggDF, limLow = 0, limHigh = 2
#       , title = paste0(gene, "\n", grouping)
#     )
#     return(ggFp)
#   })
#   ggTsne <- TSNE_Plot(so)
#   ggL <- append(list(ggTsne), ggL)
#   # extract the legend from one of the plots
#   legend <- get_legend(ggL[[2]])
#   # Remove legends from plots
#   ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
#   # plot_grid combine tSNE graphs
#   pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
#   # add the legend to the row we made earlier. Give it one-third of the width
#   # of one plot (via rel_widths).
#   pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
#   # now add the title
#   title <- ggdraw() + draw_label(paste0(graphCodeTitle
#     , "\n\nGenes (5) that correlate highly to top PCs (4)"
#     , "\ntSNE plot of reclustered cells, each point is a cell"
#     , "\nColor indicates mean normalized centered scaled expression"
#     , "\n"
#     , "\nNormalized, mean centered, variance scaled across full dataset"
#     , "\nRemove cells < 200 genes detected in cluster"
#     , "\nRemove genes detected in < 3 cells in cluster"
#     , "\nOnly recluster clusters >= 100 cells"
#     , "\nSeurat variable genes used for clustering"
#     , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
#     , "\n"))
#   # rel_heights values control title margins
#   plot_grid(title, pg, ncol = 1
#     , rel_heights = c(length(ggL)*0.2, length(ggL)))
#   ggsave(paste0(
#     outGraph, "PChighGenes_FeaturePlot_NormFullCentScale_Cluster", cl, ".png")
#     , width = 16, height = length(ggL), limitsize = FALSE)
#   
#   # ## Feature plot - normalized, mean centered scaled on cluster
#   # # Loop through and plot each group of genes
#   # ggL <- apply(df, 1, function(v1) {
#   #   gene <- v1[["Gene.Symbol"]]
#   #   grouping <- v1[["Grouping"]]
#   #   print(gene)
#   #   ggDF <- Mean_Expression(gene, so@scale.data)
#   #   ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
#   #   ggFp <- Feature_Plot(ggDF
#   #     , title = paste0(gene, "\n", grouping)
#   #   )
#   #   return(ggFp)
#   # })
#   # ggTsne <- TSNE_Plot(so)
#   # ggL <- append(list(ggTsne), ggL)
#   # # extract the legend from one of the plots
#   # legend <- get_legend(ggL[[2]])
#   # # Remove legends from plots
#   # ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
#   # # plot_grid combine tSNE graphs
#   # pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
#   # # add the legend to the row we made earlier. Give it one-third of the width
#   # # of one plot (via rel_widths).
#   # pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
#   # # now add the title
#   # title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   #   , "\n\nGenes (5) that correlate highly to top PCs (4)"
#   #   , "\ntSNE plot of reclustered cells, each point is a cell"
#   #   , "\nColor indicates mean normalized centered scaled expression"
#   #   , "\n"
#   #   , "\nNormalized, mean centered, variance scaled across cluster"
#   #   , "\nRemove cells < 200 genes detected in cluster"
#   #   , "\nRemove genes detected in < 3 cells in cluster"
#   #   , "\nOnly recluster clusters >= 100 cells"
#   #   , "\nSeurat variable genes used for clustering"
#   #   , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
#   #   , "\n"))
#   # # rel_heights values control title margins
#   # plot_grid(title, pg, ncol = 1
#   #   , rel_heights = c(length(ggL)*0.2, length(ggL)))
#   # ggsave(paste0(
#   #   outGraph, "KnownMarks_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
#   #   , width = 16, height = length(ggL), limitsize = FALSE)
#   
# })
################################################################################

### Feature plot Luis known markers

# Mean expression of marker gene groups for each cluster2
lapply(names(lso), function(cl) {
  
  print(cl)
  so <- lso[[cl]]
  print(head(so@dr$tsne@cell.embeddings))
  
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  ## Luis marker genes
  
  # Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
  kmDF <- kmDF[! is.na(kmDF$Grouping), ]
  kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
  kmDFL <- split(kmDF, kmDF$Grouping)
  
  # Feature plot
  # Normalized, mean centered scaled
  # Loop through and plot each group of genes
  ggL <- lapply(names(kmDFL), function(grouping) {
    print(grouping)
    genes <- kmDF$Gene.Symbol[kmDF$Grouping == grouping]
    print("Mean expression")
    ggDF <- Mean_Expression(ggDF, genes, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = -1.5, limHigh = 1.5)
    print("Feature plot")
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0("\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized centered scaled expression"
    , "\n"
    , "\nNormalized, mean centered, variance scaled across full dataset"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.1, length(ggL)))
  
  ggsave(paste0(
    outGraph, "KnownMarksGroups_FeaturePlot_NormCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
})

## Marker genes for each cluster2 of cluster 0-4
lapply(names(lso)[1:10], function(cl) {
  
  print(cl)
  so <- lso[[cl]]
  print(head(so@dr$tsne@cell.embeddings))
  
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  ## Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  df <- kmDF[kmDF$Grouping %in% c("IP", "Neuron", "GABAergic", "Glutamatergic"
    , "Dopaminergic Neuron", "Cortical Migratin", "Excitatory Deep Layer Cortical"
    , "Excitatory Upper Layer Cortical", "Subplate"), ]
  ldf <- split(df, df$Gene.Symbol)
  
  ## Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across full dataset"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormFullCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
  ## Feature plot - normalized, mean centered scaled on cluster
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, so@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across cluster"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
})

## Marker genes for each cluster2 of cluster 5-6
lapply(names(lso)[c(4, 11, 12)], function(cl) {
  
  print(cl)
  so <- lso[[cl]]
  print(head(so@dr$tsne@cell.embeddings))
  
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  ## Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  df <- kmDF[kmDF$Grouping %in% c("GABAergic interneuron", "GABAergic Neuron"
    , "Neuron", "Cajal-Retzius", "Cortical Migrating"), ]
  ldf <- split(df, df$Gene.Symbol)
  
  ## Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across full dataset"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormFullCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
  ## Feature plot - normalized, mean centered scaled on cluster
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, so@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across cluster"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
})

## Marker genes for each cluster2 of cluster 2, 7, 8, 9, 10
lapply(names(lso)[c(5, 13, 14, 15, 16)], function(cl) {
  
  print(cl)
  so <- lso[[cl]]
  print(head(so@dr$tsne@cell.embeddings))
  
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  ## Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  df <- kmDF[kmDF$Grouping %in% c("vRG", "oRG"
    , "RG", "IP", "Glia", "Neuron"), ]
  ldf <- split(df, df$Gene.Symbol)
  
  ## Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across full dataset"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormFullCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
  ## Feature plot - normalized, mean centered scaled on cluster
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    grouping <- v1[["Grouping"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, so@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0(gene, "\n", grouping)
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nPublished marker genes, cluster round2 for cluster: ", cl
    , "\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean normalized expression"
    , "\nNormalized, mean centered, variance scaled across cluster"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "KnownMarks_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
})
################################################################################

### Feature plot Molyneaux

# Loop through cluster Seurat objects
lapply(names(lso), function(cl) {
  
  so <- lso[[cl]]
  
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
  
  # ## Subset to marker genes of interest for Luis' excel file
  # # Cleanup marker data frame
  # df <- kmDF[kmDF$Grouping %in% c("IP", "Neuron", "GABAergic", "Glutamatergic"
  #   , "Dopaminergic Neuron", "Cortical Migratin", "Excitatory Deep Layer Cortical"
  #   , "Excitatory Upper Layer Cortical", "Subplate"), ]
  # ldf <- split(df, df$Gene.Symbol)
  
  ## Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- apply(mmDF, 1, function(v1) {
    gene <- v1[["hgnc_symbol"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, so@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = -1.5, limHigh = 1.5)
    ggFp <- FeaturePlot_Graph_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0("\n", gene, " - ", v1[["mgi_symbol"]])
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nMolyneaux layer genes, cluster round2 for cluster: ", cl
    , "\n\ntSNE plot of reclustered cells, each point is a cell"
    , "\nColor indicates mean of normalized expression, mean centered, variance scaled across full dataset"
    , "\n"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "MolyneauxMarks_FeaturePlot_NormFullCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
  ## Feature plot - normalized, not mean centered scaled on cluster
  # Loop through and plot each group of genes
  ggL <- apply(mmDF, 1, function(v1) {
    gene <- v1[["hgnc_symbol"]]
    print(gene)
    exDF <- noCentExM[row.names(noCentExM) %in% row.names(so@scale.data)
      , colnames(noCentExM) %in% colnames(so@scale.data)]
    ggDF <- Mean_Expression(ggDF, gene, exDF)
    ggDF <- Set_Limits(ggDF, limLow = -1, limHigh = 3)
    ggFp <- FeaturePlot_Graph(ggDF, limLow = -1, limHigh = 3
      , title = paste0("\n", gene, " - ", v1[["mgi_symbol"]])
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(so)
  ggL <- append(list(ggTsne), ggL)
  # extract the legend from one of the plots
  legend <- get_legend(ggL[[2]])
  # Remove legends from plots
  ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'v', axis = 'l')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nMolyneaux layer genes, cluster round2 for cluster: ", cl
    , "\nColor indicates mean normalized expression"
    , "\n"
    , "\nRemove cells < 200 genes detected in cluster"
    , "\nRemove genes detected in < 3 cells in cluster"
    , "\nOnly recluster clusters >= 100 cells"
    , "\nSeurat variable genes used for clustering"
    , "\ntSNE PCA 1-40, cluster round 2 tSNE PCA 1-30"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "Molyneaux_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)
  
})
################################################################################