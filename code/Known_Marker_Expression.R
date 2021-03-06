# Damon Polioudakis
# 2017-02-21
# Plot known marker lists as heatmaps, violin plots, feature plots

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
source("Function_Library.R")

# Set variable to gene of interest

## Inputs

# Seurat
# PC 1-40
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# Remove CC genes from variable gene list used for clustering
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# # Regress CC score
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccTirosh <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")

# Lake 2017 cluster DE tables
lake_DE_DF <- read.csv(
  "../source/Lake_2018_TS3_Cluster_DE_Genes.csv"
  , skip = 4, header = TRUE, fill = TRUE
)

# Nowakowski cluster DE table
nowakowski_DE_DF <- read.csv(
  "../nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)

# # Molyneaux markers
# mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)
#
# # Markers from Lake 2016
# lmDF <- read.csv("../source/Lake_2016_ST5.csv", skip = 4, header = TRUE
#   , fill = TRUE)

# Luis metaMat results
mmapDF <- read.csv("../source/Gene_Lists/Overlapped-Genes.csv", header = TRUE)

## Variables
graphCodeTitle <- "Known_Marker_Expression.R"
outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Known_Marker_Expression_"
# outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrainCC_PC1to40/Known_Marker_Expression_"
# outGraph <- "../analysis/graphs/Known_Marker_Expression/DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_RemCC_PC1to40/Known_Marker_Expression_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

Format_Cell_Cycle_Genes_Table <- function(){
  print("Format_Cell_Cycle_Genes_Table")
  ccDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(ccDF) <- c("Grouping", "Gene.Symbol")
  ccDF$Gene.Symbol <- gsub(" *", "", ccDF$Gene.Symbol)
  ccDF <- ccDF[! ccDF$Gene.Symbol == "", ]
  ccDF <- ccDF[! is.na(ccDF$Grouping), ]
  ccDF$Grouping <- gsub(" *", "", ccDF$Grouping)
  ccDF$Grouping <- factor(ccDF$Grouping, levels = unique(ccDF$Grouping))
  return(ccDF)
}
################################################################################

### Format

## Luis markers - cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))

## Lake
lake_DE_DF <- lake_DE_DF[ ,1:7]
# DE filter higher
lake_DE_DF <- lake_DE_DF[
  lake_DE_DF$"Average.Difference..log.fold.change." > 0.25
  , c("Gene", "Average.Difference..log.fold.change.", "Cluster")]

## Nowakowski
# DE filter higher
nowakowski_DE_DF <- nowakowski_DE_DF[
  nowakowski_DE_DF$avg_diff > 0.75, c("gene", "avg_diff", "cluster")]
################################################################################

### Plots for paper

Plot_Paper_Figures <- function(){

  print("Plot_Paper_Figures")

  ## Subset of marker genes heatmap
  genesL <- list(
    vRG = c("CRYAB")
    , oRG = c("HOPX")
    , RG = c("VIM", "PAX6", "HES1", "SOX2")
    , IPC = c("EOMES")
    , Neuron = c("STMN2", "NEUROD6")
    , "Excitatory deep layer" = c("BCL11B", "TBR1", "SOX5")
    , "Excitatory upper layer" = c("SATB2")
    , Interneuron = c("DLX1", "DLX2", "SST", "CALB2")
    , OPC = c("OLIG1", "OLIG2")
    , Endothelial = c("ITM2A", "CLDN5")
    , Pericyte = c("RGS5")
    , Microglia = c("AIF1", "CX3CR1")
  )
  df1 <- melt(genesL)
  # Expression heatmap, cells ordered by cluster
  # Normalized, mean centered and scaled
  geneGroupDF <- data.frame(Gene = df1$value, Group = df1$L1)
  cellID_clusterID <- centSO@ident
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
  )
  ggsave(paste0(outGraph
      , "KnownMarksSubset_HeatmapSetColWidths_NormalizedCenteredScaled_paper.png"
    )
    , width = 10, height = 7
  )


  ## Subset of marker genes feature plots  - individual genes

  # Genes
  genesL <- list(
    vRG = c("CRYAB")
    , oRG = c("HOPX")
    , RG = c("VIM", "PAX6", "HES1", "SOX2")
    , IPC = c("EOMES")
    , Neuron = c("STMN2", "NEUROD6")
    , "Excitatory deep layer" = c("BCL11B", "TBR1", "SOX5")
    , "Excitatory upper layer" = c("SATB2")
    , Interneuron = c("DLX1", "DLX2", "SST", "CALB2")
    , OPC = c("OLIG1", "OLIG2")
    , Endothelial = c("ITM2A", "CLDN5")
    , Pericyte = c("RGS5")
    , Microglia = c("AIF1", "CX3CR1")
  )
  gene_group_DF <- melt(genesL)
  # Collect tSNE values for ggplot
  tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Normalized centered scaled
  ggL <- FeaturePlot(
    genes = gene_group_DF$value
    , tsneDF = tsneDF
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = NULL
    , centScale = TRUE
  )
  ggL <- lapply(ggL, function(gg){
    gg <- gg + ggplot_set_theme_publication_nolabels +
      theme(legend.position = "none")
    return(gg)
  })
  Plot_Grid(
    ggPlotsL = ggL[-1], ncol = 6, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nExpression of marker genes"
      , "\nNormalized centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph
      , "KnownMarksSubset_FeaturePlotIndividual_NormalizedCenteredScaled_paper.png"
    )
    , width = 24, height = 16, limitsize = FALSE, dpi = 150
  )

  ## Subset of marker genes feature plots - mean expression

  # Genes
  gene_group_DF <- kmDF[kmDF$Grouping %in% c("vRG", "oRG", "RG", "IP", "Neuron"
    , "Excitatory Deep Layer Cortical", "Excitatory Upper Layer Cortical"
    , "GABAergic interneuron"), ]
  # Collect tSNE values for ggplot
  tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Normalized centered scaled
  ggL <- FeaturePlot(
    genes = gene_group_DF$Gene.Symbol
    , tsneDF = tsneDF
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = gene_group_DF$Grouping
    , centScale = TRUE
  )
  ggL <- lapply(ggL, function(gg){
    gg <- gg + ggplot_set_theme_publication
    return(gg)
  })
  ggL[[1]] <- ggL[[1]] + theme(legend.position = "none")
  Plot_Grid(
    ggPlotsL = ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nExpression of marker genes"
      , "\nNormalized centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph
      , "KnownMarksSubset_FeaturePlot_NormalizedCenteredScaled_paper.png"
    )
    , width = 12, height = 9, limitsize = FALSE
  )

  ## Cell cycle genes feature plots

  # Genes
  gene_group_DF <- Format_Cell_Cycle_Genes_Table()

  # Collect tSNE values for ggplot
  tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Normalized centered scaled
  ggL <- FeaturePlot(
    genes = gene_group_DF$Gene.Symbol
    , tsneDF = tsneDF
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = gene_group_DF$Grouping
    , centScale = TRUE
  )
  ggL <- lapply(ggL, function(gg){
    gg <- gg + ggplot_set_theme_publication
    return(gg)
  })
  ggL[[1]] <- ggL[[1]] + theme(legend.position = "none")
  Plot_Grid(
    ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nExpression of Macosko 2015 cell cycle marker genes"
      , "\nNormalized centered scaled expression"
      , "\n")
  )
  ggsave(paste0(outGraph
      , "CellCycleSubset_FeaturePlot_NormalizedCenteredScaled_paper.png"
    )
    , width = 12, height = 7, limitsize = FALSE
  )
}
################################################################################

### Luis marker genes

Plot_Luis_Markers <- function(){

  print("Plot_Luis_Markers")

  ## Feature plots

  # Feature plot
  # Normalized, no mean centering scaling
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  ggL <- FeaturePlot(
    genes = kmDF$Gene.Symbol
    , tsneDF = ggDF
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -1
    , limHigh = 3
    , geneGrouping = kmDF$Grouping
    , centScale = FALSE)
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "KnownMarks_FeaturePlot_Normalized.png")
    , width = 12, height = length(ggL)*1+4)

  # Feature plot
  # Normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  ggL <- FeaturePlot(
    genes = kmDF$Gene.Symbol
    , tsneDF = ggDF
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = kmDF$Grouping
    , centScale = TRUE)
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of known marker genes"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "KnownMarks_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1+4)


  ## Violin plots of genes by cluster

  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(outGraph, "KnownMarks_ViolinPlot_NormalizedCenteredScaled.png")
    , width = 13, height = 26)

  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "KnownMarks_ViolinPlot_Normalized.png")
    , width = 13, height = 26)


  ## Heatmaps

  # Expression heatmap, cells ordered by cluster
  # Normalized, mean centered and scaled
  geneGroupDF <- data.frame(Gene = df1$value, Group = df1$L1)
  cellID_clusterID <- centSO@ident
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = geneGroupDF
    , exprM = centSO@scale.data
    , cellID_clusterID = cellID_clusterID
  )
  ggsave(paste0(outGraph
      , "KnownMarksSubset_HeatmapSetColWidths_NormalizedCenteredScaled_paper.png"
    )
    , width = 10, height = 7
  )

  # Heatmap
  # Normalized, no mean centering scaling
  geneGroupDF <- data.frame(Gene = kmDF$Gene.Symbol, Group = kmDF$Grouping)
  ggL <- Heatmaps_By_Cluster_Combined(
    geneGroupDF = geneGroupDF, exprM = noCentExM, seuratO = centSO
    , lowerLimit = 0, upperLimit = 3, geneOrder = FALSE
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:16)
  )
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # now add the title
  title <- paste0(graphCodeTitle
    , "\n\nExpression of published marker genes by cluster"
    , "\nx-axis: Marker genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
    , "\n")
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "KnownMarks_Heatmap_Normalized.png")
    , width = 16, height = 80, limitsize = FALSE)

  # Heatmap
  # Normalized, mean centered and scaled
  geneGroupDF <- data.frame(Gene = kmDF$Gene.Symbol, Group = kmDF$Grouping)
  ggL <- Heatmaps_By_Cluster_Combined(
    geneGroupDF = geneGroupDF, exprM = centSO@scale.data, seuratO = centSO
    , lowerLimit = -1.5, upperLimit = 1.5, geneOrder = FALSE
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:16)
  )
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # now add the title
  title <- paste0(graphCodeTitle
    , "\n\nExpression of published marker genes by cluster"
    , "\nx-axis: Marker genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered, variance scaled"
    , "\n")
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "KnownMarks_Heatmap_NormalizedCenteredScaled.png")
    , width = 16, height = 80, limitsize = FALSE)
}
################################################################################

### Cell cycle genes

Plot_Cell_Cycle_Genes <- function(){

  print("Plot_Cell_Cycle_Genes")

  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

  # Subset to marker genes of interest for Luis' excel file
  # Cleanup marker data frame
  ccDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(ccDF) <- c("Grouping", "Gene.Symbol")
  ccDF$Gene.Symbol <- gsub(" *", "", ccDF$Gene.Symbol)
  ccDF <- ccDF[! ccDF$Gene.Symbol == "", ]
  ccDF <- ccDF[! is.na(ccDF$Grouping), ]
  ccDF$Grouping <- gsub(" *", "", ccDF$Grouping)
  ccDF$Grouping <- factor(ccDF$Grouping, levels = unique(ccDF$Grouping))
  ccDFL <- split(ccDF, ccDF$Grouping)

  # We can segregate this list into markers of G2/M phase and markers of S
  # phase
  ccDFL[["S.Tirosh"]] <- data.frame(
    Grouping = "S.Tirosh", Gene.Symbol = ccTirosh[1:43])
  ccDFL[["G2M.Trisosh"]] <- data.frame(
    Grouping = "G2M.Tirosh", Gene.Symbol = ccTirosh[44:98])

  # Feature plot
  # Normalized, no mean centering scaling
  # Loop through and plot each group of genes
  ggL <- lapply(names(ccDFL), function(grouping) {
    print(grouping)
    genes <- ccDF$Gene.Symbol[ccDF$Grouping == grouping]
    ggDF <- Mean_Expression(ggDF, genes, noCentExM)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 3)
    ggFp <- Feature_Plot(ggDF, title = paste0("\n", grouping)
      , limLow = 0, limHigh = 3)
    return(ggFp)
  })
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  ggTsne <- TSNE_Plot(centSO)
  ggL <- append(list(ggTsne), ggL)
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of cell cycle genes from Macosko"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
  ggsave(paste0(outGraph, "CellCycle_FeaturePlot_Normalized.png")
    , width = 14, height = length(ggL)*1.5)

  # Feature plot
  # Normalized, mean centered scaled
  # Loop through and plot each group of genes
  ggL <- lapply(names(ccDFL), function(grouping) {
    print(grouping)
    genes <- ccDF$Gene.Symbol[ccDF$Grouping == grouping]
    ggDF <- Mean_Expression(ggDF, genes, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = -1.5, limHigh = 1.5)
    ggFp <- Feature_Plot_CentScale(ggDF, title = paste0("\n", grouping)
      , limLow = -1.5, limHigh = 1.5)
    return(ggFp)
  })
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  ggTsne <- TSNE_Plot(centSO)
  ggL <- append(list(ggTsne), ggL)
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of cell cycle genes from Macosko"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
  ggsave(paste0(outGraph, "CellCycle_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 14, height = length(ggL)*1.5)

  # Heatmap
  # Normalized, no mean centering scaling
  ggDF <- merge(ccDF[c("Gene.Symbol", "Grouping")], noCentExM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  Heatmaps_By_Cluster_Combined(
    ggDF = ggDF, seuratO = centSO, lowerLimit = 0, upperLimit = 3
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
    , title = paste0(graphCodeTitle
      , "\n\nExpression of cell cycle genes (Macosko) by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "CellCycle_Heatmap_Normalized.png")
    , width = 16, height = 44)

  # Normalized, mean centered and scaled
  ggDF <- merge(ccDF[c("Gene.Symbol", "Grouping")], centSO@scale.data
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  Heatmaps_By_Cluster_Combined(
    ggDF = ggDF, seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
    , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
    , title = paste0(graphCodeTitle
      , "\n\nExpression of cell cycle genes (Macosko) by cluster"
      , "\nx-axis: Marker genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
  )
  ggsave(paste0(outGraph, "CellCycle_Heatmap_NormalizedCenteredScaled.png")
    , width = 16, height = 44)
}
################################################################################

### Molyneaux layer marker genes

Plot_Molyneaux_Layer_Markers <- function(){

  print("Plot_Molyneaux_Layer_Markers")

  ## Feature plot

  # Collect tSNE values for ggplot
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

  ## Feature plot - normalized, mean centered scaled on full dataset
  # Loop through and plot each group of genes
  ggL <- apply(mmDF, 1, function(v1) {
    gene <- v1[["hgnc_symbol"]]
    print(gene)
    ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = -1.5, limHigh = 1.5)
    ggFp <- Feature_Plot_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
      , title = paste0("\n", gene, " - ", v1[["mgi_symbol"]])
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(centSO)
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
    , "\n\nExpression of layer markers from Molyneaux review"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n\nRemove cells < 200 genes detected"
    , "\nRemove cells > 3192 (3 SD) genes detected"
    , "\nRemove genes detected in < 3 cells"
    , "\nNormalize expression"
    , "\nRegress out nUMI, donor, library lab"
    , "\ntSNE PC 1-40"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "MolyneauxMarks_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 14, height = length(ggL), limitsize = FALSE)

  ## Feature plot - normalized, not mean centered scaled on cluster
  # Loop through and plot each group of genes
  ggL <- apply(df, 1, function(v1) {
    gene <- v1[["Gene.Symbol"]]
    print(gene)
    exDF <- noCentExM[row.names(noCentExM) %in% row.names(centSO@scale.data)
      , colnames(noCentExM) %in% colnames(centSO@scale.data)]
    ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
    ggDF <- Set_Limits(ggDF, limLow = 0, limHigh = 2)
    ggFp <- Feature_Plot(ggDF
      , title = paste0("\n", gene, " - ", v1[["mgi_symbol"]])
    )
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(centSO)
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
    , "\n\nExpression of layer markers from Molyneaux review"
    , "\nNormalized expression"
    , "\n\nRemove cells < 200 genes detected"
    , "\nRemove cells > 3192 (3 SD) genes detected"
    , "\nRemove genes detected in < 3 cells"
    , "\nNormalize expression"
    , "\nRegress out nUMI, donor, library lab"
    , "\ntSNE PC 1-40"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1
    , rel_heights = c(length(ggL)*0.05, length(ggL)))
  ggsave(paste0(
    outGraph, "PC1-40_Molyneaux_FeaturePlot_NormClusterCentScale_Cluster", cl, ".png")
    , width = 14, height = length(ggL), limitsize = FALSE)

  ## Violin plots of genes by cluster

  Gene_Expression_By_Cluster_ViolinPlot <- function(genes, exprM, seuratO, clusterIDs
    , geneOrder = NULL){
    # Normalized, no mean centering scaling
    ggDF <- merge(genes, exprM
      , by.x = 1, by.y = "row.names", all.x = TRUE)
    row.names(ggDF) <- ggDF$x
    ggDF <- ggDF[ ,-1]
    ggDF <- t(ggDF)
    ggDF <- as.data.frame(ggDF)
    # Add cluster ID
    idx <- match(names(centSO@ident), row.names(ggDF))
    print(head(idx, 20))
    ggDF$CLUSTER <- centSO@ident[idx]
    ggDF <- melt(ggDF)
    # Set gene order
    ggDF$variable <- factor(ggDF$variable, levels = genes)
    print(str(ggDF))
    # Violin plots of genes by cluster
    png(paste0(outGraph, "ClusterViolinPlot.png")
      , width = 13, height = 6, units = "in", res = 300)
    ggplot(ggDF, aes(x = CLUSTER, y = value)) +
      geom_violin(aes(fill = CLUSTER)) +
      geom_jitter(size = 0.05, height = 0, alpha = 0.1) +
      facet_wrap(~variable, scales = "free", ncol = 5) +
      theme(legend.position = "none") +
      ylab("Normalized expression") +
      xlab("Clusters")
  }

  genes <- as.character(mmDF$hgnc_symbol)
  genes <- factor(genes, levels = genes)
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = genes
    , exprM = noCentExM
    , seuratO = centSO
    , clusterIDs = centSO@ident) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression - read depth, regressed nUMI, brain, lablibrary"
      , "\n"))
  ggsave(paste0(outGraph, "PC1-40_Molyneaux_ViolinPlot_Norm.png"), height = 30)


  ScatterPlot_Expression <- function (genes, exprM, limLow, limHigh){
    # Normalized, no mean centering scaling
    ggDF <- merge(genes, exprM
      , by.x = 1, by.y = "row.names", all.x = TRUE)
    row.names(ggDF) <- ggDF$x
    ggDF <- ggDF[ ,-1]
    ggDF <- t(ggDF)
    ggDF <- as.data.frame(ggDF)
    # Order by genes
    ggDF <- ggDF[ ,match(genes, colnames(ggDF))]
    # # Limits
    # ggDF[ggDF > limHigh] <- limHigh
    # ggDF[ggDF < limLow] <- limLow
    print(head(ggDF))
    # PAX6 STMN2 EOMES
    gg <- ggplot(ggDF, aes_string(x = colnames(ggDF)[1], y = colnames(ggDF)[2])) +
      geom_point(size = 0.5, alpha = 0.5, aes_string(color = colnames(ggDF)[3])) +
      # scale_color_distiller(name = "Normalized\nexpression", type = "div"
      #   , palette = 5, direction = -1) +
      scale_color_viridis()
    # # Limits
    # ggDF[ggDF > limHigh] <- limHigh
    # ggDF[ggDF < limLow] <- limLow
    # print(head(ggDF))
    # # PAX6 STMN2 EOMES
    # p1 <- ggplot(ggDF, aes_string(x = colnames(ggDF)[1], y = colnames(ggDF)[2])) +
    #   geom_point(size = 0.5, alpha = 0.5, aes_string(color = colnames(ggDF)[3])) +
    #   coord_cartesian(xlim = c(limLow, limHigh), ylim = c(limLow, limHigh)) +
    #   # scale_color_distiller(name = "Normalized\nexpression", type = "div"
    #   #   , palette = 5, direction = -1) +
    #   scale_color_viridis(limits = c(limLow, limHigh))
    return(gg)
  }

  ScatterPlot_Expression(c("TLE1", "TLE3", "MDGA1"), exprM = noCentExM)
  ggsave(paste0(outGraph, "PC1-40_Molyneaux_Scatter_TLE1_TLE3_Norm.png")
    , width = 7, height = 7)

  ScatterPlot_Expression(c("BCL11B", "SOX5", "TLE4"), exprM = noCentExM)
  ggsave(paste0(outGraph, "PC1-40_Molyneaux_Scatter_BCL11B_SOX5_Norm.png")
    , width = 7, height = 7)

  ScatterPlot_Expression(c("TLE1", "BCL11B", "TLE4"), exprM = noCentExM)
  ggsave(paste0(outGraph, "PC1-40_Molyneaux_Scatter_TLE1_BCL11B_Norm.png")
    , width = 7, height = 7)

  ScatterPlot_Expression(c("TLE3", "BCL11B", "TLE4"), exprM = noCentExM)
  ggsave(paste0(outGraph, "PC1-40_Molyneaux_Scatter_TLE3_BCL11B_Norm.png")
    , width = 7, height = 7)
}
################################################################################

### Lake cluster enriched

Plot_Lake2017_Enriched <- function(){

  print("Plot_Lake2017_Enriched")

  ## Feature plots

  # Feature plot
  # Normalized, no mean centering scaling
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = lake_DE_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = lake_DE_DF$Cluster
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Lake 2017 cluster enriched"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "LakeEnriched_FeaturePlot_Normalized.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = lake_DE_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = as.matrix(centSO@scale.data)
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = lake_DE_DF$Cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Lake 2017 cluster enriched"
    , "\nCluster round 1 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "LakeEnriched_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)

}
################################################################################

### Nowakowski cluster enriched

Plot_Nowakowski_Enriched <- function(){

  print("Plot_Nowakowski_Enriched")

  # Feature plot
  # Normalized, no mean centering scaling
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = nowakowski_DE_DF$gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = nowakowski_DE_DF$cluster
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Nowakowski cluster enriched"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph, "NowakowskiEnriched_FeaturePlot_Normalized.png")
    , width = 12, height = length(ggL)*1)

  # Feature plot
  # Normalized, mean centered scaled
  # Collect tSNE values for ggplot
  ggL <- FeaturePlot(
    genes = nowakowski_DE_DF$gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = nowakowski_DE_DF$cluster
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(graphCodeTitle
    , "\n\nExpression of Nowakowski cluster enriched"
    , "\nNormalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(paste0(outGraph
    , "NowakowskiEnriched_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = length(ggL)*1)

}
################################################################################

### Luis layer marker overlaps

Plot_Marker_Intersections <- function(){

  print("Plot_Marker_Intersections")

  # Violin plots of expression
  # Normalized
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% c(
    "Excitatory Upper Layer Cortical", "Excitatory Deep Layer Cortical")]
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = genes
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , geneOrder = genes
  ) +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(outGraph, "LuisLayerMarks_ViolinPlot_Normalized.png")
    , width = 13, height = 20)


  ## Number of pairwise intersections of upper and lower layer markers

  uplGenes <- kmDF$Gene.Symbol[kmDF$Grouping %in% c(
    "Excitatory Upper Layer Cortical")]
  lowlGenes <- kmDF$Gene.Symbol[kmDF$Grouping %in% c(
    "Excitatory Deep Layer Cortical")]
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% c(
    "Excitatory Upper Layer Cortical", "Excitatory Deep Layer Cortical")]

  m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
  l1 <- apply(m1, 2, function(col) row.names(m1)[col])

  ldf <- lapply(l1, function(genes) {
    v2 <- genes[genes %in% uplGenes]
    v3 <- genes[genes %in% lowlGenes]
    expand.grid(v2, v3)
  })
  df1 <- do.call("rbind", ldf)

  # Format for barplot
  df2 <- data.frame(table(paste(df1[ ,1], df1[ ,2])))

  # Format for matrix
  df3 <- dcast(df1, Var1 ~ Var2)
  row.names(df3) <- df3$Var1
  df3 <- df3[ ,-1]
  v1 <- rowSums(m1)
  # Union of each gene to normalize for frequency
  m2 <- matrix(NA, nrow(df3), ncol(df3))
  colnames(m2) <- colnames(df3)
  rownames(m2) <- rownames(df3)
  for (cName in colnames(m2)) {
    for (rName in rownames(m2)) {
      m2[rName,cName] <- v1[[cName]] + v1[[rName]]
    }
  }
  # Normalize for frequency
  df4 <- round((df3 / m2) * 100, 1)
  df4$Gene <- row.names(df4)
  df4 <- melt(df4)

  # Format for ggplot
  df3$Gene <- row.names(df3)
  df3 <- melt(df3)

  # number of cells each marker expressed in
  df5 <- data.frame(Number = rowSums(m1))
  df5$Gene <- factor(row.names(df5), levels = row.names(df5))

  # Plot number of cells each marker expressed in
  ggplot(df5, aes(x = Gene, y = Number)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nNumber of cells expressing layer markers"))
  ggsave(paste0(outGraph, "LuisLayerMarks_Number_Barplot.png")
    , width = 5, height = 5)

  # Plot counts of intersections
  ggplot(df3, aes(x = Gene, y = variable, fill = value)) +
    geom_tile() +
    geom_text(data = df3, aes(x = Gene, y = variable, label = value), color = "black") +
    scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "Number") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Upper layer markers") +
    ylab("Lower layer markers") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nNumber of cells expressing both upper and lower layer markers"))
  ggsave(paste0(outGraph, "LuisLayerMarks_NumberIntersect_Heatmap.png")
    , width = 7, height = 7)

  # Plot counts of intersections normalized for frequency
  ggplot(df4, aes(x = Gene, y = variable, fill = value)) +
    geom_tile() +
    geom_text(data = df4, aes(x = Gene, y = variable, label = value), color = "black") +
    scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "Percent") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Upper layer markers") +
    ylab("Lower layer markers") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nCells expressing both upper and lower layer markers"
      , "\n(Intersection / Union) * 100")
      , "\nExample: 24.3% of cells expressing either POU3F2 or SOX5 express both")
  ggsave(paste0(outGraph, "LuisLayerMarks_PercentIntersect_Heatmap.png")
    , width = 7, height = 7)

  # Barplot of numbers of intersections
  ggplot(df2, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ylab("Combinations of upper and lower layer markers") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nNumber of occurances of upper and lower layer markers expressed in same cell"
      , "\nExpression filter: > 0.5 log normalized expression"))
  ggsave(paste0(outGraph, "LuisLayerMarks_NumberIntersect.png")
    , width = 7, height = 7)
}
################################################################################

### Run

Main_Function(){
  Plot_Paper_Figures()
  Plot_Luis_Markers()
  Plot_Cell_Cycle_Genes()
  # Plot_Molyneaux_Layer_Markers()
  # Plot_Marker_Intersections()
  Plot_Lake2017_Enriched()
  Plot_Nowakowski_Enriched()
  # Plot_Marker_Intersections()
}

Main_Function()
################################################################################
