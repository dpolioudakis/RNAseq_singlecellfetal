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
# require(pheatmap)
require(RColorBrewer)
require(viridis)
# require(fdrtool)
require(monocle)
require(WGCNA)
require(gridExtra)
require(ggpubr)
require(ggbeeswarm)
source("Function_Library.R")
source("Seurat_Cluster_Cycling_vRGoRG_Functions.R")

options(stringsAsFactors = FALSE)

## Inputs

# Seurat object
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TESTcluster0278910_seuratO.Robj")
# centSO <- ssCentSO; rm(ssCentSO)
# noCentExM <- ssNoCentExM; rm(ssNoCentExM)

# Cluster DE genes
cluster_DE_DF <- read.table(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , header = TRUE)

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccGenes <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
sGenes <- ccGenes[1:43]
g2mGenes <- ccGenes[44:98]

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_Cluster_Cycling_vRGoRG_DS2-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Load analyzed data

# DE cell types
deDF <- read.csv(paste0(outTable, "DE_CellTypes.csv"))

# DE mixed marker cells
transition_state_DE_DF <- read.csv(
  paste0(outTable, "DE_CellTransitionStates.csv"))

# Eigengenes of DE genes
# ME_CellType_L
load(file = paste0(outData, "ME_CellType.rdata"))

# Eigengenes of cell type enriched genes
# ME_CellTypeEnriched_L
load(file = paste0(outData, "ME_CellTypeEnriched.rdata"))
################################################################################

### Plot eigengene of DE genes

Calculate_Quantile_Limits <- function(df, value_col, group_col){
# The lower and upper hinges correspond to the first and third quartiles (the
# 25th and 75th percentiles). The upper whisker extends from the hinge to the
# largest value no further than 1.5 * IQR from the hinge (where IQR is the
# inter-quartile range, or distance between the first and third quartiles). The
# lower whisker extends from the hinge to the smallest value at most 1.5 * IQR
# of the hinge

  quantile_DF <- aggregate(df[ ,value_col]~df[ ,group_col], df
    , function(x){quantile(x, c(0.25, 0.75))}
  )
  # Flatten dataframe
  quantile_DF <- cbind(quantile_DF[,1], quantile_DF[,2])
  quantile_DF <- as.data.frame(quantile_DF)
  quantile_DF[,2] <- as.numeric(quantile_DF[,2])
  quantile_DF[,3] <- as.numeric(quantile_DF[,3])
  iqr <- quantile_DF[,3] - quantile_DF[,2]
  upper_whiskers <- quantile_DF[,3] + (iqr * 1.5)
  lower_whiskers <- quantile_DF[,2] - (iqr * 1.5)

  # Limits
  limits_L <- list(
    Minimum = min(df[ ,value_col][df[ ,value_col] > lower_whiskers])
    , Maximum = max(df[ ,value_col][df[ ,value_col] < upper_whiskers])
  )
  return(limits_L)
}

Plot_ME_CellType_Genes <- function(
  me_markerFlag_DF, color_1 = NA, color_2 = NA){
  # browser()

  print("Plot_ME_CellType_Genes")

  me_markerFlag_DF$variable <- factor(me_markerFlag_DF$variable
    , levels = c("RG eigengene", "IP eigengene", "Neuron eigengene"))
  me_markerFlag_DF$variable <- droplevels(me_markerFlag_DF$variable)

  me_markerFlag_DFL <- split(me_markerFlag_DF
    , list(me_markerFlag_DF$Fold_Change_Cutoff
    , me_markerFlag_DF$variable)
  )

  ggL <- lapply(names(me_markerFlag_DFL), function(name){

    print(name)

    me_markerFlag_DF <- me_markerFlag_DFL[[name]]

    # Fill color of boxplots
    fill_color <- ifelse(grepl("RG", name), "#8dd3c7"
      , ifelse(grepl("IP", name), "#bebada"
        , ifelse(grepl("Neuron", name), "#fb8072", NA
    )))

    # For mixed marker cells subset to cycling cells only
    idx <- me_markerFlag_DF$Cell_Subset_025 %in% c(
        "RG IP", "RG Neuron", "IP Neuron") &
      me_markerFlag_DF$CellID %in% names(centSO@ident)[
        centSO@ident %in% 0]
    me_markerFlag_DF <- me_markerFlag_DF[! idx, ]

    # Format for ggplot
    me_markerFlag_DF$Cell_Subset <- factor(me_markerFlag_DF$Cell_Subset
      , levels = c("RG", "RG IP", "IP", "RG Neuron", "IP Neuron", "Neuron"))
    me_markerFlag_DF$variable <- factor(me_markerFlag_DF$variable
      , levels = c("RG eigengene", "IP eigengene", "Neuron eigengene"))
    # Determine ggplot limits from whiskers
    limits_L <- Calculate_Quantile_Limits(
      me_markerFlag_DF, value_col = "value", group_col = "Cell_Subset"
    )
    # outliers_removed_DF <- GGplot_Remove_Outliers(
    #   me_markerFlag_DF, value_col = "value", group_col = "Cell_Subset"
    # )

    print("Plotting ME boxplot...")
    gg <- ggplot(me_markerFlag_DF, aes(x = Cell_Subset, y = value)) +
      geom_boxplot(aes(fill = variable), outlier.shape = NA) +
      scale_fill_manual(values = fill_color) +
      coord_cartesian(
        ylim = c(limits_L$Minimum * 1.05, limits_L$Maximum * 1.8)) +
      stat_compare_means(value ~ Cell_Subset, data = me_markerFlag_DF
        , comparisons = list(c(1,2), c(1,3), c(2,3))
        , method = "t.test", p.adjust.method = "none", label = "p.signif"
        , label.y = c(limits_L$Maximum*c(1.2, 1.4, 1.6))
      ) +
      ggplot_set_theme_publication +
      theme(legend.position = "none") +
      theme(text = element_text(size = 12, colour = "black")) +
      ylab("Eigengene value") +
      xlab("Cell subset") +
      ggtitle(paste0("DE log2 fold change cutoff: ", name
        , "\nNumber of DE genes: ", me_markerFlag_DF$Number_DE_Genes))
      return(gg)
  })
  # ggL[4:6] <- lapply(ggL[4:6], function(gg){
  #   gg + scale_fill_manual(values = color_2)
  # })
  return(ggL)
}

# Plotting
Plot_ME_CellType_Genes_Run <- function(){

  # RG IP cluster 8
  gg1L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["RG_to_IP_8"]]
  )
  Plot_Grid(gg1L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE RG vs IP genes in"
      , "\ncluster 8 RG+, IP+, RG+IP+, Neuron-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPpNn_Cluster8.pdf")
    , height = 8, width = 12)

  # RG IP cluster 10
  gg2L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["RG_to_IP_10"]]
  )
  Plot_Grid(gg2L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE RG vs IP genes in"
      , "\ncluster 10 RG+, IP+, RG+IP+, Neuron-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPpNn_Cluster10.pdf")
    , height = 8, width = 12)

  # RG neuron cluster 8
  gg3L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["RG_to_Neuron_08"]]
  )
  Plot_Grid(gg3L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE RG vs Neuron genes in"
      , "\ncluster 0,8 RG+, Neuron+, RG+Neuron+, IP-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPnNp_Cluster08.pdf")
    , height = 8, width = 12)

  # RG neuron cluster 10
  gg4L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["RG_to_Neuron_010"]]
  )
  Plot_Grid(gg4L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE RG vs Neuron genes in"
      , "\ncluster 0,10 RG+, Neuron+, RG+Neuron+, IP-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGpIPnNp_Cluster010.pdf")
    , height = 8, width = 12)

  # IP neuron cluster 8
  gg5L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["IP_to_Neuron_08"]]
  )
  Plot_Grid(gg5L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE IP vs Neuron genes in"
      , "\ncluster 0,8 IP+, Neuron+, IP+Neuron+, RG-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGnIPpNp_Cluster08.pdf")
    , height = 8, width = 12)

  # IP neuron cluster 10
  gg6L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellType_L[["IP_to_Neuron_010"]]
  )
  Plot_Grid(gg6L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of DE IP vs Neuron genes in"
      , "\ncluster 0,10 IP+, Neuron+, IP+Neuron+, RG-")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_RGnIPpNp_Cluster010.pdf")
    , height = 8, width = 12)

  # Paper
  Append_List <- function(list_of_list_objects){
    ggL <- append(list_of_list_objects[1], list_of_list_objects[2])
    for(i in 3:length(list_of_list_objects)){
      gg <- list_of_list_objects[i]
      ggL <- append(ggL, gg)
    }
    return(ggL)
  }
  ggL <- Append_List(list_of_list_objects = c(
    # S phase
    gg1L[2]
    , gg3L[2]
    , gg5L[2]
    , gg1L[5]
    , gg3L[5]
    , gg5L[5]
    # G2/M
    , gg2L[2]
    , gg4L[2]
    , gg6L[2]
    , gg2L[5]
    , gg4L[5]
    , gg6L[5]
  ))
  # ggL <- lapply(ggL, function(gg){
  #   gg <- gg + ggtitle("") + xlab("")
  #   return(gg)
  # })
  Plot_Grid(ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nME of DE genes")
  )
  ggsave(paste0(outGraph, "DE_ME_boxplot_paper.pdf")
    , height = 12, width = 10)
}
Plot_ME_CellType_Genes_Run()

# Violin plots

Plot_ME_CellType_Genes_Violin <- function(
  me_markerFlag_DF, color_1 = NA, color_2 = NA){

  print("Plot_ME_CellType_Genes_Violin")

  me_markerFlag_DF$variable <- factor(me_markerFlag_DF$variable
    , levels = c("RG eigengene", "IP eigengene", "Neuron eigengene"))
  me_markerFlag_DF$variable <- droplevels(me_markerFlag_DF$variable)

  me_markerFlag_DFL <- split(me_markerFlag_DF
    , list(me_markerFlag_DF$Fold_Change_Cutoff
    , me_markerFlag_DF$variable)
  )

  ggL <- lapply(names(me_markerFlag_DFL), function(name){

    print(name)

    me_markerFlag_DF <- me_markerFlag_DFL[[name]]

    # Fill color of boxplots
    fill_color <- ifelse(grepl("RG", name), "#8dd3c7"
      , ifelse(grepl("IP", name), "#bebada"
        , ifelse(grepl("Neuron", name), "#fb8072", NA
    )))

    # For mixed marker cells subset to cycling cells only
    idx <- me_markerFlag_DF$Cell_Subset %in% c(
        "RG IP", "RG Neuron", "IP Neuron") &
      me_markerFlag_DF$CellID %in% names(centSO@ident)[
        centSO@ident %in% 0]
    me_markerFlag_DF <- me_markerFlag_DF[! idx, ]

    # Format for ggplot
    me_markerFlag_DF$Cell_Subset <- factor(me_markerFlag_DF$Cell_Subset
      , levels = c("RG", "RG IP", "IP", "RG Neuron", "IP Neuron", "Neuron"))
    me_markerFlag_DF$variable <- factor(me_markerFlag_DF$variable
      , levels = c("RG eigengene", "IP eigengene", "Neuron eigengene"))

    # Determine ggplot limits from whiskers
    limits_L <- Calculate_Quantile_Limits(
      me_markerFlag_DF, value_col = "value", group_col = "Cell_Subset"
    )
    # outliers_removed_DF <- GGplot_Remove_Outliers(
    #   me_markerFlag_DF, value_col = "value", group_col = "Cell_Subset"
    # )

    print("Plotting ME boxplot...")
    gg <- ggplot(me_markerFlag_DF, aes(x = Cell_Subset, y = value)) +
      geom_violin(aes(fill = variable)) +
      scale_fill_manual(values = fill_color) +
      geom_quasirandom(size = 0.1) +
      coord_cartesian(
        ylim = c(limits_L$Minimum * 1.05, limits_L$Maximum * 1.8)) +
      # stat_compare_means(value ~ Cell_Subset, data = me_markerFlag_DF
      #   , comparisons = list(c(1,2), c(1,3), c(2,3))
      #   , method = "t.test", p.adjust.method = "none", label = "p.signif"
      #   , label.y = c(limits_L$Maximum*c(1.2, 1.4, 1.6))
      # ) +
      ggplot_set_theme_publication +
      theme(legend.position = "none") +
      theme(text = element_text(size = 12, colour = "black")) +
      ylab("Eigengene value") +
      xlab("Cell subset") +
      ggtitle(paste0("DE log2 fold change cutoff: ", name
        , "\nNumber of DE genes: ", me_markerFlag_DF$Number_DE_Genes))
      return(gg)
  })
  return(ggL)
}

Plot_ME_CellType_Genes_Violin_Run <- function(){
  # RG IP cluster 8
  gg1L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["RG_to_IP_8"]]
    , color_1 = "#8dd3c7", color_2 = "#bebada"
  )
  # RG IP cluster 10
  gg2L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["RG_to_IP_10"]]
    , color_1 = "#8dd3c7", color_2 = "#bebada"
  )
  # RG neuron cluster 8
  gg3L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["RG_to_Neuron_08"]]
    , color_1 = "#8dd3c7", color_2 = "#fb8072"
  )
  # RG neuron cluster 10
  gg4L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["RG_to_Neuron_010"]]
    , color_1 = "#8dd3c7", color_2 = "#fb8072"
  )
  # IP neuron cluster 8
  gg5L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["IP_to_Neuron_08"]]
    , color_1 = "#bebada", color_2 = "#fb8072"
  )
  # IP neuron cluster 10
  gg6L <- Plot_ME_CellType_Genes_Violin(
    me_markerFlag_DF = ME_CellType_L[["IP_to_Neuron_010"]]
    , color_1 = "#bebada", color_2 = "#fb8072"
  )
  ggL <- Append_List(list_of_list_objects = c(
    # S phase
    gg1L[2]
    , gg3L[2]
    , gg5L[2]
    , gg1L[5]
    , gg3L[5]
    , gg5L[5]
    # G2/M
    , gg2L[2]
    , gg4L[2]
    , gg6L[2]
    , gg2L[5]
    , gg4L[5]
    , gg6L[5]
  ))
  Plot_Grid(ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nME of DE genes")
  )
  ggsave(paste0(outGraph, "DE_ME_violin.png")
    , height = 12, width = 12)
}
Plot_ME_CellType_Genes_Violin_Run()


# Plotting ME of cell type enriched genes
Plot_ME_CellTypeEnriched_Genes_Run <- function(){
  # RG IP cluster 8
  gg1L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["RG_to_IP_8"]]
  )
  Plot_Grid(gg1L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of RG or IP cell type enriched genes in"
      , "\ncluster 8 RG+, IP+, RG+IP+, Neuron-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGpIPpNn_Cluster8.pdf")
    , height = 8, width = 12)

  # RG IP cluster 10
  gg2L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["RG_to_IP_10"]]
  )
  Plot_Grid(gg2L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of RG vs IP cell type enriched genes in"
      , "\ncluster 10 RG+, IP+, RG+IP+, Neuron-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGpIPpNn_Cluster10.pdf")
    , height = 8, width = 12)

  # RG neuron cluster 8
  gg3L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["RG_to_Neuron_08"]]
  )
  Plot_Grid(gg3L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of RG vs Neuron cell type enriched genes in"
      , "\ncluster 0,8 RG+, Neuron+, RG+Neuron+, IP-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGpIPnNp_Cluster08.pdf")
    , height = 8, width = 12)

  # RG neuron cluster 10
  gg4L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["RG_to_Neuron_010"]]
  )
  Plot_Grid(gg4L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of RG vs Neuron cell type enriched genes in"
      , "\ncluster 0,10 RG+, Neuron+, RG+Neuron+, IP-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGpIPnNp_Cluster010.pdf")
    , height = 8, width = 12)

  # IP neuron cluster 8
  gg5L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["IP_to_Neuron_08"]]
  )
  Plot_Grid(gg5L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of IP vs Neuron cell type enriched genes in"
      , "\ncluster 0,8 IP+, Neuron+, IP+Neuron+, RG-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGnIPpNp_Cluster08.pdf")
    , height = 8, width = 12)

  # IP neuron cluster 10
  gg6L <- Plot_ME_CellType_Genes(
    me_markerFlag_DF = ME_CellTypeEnriched_L[["IP_to_Neuron_010"]]
  )
  Plot_Grid(gg6L, ncol = 3, rel_height = 0.2
    , title = paste0(graphCodeTitle
      , "\n\nME of IP vs Neuron cell type enriched genes in"
      , "\ncluster 0,10 IP+, Neuron+, IP+Neuron+, RG-")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_RGnIPpNp_Cluster010.pdf")
    , height = 8, width = 12)

  # Paper
  Append_List <- function(list_of_list_objects){
    ggL <- append(list_of_list_objects[1], list_of_list_objects[2])
    for(i in 3:length(list_of_list_objects)){
      gg <- list_of_list_objects[i]
      ggL <- append(ggL, gg)
    }
    return(ggL)
  }
  ggL <- Append_List(list_of_list_objects = c(
    # S phase
    gg1L[2]
    , gg3L[2]
    , gg5L[2]
    , gg1L[5]
    , gg3L[5]
    , gg5L[5]
    # G2/M
    , gg2L[2]
    , gg4L[2]
    , gg6L[2]
    , gg2L[5]
    , gg4L[5]
    , gg6L[5]
  ))
  Plot_Grid(ggL, ncol = 3, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\nME of cell type enriched genes")
  )
  ggsave(paste0(outGraph, "CellTypeEnriched_ME_boxplot_paper.pdf")
    , height = 12, width = 12)
}
Plot_ME_CellTypeEnriched_Genes_Run()

## Average Expression
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[
#   df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
# cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
# df2 <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]
# df2 <- melt(df2)
# df2$DE <- "Neither"
# df2$DE[df2$Var1 %in% row.names(deDF)[
#   deDF$Log2_FC_Group1_vs_Group2 > 1 & deDF$FDR < 0.05]] <- "RG"
# df2$DE[df2$Var1 %in% row.names(deDF)[
#   deDF$Log2_FC_Group1_vs_Group2 < -1 & deDF$FDR < 0.05]] <- "IP"
# df2$CellID <- df2$Var2
# df2 <- Marker_Expression_Flag(df2, df1)
# aggregate(value~DE+Cell_Subset, df2, mean)
# ggplot(df2, aes(x = Cell_Subset, y = value, fill = DE)) +
#   geom_boxplot() +
#   geom_violin() +
#   ylim(-0.5, 2)
# ggsave(paste0(outGraph, "DE_boxplot.png"))
################################################################################

### Plot DE signatures in mix state cells

## Percent of DE signature genes in mix state DE

Percent_Signature_DE_Genes <- function(
  comparison_1
  , comparison_2
  , fold_change = 0.25
  , fdr = 0.05
  , transition_state
  , phase){

  # DE genes major cell classes
  genes1 <- deDF$Gene[
    deDF$Log2_FC_Group1_vs_Group2 > 0.25 &
    deDF$FDR < 0.05 &
    deDF$Comparison == comparison_1
    ]
  # DE genes mixed marker cells
  genes2_DF <- transition_state_DE_DF[
    transition_state_DE_DF$Gene %in% genes1 &
    transition_state_DE_DF$Comparison == comparison_2
    , ]
  genes2_DF$FDR <- p.adjust(genes2_DF$Pvalue, method = "BH")
  genes2 <- genes2_DF$Gene[
    genes2_DF$Log2_FC_Group1_vs_Group2 > fold_change &
    genes2_DF$FDR < fdr
    ]

  # DE genes major cell classes
  genes3 <- deDF$Gene[
    deDF$Log2_FC_Group1_vs_Group2 < -0.25 &
    deDF$FDR < 0.05 &
    deDF$Comparison == comparison_1
    ]
  # DE genes mixed marker cells
  genes4_DF <- transition_state_DE_DF[
    transition_state_DE_DF$Gene %in% genes3 &
    transition_state_DE_DF$Comparison == comparison_2
    , ]
  genes4_DF$FDR <- p.adjust(genes4_DF$Pvalue, method = "BH")
  genes4 <- genes4_DF$Gene[
    genes4_DF$Log2_FC_Group1_vs_Group2 < fold_change &
    genes4_DF$FDR < fdr
    ]

  print(c(length(genes1), length(genes2), length(genes3), length(genes4)))
  # Calculate percent
  percents <- c(
    length(intersect(genes1, genes2))/length(genes1) * 100
    , length(intersect(genes3, genes4))/length(genes3) * 100
  )
  # Calculate porportion pvalue
  pvalues <- c(
    prop.test(length(intersect(genes1, genes2)), length(genes1), p = 0.5, correct = FALSE)$p.value
    , prop.test(length(intersect(genes3, genes4)), length(genes3), p = 0.5, correct = FALSE)$p.value
  )

  # Format
  percent_signature_DF <- data.frame(
    Percent_1 = percents
    , Percent_2 = 100-percents
    , Pvalue = pvalues
    , Signature = c(
      deDF$DE_Group[
        deDF$Comparison == comparison_1 &
        deDF$Log2_FC_Group1_vs_Group2 > 0][1]
      , deDF$DE_Group[
        deDF$Comparison == comparison_1 &
        deDF$Log2_FC_Group1_vs_Group2 < 0][1]
      )
    , Transition_State = transition_state
    , Phase = phase
  )
  return(percent_signature_DF)
}

# Percent of DE signature genes in mix state DE
ggDF <- rbind(
  # RG IP
  S_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
  )
  # RG Neuron
  , S_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
  )
  # IP neuron
  , S_phase = Percent_Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , fold_change = 0
    , fdr = 1.1
    , transition_state = "IP vs IPNeuron"
    , phase = "G2M phase"
  )
)
ggDF$Signature <- factor(ggDF$Signature, levels = c("RG", "IP", "Neuron"))
ggDF$Phase <- factor(ggDF$Phase, levels = c("S phase", "G2M phase"))
ggDF$Transition_State <- factor(ggDF$Transition_State
  , levels = c("RG vs RGIP", "RG vs RGNeuron", "IP vs IPNeuron")
)

# Plot
Plot_PercentDE_Cell_Type_DE_Genes <- function(){
  print("Plot_PercentDE_Cell_Type_DE_Genes")
  # browser()
  # Format for ggplot
  idx <- seq(1, nrow(ggDF), 2)
  ggDF$Percent_1[idx] <- ggDF$Percent_1[idx] * -1
  idx <- seq(2, nrow(ggDF), 2)
  ggDF$Percent_2[idx] <- ggDF$Percent_2[idx] * -1
  ggDF$Signature <- factor(ggDF$Signature, levels = c("Neuron", "IP", "RG"))
  ggDF <- melt(ggDF, id.vars = c("Signature", "Transition_State", "Phase", "Pvalue"))
  ggDF$Transition_State_DE <-
    gsub(".* vs ", "", ggDF$Transition_State)
  idx <- c(seq(1,12,2), seq(14,24,2))
  ggDF$Transition_State_DE[idx] <-
    gsub(" vs .*", "", ggDF$Transition_State[idx])
  ggDF$Transition_State_DE <- factor(ggDF$Transition_State_DE
    , levels = c("RGIP", "RGNeuron", "IPNeuron", "RG", "IP"))
  # Plot
  ggplot(ggDF, aes(x = Signature, y = value
    , fill = Transition_State_DE)) +
    facet_wrap(~Transition_State+Phase, scales = "free", ncol = 2) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(y = 95, label = signif(Pvalue, 2)), size = 2) +
    coord_flip() +
    ylim(c(-100,100)) +
    scale_fill_manual(
      breaks = c("RG", "IP", "Neuron", "RGIP", "RGNeuron", "IPNeuron")
      , values = c(
        "RG" = "#8DD1C6"
        , "IP" = "#BFBBDA"
        # , "Neuron" = "#F47F73"
        , "RGIP" = "#A8C6D0"
        , "RGNeuron" = "#C9AC9A"
        , "IPNeuron" = "#CB9898"
        )
      ) +
    ggplot_set_theme_publication +
    ggtitle(paste0(graphCodeTitle
      , "\nPercent positive or negative fold change in mixed state cells"
      , "\nof top cell type DE genes"))
}
Plot_PercentDE_Cell_Type_DE_Genes()
ggsave(paste0(outGraph, "DE_Mixed_Percent_barplot.pdf")
  , width = 7, height = 5)


## DE of DE signature genes in mix state DE

Signature_DE_Genes <- function(
  comparison_1, comparison_2, transition_state, phase, fold_change_cutoff, DE_group_label){

  if (fold_change_cutoff > 0) {
    fold_change_cutoff_idx <- deDF$Log2_FC_Group1_vs_Group2 > fold_change_cutoff
  } else if (fold_change_cutoff < 0) {
    fold_change_cutoff_idx <- deDF$Log2_FC_Group1_vs_Group2 < fold_change_cutoff
  }
  genes1 <- deDF$Gene[
    fold_change_cutoff_idx &
    deDF$FDR < 0.05 &
    deDF$Comparison == comparison_1
    ]
  subset_transition_state_DE_DF <- transition_state_DE_DF[
    transition_state_DE_DF$Gene %in% genes1 &
    transition_state_DE_DF$Comparison == comparison_2,
    ]

  subset_transition_state_DE_DF <- subset_transition_state_DE_DF[
    order(subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2), ]
  subset_transition_state_DE_DF$Gene <- factor(
    subset_transition_state_DE_DF$Gene
    , levels = subset_transition_state_DE_DF$Gene
  )

  subset_transition_state_DE_DF$DE_Enriched <- as.factor(DE_group_label)
  subset_transition_state_DE_DF$Transition_State <- as.factor(transition_state)
  subset_transition_state_DE_DF$Phase <- as.factor(phase)

  return(subset_transition_state_DE_DF)
}

ggDFL <- list(
  Signature_DE_Genes(
      comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "RG"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "RG"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "IP"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "IP"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "RG"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "RG"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "Neuron"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "Neuron"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "IP"
  )
  , Signature_DE_Genes(
      comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "G2M phase"
    , fold_change_cutoff = 0.25
    , DE_group_label = "IP"
  )
  , Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "Neuron"
  )
  , Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "G2M phase"
    , fold_change_cutoff = -0.25
    , DE_group_label = "Neuron"
  )
)
ggL <- lapply(ggDFL, function(ggDF){
  ggplot(ggDF, aes(x = Gene, y = Log2_FC_Group1_vs_Group2)) +
    geom_bar(stat = "identity") +
    ylab("Log2 fold change") +
    xlab("Genes") +
    coord_cartesian(ylim = c(-1, 1)) +
    ggplot_set_theme_publication +
    theme(axis.text.x = element_blank()
      , axis.ticks.x = element_blank()
    ) +
    ggtitle(paste0(
      "\nTransition state: ", ggDF$Transition_State
      , "\nCell cycle phase: ", ggDF$Phase
      , "\nDE Signature: ", ggDF$DE_Enriched
    ))
})
Plot_Grid(ggL, ncol = 4, title = "DE Mixed")
ggsave(paste0(outGraph, "DE_Mixed_barplot.png"), width = 11, height = 9)


## Fold change of DE signature genes in mixed state cells barplot

Fold_Change_Signature_DE_Genes <- function(
  comparison_1, comparison_2, transition_state, phase, fold_change_cutoff, cell_types, transition_states, signature){

  # Subset to cell type DE genes and calculate mean fold change
  if (fold_change_cutoff > 0) {
    fold_change_cutoff_idx <- deDF$Log2_FC_Group1_vs_Group2 > fold_change_cutoff
  } else if (fold_change_cutoff < 0) {
    fold_change_cutoff_idx <- deDF$Log2_FC_Group1_vs_Group2 < fold_change_cutoff
  }
  subset_deDF <- deDF[
    fold_change_cutoff_idx &
    deDF$FDR < 0.05 &
    deDF$Comparison == comparison_1
    , ]
  genes1 <- subset_deDF$Gene
  # Mean fold change
  mn_cell_types <- mean(subset_deDF$Log2_FC_Group1_vs_Group2)
  # Change negative fold changes to positive
  if (mn_cell_types < 0){mn_cell_types <- mn_cell_types * -1}

  # Subset transition state cells to cell type DE genes and calculate mean fold change
  subset_transition_state_DE_DF <- transition_state_DE_DF[
    transition_state_DE_DF$Gene %in% genes1 &
    transition_state_DE_DF$Comparison == comparison_2
    , ]
  # Mean fold change
  mn_transition_state <-
    mean(subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2)
  mn_transition_state <- mn_transition_state * -1
  # # Change negative fold changes to positive
  # if (mn_transition_state < 0){
  #   mn_transition_state <- mn_transition_state * -1}

  # Format
  percent_DF <- data.frame(
    Mean_Transition_State = mn_transition_state
    , Mean_Cell_Types = mn_cell_types
    , Percent_Fold_Change = mn_transition_state / mn_cell_types * 100
    , Signature = signature
    , Cell_Types = comparison_1
    , Transition_States = transition_states
    , Phase = phase
  )
  percent_DF$Signature <- factor(percent_DF$Signature
    , levels = rev(c("RG", "IP", "Neuron")))
  percent_DF$Phase <- factor(percent_DF$Phase
    , levels = c("S phase", "G2M phase"))
  percent_DF$Transition_States <- factor(percent_DF$Transition_States
    , levels = c("RG vs RGIP", "RG vs RGNeuron", "IP vs IPNeuron"))
  percent_DF$Cell_Types <- factor(percent_DF$Cell_Types
    , levels = c("RG_vs_IP", "RG_vs_Neuron", "IP_vs_Neuron"))
  return(percent_DF)
}

# Fold change of cell type enriched genes in mixed state cells versus cell types
percent_fold_change_DFL <- list(
  # RG IP
  Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , fold_change_cutoff = 0.25
    , transition_states = "RG vs RGIP"
    , phase = "S phase"
    , signature = "RG"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , fold_change_cutoff = -0.25
    , transition_states = "RG vs RGIP"
    , phase = "S phase"
    , signature = "IP"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , fold_change_cutoff = 0.25
    , transition_states = "RG vs RGIP"
    , phase = "G2M phase"
    , signature = "RG"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , fold_change_cutoff = -0.25
    , transition_states = "RG vs RGIP"
    , phase = "G2M phase"
    , signature = "IP"
  )
  # RG neuron
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , fold_change_cutoff = 0.25
    , transition_states = "RG vs RGNeuron"
    , phase = "S phase"
    , signature = "RG"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , fold_change_cutoff = -0.25
    , transition_states = "RG vs RGNeuron"
    , phase = "S phase"
    , signature = "Neuron"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , fold_change_cutoff = 0.25
    , transition_states = "RG vs RGNeuron"
    , phase = "G2M phase"
    , signature = "RG"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , fold_change_cutoff = -0.25
    , transition_states = "RG vs RGNeuron"
    , phase = "G2M phase"
    , signature = "Neuron"
  )
  # IP neuron
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , fold_change_cutoff = 0.25
    , transition_states = "IP vs IPNeuron"
    , phase = "S phase"
    , signature = "IP"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , fold_change_cutoff = -0.25
    , transition_states = "IP vs IPNeuron"
    , phase = "S phase"
    , signature = "Neuron"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , fold_change_cutoff = 0.25
    , transition_states = "IP vs IPNeuron"
    , phase = "G2M phase"
    , signature = "IP"
  )
  , Fold_Change_Signature_DE_Genes(
    comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , fold_change_cutoff = -0.25
    , transition_states = "IP vs IPNeuron"
    , phase = "G2M phase"
    , signature = "Neuron"
  )
)
percent_fold_change_DF <- do.call("rbind", percent_fold_change_DFL)
# Plot
ggplot(percent_fold_change_DF, aes(x = Signature, y = Percent_Fold_Change
  , fill = Signature)) +
  facet_wrap(~Cell_Types+Transition_States+Phase, scales = "free", ncol = 2) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  ylim(c(-100,100)) +
  scale_fill_manual(
    breaks = c("RG", "IP", "Neuron")
    , values = c(
      "RG" = "#8DD1C6"
      , "IP" = "#BFBBDA"
      , "Neuron" = "#F47F73"
      )
    ) +
  ggplot_set_theme_publication +
  # scale_fill_manual(values = c(Neuron = "#fb8072", IP = "#bebada", RG = "#8dd3c7")) +
  ggtitle(paste0(graphCodeTitle
    , "\nDE fold change in mixed state cells versus cell types"
    , "\nof top cell type enriched genes"))
ggsave(paste0(outGraph
  , "DE_Mixed_PercentFoldChange_barplot.pdf")
  , width = 7, height = 5)
################################################################################

### DE of top eniched genes in cell type

Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes <- function(
  seurat_clusters, DE_comparison_2, transition_state, phase, fold_change_cutoff, cell_type_label){

  # Subset genes to top enriched for Seurat clusters of interest
  genes1 <- cluster_DE_DF$Gene[cluster_DE_DF$Cluster == seurat_clusters]
  # Subset by genes and mixed state comparison
  subset_transition_state_DE_DF <- transition_state_DE_DF[
    transition_state_DE_DF$Gene %in% genes1 &
    transition_state_DE_DF$Comparison == DE_comparison_2,
    ]
  # Format
  subset_transition_state_DE_DF <- subset_transition_state_DE_DF[
    order(subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2), ]
  subset_transition_state_DE_DF$Gene <- factor(
    subset_transition_state_DE_DF$Gene
    , levels = subset_transition_state_DE_DF$Gene
  )
  subset_transition_state_DE_DF$DE_Enriched <- as.factor(cell_type_label)
  subset_transition_state_DE_DF$Transition_State <- as.factor(transition_state)
  subset_transition_state_DE_DF$Phase <- as.factor(phase)

  return(subset_transition_state_DE_DF)
}

subset_transition_state_DE_DFL <- list(
  Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(7,9)
    , DE_comparison_2 = "RG_vs_RGIP_Sphase"
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
    , cell_type_label = "RG"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(7,9)
    , DE_comparison_2 = "RG_vs_RGIP_G2Mphase"
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
    , cell_type_label = "RG"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(2)
    , DE_comparison_2 = "RG_vs_RGIP_Sphase"
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
    , cell_type_label = "IP"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(2)
    , DE_comparison_2 = "RG_vs_RGIP_G2Mphase"
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
    , cell_type_label = "IP"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(7,9)
    , DE_comparison_2 = "RG_vs_RGNeuron_Sphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
    , cell_type_label = "RG"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(7,9)
    , DE_comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
    , cell_type_label = "RG"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(0)
    , DE_comparison_2 = "RG_vs_RGNeuron_Sphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
    , cell_type_label = "Neuron"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(0)
    , DE_comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
    , cell_type_label = "Neuron"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(2)
    , DE_comparison_2 = "IP_vs_IPNeuron_Sphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
    , cell_type_label = "IP"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(2)
    , DE_comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "G2M phase"
    , cell_type_label = "IP"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(0)
    , DE_comparison_2 = "IP_vs_IPNeuron_Sphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
    , cell_type_label = "Neuron"
  )
  , Subset_Transition_State_DE_By_Cell_Type_Enriched_Genes(
      seurat_clusters = c(0)
    , DE_comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , transition_state = "IP vs IPNeuron"
    , phase = "G2M phase"
    , cell_type_label = "Neuron"
  )
)
subset_transition_state_DE_DF <-
  do.call("rbind", subset_transition_state_DE_DFL)

## Barplot
ggL <- lapply(subset_transition_state_DE_DFL, function(ggDF){
  ggplot(ggDF, aes(x = Gene, y = Log2_FC_Group1_vs_Group2)) +
    geom_bar(stat = "identity") +
    ylab("Log2 fold change") +
    xlab("Genes") +
    coord_cartesian(ylim = c(-1, 1)) +
    ggplot_set_theme_publication +
    theme(axis.text.x = element_blank()
      , axis.ticks.x = element_blank()
    ) +
    ggtitle(paste0(
      "\nTransition state: ", ggDF$Transition_State
      , "\nCell cycle phase: ", ggDF$Phase
      , "\nDE Signature: ", ggDF$DE_Enriched
    ))
})
Plot_Grid(ggL, ncol = 4, title = "DE in mixed state cells of top cell type enriched genes")
ggsave(paste0(outGraph, "TopCellTypeEnriched_Mixed_barplot.png")
  , width = 11, height = 18)


## Percentage of enriched genes DE in mixed state cells barplot

PercentDE_Cell_Type_Enriched_Genes <- function(){
  percent_DE_genes_LDF <- lapply(subset_transition_state_DE_DFL
    , function(subset_transition_state_DE_DF){
    print("PercentDE_Cell_Type_Enriched_Genes")
    v1 <- subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2 > 0
    pval1 <- prop.test(sum(v1), length(v1), p = 0.5, correct = FALSE)$p.value
    v2 <- subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2 < 0
    pval2 <- prop.test(sum(v2), length(v2), p = 0.5, correct = FALSE)$p.value
    percent_DF <- as.data.frame(rbind(Percent_Of_Table(v1), Percent_Of_Table(v2)))
    percent_DF$Pvalue <- c(pval1, pval2)
    percent_DF$Signature <- subset_transition_state_DE_DF$DE_Enriched[1]
    percent_DF$Transition_State <-
      subset_transition_state_DE_DF$Transition_State[1]
    percent_DF$Phase <- subset_transition_state_DE_DF$Phase[1]
    return(percent_DF)
  })
  percent_DE_genes_DF <- do.call("rbind", percent_DE_genes_LDF)
  return(percent_DE_genes_DF)
}

Plot_PercentDE_Cell_Type_Enriched_Genes <- function(){
  print("Plot_PercentDE_Cell_Type_Enriched_Genes")
  # Calculate percent
  percent_DE_genes_DF <- PercentDE_Cell_Type_Enriched_Genes()
  # Format for ggplot
  percent_DE_genes_DF$Transition_State_DE <-
    gsub(".* vs ", "", percent_DE_genes_DF$Transition_State)
  idx <- seq(1, nrow(percent_DE_genes_DF)-1, 2)
  percent_DE_genes_DF$Transition_State_DE <-
    gsub(" vs .*", "", percent_DE_genes_DF$Transition_State)
  percent_DE_genes_DF$Transition_State_DE[idx] <-
    gsub(".* vs ", "", percent_DE_genes_DF$Transition_State[idx])
  percent_DE_genes_DF$Transition_State_DE <- factor(
    percent_DE_genes_DF$Transition_State_DE
    , levels = rev(c("RG", "IP", "Neuron", "RGIP", "RGNeuron", "IPNeuron"))
  )
  colnames(percent_DE_genes_DF)[c(1,2)] <- c("Percent_1", "Percent_2")
  idx <- seq(2, nrow(percent_DE_genes_DF), 2)
  percent_DE_genes_DF$Percent_1[idx] <- percent_DE_genes_DF$Percent_1[idx] * -1
  percent_DE_genes_DF$Signature <- factor(
    percent_DE_genes_DF$Signature, levels = c("Neuron", "IP", "RG"))
  # Plot
  ggplot(percent_DE_genes_DF, aes(x = Signature, y = Percent_1
    , fill = Transition_State_DE)) +
    facet_wrap(~Transition_State+Phase, scales = "free", ncol = 2) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(y = 95, label = signif(Pvalue, 2)), size = 2) +
    coord_flip() +
    scale_fill_manual(
      breaks = c("RG", "IP", "Neuron", "RGIP", "RGNeuron", "IPNeuron")
      , values = c(
        "RG" = "#8DD1C6"
        , "IP" = "#BFBBDA"
        # , "Neuron" = "#F47F73"
        , "RGIP" = "#A8C6D0"
        , "RGNeuron" = "#C9AC9A"
        , "IPNeuron" = "#CB9898"
        )
      ) +
    ylim(c(-100,100)) +
    ggplot_set_theme_publication +
    ggtitle(paste0(graphCodeTitle
      , "\nPercent positive or negative fold change in mixed state cells"
      , "\nof top cell type enriched genes"))
}
Plot_PercentDE_Cell_Type_Enriched_Genes()
ggsave(paste0(outGraph, "TopCellTypeEnriched_Mixed_Percent_barplot.pdf")
  , width = 7, height = 5)

## Fold change of enriched genes DE in mixed state cells barplot

Fold_Change_Cell_Type_Enriched_Genes <- function(
  signature, seurat_clusters, cell_types, transition_states, phase){

    print("Fold_Change_Cell_Type_Enriched_Genes")

    # Subset mixed state DE
    subset_transition_state_DE_DF <- subset_transition_state_DE_DF[
      with(subset_transition_state_DE_DF
        , Transition_State == transition_states
        & DE_Enriched == signature
        & Phase == phase
      ), ]
    # Mean fold change
    mn_transition_state <-
      mean(subset_transition_state_DE_DF$Log2_FC_Group1_vs_Group2)
    # Change negative fold changes to positive
    if (mn_transition_state < 0){
      mn_transition_state <- mn_transition_state * -1}

    # Subset genes to top enriched for Seurat clusters of interest
    genes1 <- cluster_DE_DF$Gene[cluster_DE_DF$Cluster == seurat_clusters]
    subset_DE_DF <- deDF[! is.na(deDF$DE_Group), ]
    subset_DE_DF <- subset_DE_DF[
      with(subset_DE_DF
        , Comparison == cell_types &
        DE_Group == signature &
        Gene %in% genes1
      ), ]
    # Mean fold change
    mn_cell_types <- mean(subset_DE_DF$Log2_FC_Group1_vs_Group2)
    mn_cell_types <- mn_cell_types * -1
    # # Change negative fold changes to positive
    # if (mn_cell_types < 0){mn_cell_types <- mn_cell_types * -1}

    # Format
    percent_DF <- data.frame(
      Mean_Transition_State = mn_transition_state
      , Mean_Cell_Types = mn_cell_types
      , Percent_Fold_Change = mn_transition_state / mn_cell_types * 100
      , Signature = signature
      , Cell_Types = cell_types
      , Transition_States = transition_states
      , Phase = phase
    )
    return(percent_DF)
}

# Fold change of cell type enriched genes in mixed state cells versus cell types
percent_fold_change_DFL <- list(
  # RG IP
  Fold_Change_Cell_Type_Enriched_Genes(
    signature = "RG"
    , seurat_clusters = c(7,9)
    , cell_types = "RG_vs_IP"
    , transition_states = "RG vs RGIP"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "RG"
    , seurat_clusters = c(7,9)
    , cell_types = "RG_vs_IP"
    , transition_states = "RG vs RGIP"
    , phase = "G2M phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "IP"
    , seurat_clusters = c(2)
    , cell_types = "RG_vs_IP"
    , transition_states = "RG vs RGIP"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "IP"
    , seurat_clusters = c(2)
    , cell_types = "RG_vs_IP"
    , transition_states = "RG vs RGIP"
    , phase = "G2M phase"
  )
  # RG Neuron
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "RG"
    , seurat_clusters = c(7,9)
    , cell_types = "RG_vs_Neuron"
    , transition_states = "RG vs RGNeuron"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "RG"
    , seurat_clusters = c(7,9)
    , cell_types = "RG_vs_Neuron"
    , transition_states = "RG vs RGNeuron"
    , phase = "G2M phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "Neuron"
    , seurat_clusters = c(0)
    , cell_types = "RG_vs_Neuron"
    , transition_states = "RG vs RGNeuron"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "Neuron"
    , seurat_clusters = c(0)
    , cell_types = "RG_vs_Neuron"
    , transition_states = "RG vs RGNeuron"
    , phase = "G2M phase"
  )
  # IP Neuron
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "IP"
    , seurat_clusters = c(2)
    , cell_types = "IP_vs_Neuron"
    , transition_states = "IP vs IPNeuron"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "IP"
    , seurat_clusters = c(2)
    , cell_types = "IP_vs_Neuron"
    , transition_states = "IP vs IPNeuron"
    , phase = "G2M phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "Neuron"
    , seurat_clusters = c(0)
    , cell_types = "IP_vs_Neuron"
    , transition_states = "IP vs IPNeuron"
    , phase = "S phase"
  )
  , Fold_Change_Cell_Type_Enriched_Genes(
    signature = "Neuron"
    , seurat_clusters = c(0)
    , cell_types = "IP_vs_Neuron"
    , transition_states = "IP vs IPNeuron"
    , phase = "G2M phase"
  )
)
percent_fold_change_DF <- do.call("rbind", percent_fold_change_DFL)
# Format
percent_fold_change_DF$Signature <- factor(percent_fold_change_DF$Signature
  , levels = rev(c("RG", "IP", "Neuron")))
percent_fold_change_DF$Phase <- factor(percent_fold_change_DF$Phase
  , levels = c("S phase", "G2M phase"))
percent_fold_change_DF$Transition_States <- factor(percent_fold_change_DF$Transition_States
  , levels = c("RG vs RGIP", "RG vs RGNeuron", "IP vs IPNeuron"))
percent_fold_change_DF$Cell_Types <- factor(percent_fold_change_DF$Cell_Types
  , levels = c("RG_vs_IP", "RG_vs_Neuron", "IP_vs_Neuron"))
# Plot
ggplot(percent_fold_change_DF, aes(x = Signature, y = Percent_Fold_Change
  , fill = Signature)) +
  facet_wrap(~Cell_Types+Transition_States+Phase, scales = "free", ncol = 2) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  scale_fill_manual(
    breaks = c("RG", "IP", "Neuron")
    , values = c(
      "RG" = "#8DD1C6"
      , "IP" = "#BFBBDA"
      , "Neuron" = "#F47F73"
      )
    ) +
  ylim(c(-100,100)) +
  ggplot_set_theme_publication +
  ggtitle(paste0(graphCodeTitle
    , "\nDE fold change in mixed state cells versus cell types"
    , "\nof top cell type enriched genes"))
ggsave(paste0(outGraph
  , "TopCellTypeEnriched_Mixed_PercentFoldChange_barplot.pdf")
  , width = 7, height = 5)
################################################################################

### Percent of cells passing expression filters per cluster

Mixed_Marker_By_Cluster_Percent <- function(
  exM
  , seuratO
  , highThreshold
  , lowThreshold
  , type_keep
  , cluster_order = NULL){
    # browser()
    print("Mixed_Marker_By_Cluster_Percent")
    expr_flag_DF <- Average_MarkersExp_Per_Cell(exM = exM, seuratO = seuratO)
    expr_flag_DF <- Positive_Negative_Expression_Flag(
      exDF = expr_flag_DF
      , highThreshold = highThreshold
      , lowThreshold = lowThreshold
    )
    expr_flag_DF$Type_Plot <- expr_flag_DF$TYPE
    expr_flag_DF$Type_Plot[! expr_flag_DF$Type_Plot %in% type_keep] <- "NA"
    expr_flag_DF$Type_Plot <- droplevels(expr_flag_DF$Type_Plot)
    expr_flag_DF <- aggregate(expr_flag_DF$Type_Plot
      , list(expr_flag_DF$CLUSTER), Percent_Of_Table)
    tmp_df <- as.data.frame(expr_flag_DF[ ,2])
    tmp_df$CLUSTER <- expr_flag_DF[ ,1]
    expr_flag_DF <- melt(tmp_df)
    return(expr_flag_DF)
}

Mixed_Marker_By_Cluster_Percent_Barplot <- function(
  exM
  , seuratO
  , highThreshold
  , lowThreshold
  , type_keep
  , title
  , cluster_order = NULL) {
  # browser()
  print("Mixed_Marker_By_Cluster_Percent_Barplot")
  expr_flag_DF <- Mixed_Marker_By_Cluster_Percent(
    exM = exM
    , seuratO = seuratO
    , highThreshold = highThreshold
    , lowThreshold = lowThreshold
    , type_keep = type_keep
    , cluster_order = cluster_order
  )
  # Format
  # Set cluster order
  if (! is.null(cluster_order)){
    expr_flag_DF$CLUSTER <- factor(as.character(expr_flag_DF$CLUSTER)
      , levels = as.character(cluster_order))
  }
  expr_flag_DF <- expr_flag_DF[! is.na(expr_flag_DF$CLUSTER), ]
  # Set marker plotting order
  expr_flag_DF$variable <- as.character(expr_flag_DF$variable)
  expr_flag_DF$variable[is.na(expr_flag_DF$variable)] <- "Unknown"
  expr_flag_DF$variable <- factor(expr_flag_DF$variable
    , levels = rev(
      c("RG+", "IP+", "Neuron+", "IP+ RG+", "Neuron+ RG+"
      , "Neuron+ IP+", "Interneuron+ IP+", "Endothelial+ IP+", "Unknown")
    )
  )
  # Plot
  gg <- ggplot(expr_flag_DF, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    scale_fill_manual(
      name = "Cell type markers"
      # Set legend order and colors
      , breaks = c(
        "RG+"
        , "IP+"
        , "Neuron+"
        , "IP+ RG+"
        , "Neuron+ RG+"
        , "Neuron+ IP+"
        , "Interneuron+ IP+"
        , "Endothelial+ IP+"
        , "Unknown"
      )
      , values = c(
        "RG+" = "#8DD1C6"
        , "IP+" = "#BFBBDA"
        , "Neuron+" = "#F47F73"
        , "IP+ RG+" = "#A8C6D0"
        , "Neuron+ RG+" = "#C9AC9A"
        , "Neuron+ IP+" = "#CB9898"
        , "Interneuron+ IP+" = "#A1B5D6"
        , "Endothelial+ IP+" = "#BAC8A9"
        , "Unknown" = "#E7E6E5")
      ) +
    xlab("Cluster") +
    ylab("Percent of cells") +
    ggtitle(title)
  return(gg)
}

# Plot
gg1 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+")
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+")
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.75, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+")
  , title = "Keep CC\n+ = > 0.75 normalized expression\n- = < 0.25 normalized expression"
)
gg4 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.4, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+")
  , title = "Keep CC\n+ = > 0.4 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, gg4, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of marker expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot.pdf")
  , width = 12, height = 14)

# Plot
gg1 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+", "RG+", "IP+", "Neuron+", "Neuron+ RG+ IP+")
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+", "RG+", "IP+", "Neuron+", "Neuron+ RG+ IP+")
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.75, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+", "RG+", "IP+", "Neuron+", "Neuron+ RG+ IP+")
  , title = "Keep CC\n+ = > 0.75 normalized expression\n- = < 0.25 normalized expression"
)
gg4 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.4, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+", "RG+", "IP+", "Neuron+")
  , title = "Keep CC\n+ = > 0.4 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, gg4, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of marker expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentMixedCluster2_Barplot.pdf")
  , width = 12, height = 14)

# Paper - plot with reordered clusters
gg <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "Endothelial+ IP+"
    , "Interneuron+ IP+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  , title = paste0(graphCodeTitle
    , "\n\nPercent of cells passing combinations of marker expression filters"
    , "\nKeep CC"
    , "\n+ = > 0.5 normalized expression"
    , "\n- = < 0.25 normalized expression"
  )
)
gg + ggplot_set_theme_publication
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot_paper.pdf")
  , width = 5, height = 3)

# Paper - plot with reordered clusters
gg <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  , title = paste0(graphCodeTitle
    , "\n\nPercent of cells passing combinations of marker expression filters"
    , "\nKeep CC"
    , "\n+ = > 0.5 normalized expression"
    , "\n- = < 0.25 normalized expression"
  )
)
gg + ggplot_set_theme_publication
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot_paper_supp.pdf")
  , width = 5, height = 3)

# Paper - table of percents for making differentiation diagram
percentDF <- Mixed_Marker_By_Cluster_Percent(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
percentDF <- percentDF[with(percentDF, CLUSTER == 10), ]
percentDF$Percent_NA_Removed <-
  c(NA, (percentDF$value[2:6] / sum(percentDF$value[2:6])) * 100, NA)
percentDF[ ,c(3:4)] <- round(percentDF[ ,c(3:4)], 1)
write.csv(percentDF, file = paste0(
    outTable, "PercentMixedCluster_paper_supp.csv")
  , quote = FALSE)

Mixed_Marker_By_Cluster_Numbers <- function(
  exM
  , seuratO
  , highThreshold
  , lowThreshold
  , type_keep
  , cluster_order = NULL){
    # browser()
    print("Mixed_Marker_By_Cluster_Numbers")
    expr_flag_DF <- Average_MarkersExp_Per_Cell(exM = exM, seuratO = seuratO)
    expr_flag_DF <- Positive_Negative_Expression_Flag(
      exDF = expr_flag_DF
      , highThreshold = highThreshold
      , lowThreshold = lowThreshold
    )
    expr_flag_DF$Type_Plot <- expr_flag_DF$TYPE
    expr_flag_DF$Type_Plot[! expr_flag_DF$Type_Plot %in% type_keep] <- "NA"
    expr_flag_DF$Type_Plot <- droplevels(expr_flag_DF$Type_Plot)
    numbers_DF <- with(expr_flag_DF, table(Type_Plot, CLUSTER))
    numbers_DF <- as.data.frame(numbers_DF)
    return(numbers_DF)
}
number_DF <- Mixed_Marker_By_Cluster_Numbers(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
write.csv(number_DF, file = paste0(
    outTable, "Number_Mixed_Cluster.csv")
  , quote = FALSE)
################################################################################

### Percent of cells passing expression filters of cell type enriched genes

Mixed_CellTypeEnriched_By_Cluster_Percent <- function(
  exM
  , seuratO
  , highThreshold
  , lowThreshold
  , type_keep
  , fold_change
  , cluster_order = NULL){
    # browser()
    print("Mixed_CellTypeEnriched_By_Cluster_Percent")
    expr_flag_DF <- Format_Mean_CellTypeEnriched_Expression(
      exM = exM
      , seuratO = seuratO
      , fold_change = fold_change
    )
    expr_flag_DF <- Positive_Negative_Expression_Flag(
      exDF = expr_flag_DF
      , highThreshold = highThreshold
      , lowThreshold = lowThreshold
    )
    expr_flag_DF$Type_Plot <- expr_flag_DF$TYPE
    expr_flag_DF$Type_Plot[! expr_flag_DF$Type_Plot %in% type_keep] <- "NA"
    expr_flag_DF$Type_Plot <- droplevels(expr_flag_DF$Type_Plot)
    expr_flag_DFL <- tapply(
      expr_flag_DF$Type_Plot, expr_flag_DF$CLUSTER, Percent_Of_Table
    )
    expr_flag_DF <- do.call("rbind", expr_flag_DFL)
    # expr_flag_DF <- aggregate(expr_flag_DF$Type_Plot
    #   , list(expr_flag_DF$CLUSTER), Percent_Of_Table)
    # tmp_df <- as.data.frame(expr_flag_DF[ ,2])
    # tmp_df$CLUSTER <- expr_flag_DF[ ,1]
    expr_flag_DF <- as.data.frame(expr_flag_DF)
    expr_flag_DF$CLUSTER <- rownames(expr_flag_DF)
    # Set cluster order
    if (! is.null(cluster_order)){
      expr_flag_DF$CLUSTER <- factor(as.character(expr_flag_DF$CLUSTER)
        , levels = as.character(cluster_order))
    }
    expr_flag_DF <- expr_flag_DF[! is.na(expr_flag_DF$CLUSTER), ]
    expr_flag_DF <- melt(expr_flag_DF)
    return(expr_flag_DF)
}

Mixed_CellTypeEnriched_By_Cluster_Percent_Barplot <- function(
  exM
  , seuratO
  , highThreshold
  , lowThreshold
  , fold_change
  , type_keep
  , title
  , cluster_order = NULL) {
  # browser()
  print("Mixed_CellTypeEnriched_By_Cluster_Percent_Barplot")
  expr_flag_DF <- Mixed_CellTypeEnriched_By_Cluster_Percent(
    exM = exM
    , seuratO = seuratO
    , highThreshold = highThreshold
    , lowThreshold = lowThreshold
    , type_keep = type_keep
    , cluster_order = cluster_order
    , fold_change = fold_change
  )
  # Format
  # Set marker plotting order
  expr_flag_DF$variable <- as.character(expr_flag_DF$variable)
  expr_flag_DF$variable[is.na(expr_flag_DF$variable)] <- "Unknown"
  expr_flag_DF$variable <- factor(expr_flag_DF$variable
    , levels = rev(
      c("RG+", "IP+", "Neuron+", "IP+ RG+", "Neuron+ RG+"
      , "Neuron+ IP+", "Interneuron+ IP+", "Endothelial+ IP+", "Unknown")
    )
  )
  # Plot
  gg <- ggplot(expr_flag_DF, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    scale_fill_manual(
      name = "Cell type markers"
      # Set legend order and colors
      , breaks = c(
        "RG+"
        , "IP+"
        , "Neuron+"
        , "IP+ RG+"
        , "Neuron+ RG+"
        , "Neuron+ IP+"
        , "Interneuron+ IP+"
        , "Endothelial+ IP+"
        , "Unknown"
      )
      , values = c(
        "RG+" = "#8DD1C6"
        , "IP+" = "#BFBBDA"
        , "Neuron+" = "#F47F73"
        , "IP+ RG+" = "#A8C6D0"
        , "Neuron+ RG+" = "#C9AC9A"
        , "Neuron+ IP+" = "#CB9898"
        , "Interneuron+ IP+" = "#A1B5D6"
        , "Endothelial+ IP+" = "#BAC8A9"
        , "Unknown" = "#E7E6E5")
      ) +
    xlab("Cluster") +
    ylab("Percent of cells") +
    ggtitle(title)
  return(gg)
}


Mean_CellTypeEnriched_Expression <- function(exM, clusters, fold_change ){
  print("Mean_CellTypeEnriched_Expression")
  genes <- cluster_DE_DF$Gene[
    cluster_DE_DF$Cluster %in% clusters
    & cluster_DE_DF$Log2_Fold_Change > fold_change
    ]
  mean_expr <- colMeans(
    exM[row.names(exM) %in% genes, ]
  )
  return(mean_expr)
}

Format_Mean_CellTypeEnriched_Expression <- function(exM, seuratO, fold_change){
  print("Format_Mean_CellTypeEnriched_Expression")
  # Mean expression of cell type enriched genes
  mnExDF <- data.frame(
    RG = Mean_CellTypeEnriched_Expression(
      exM = exM, clusters = c(7,9), fold_change = fold_change
    )
    , IP = Mean_CellTypeEnriched_Expression(
      exM = exM, clusters = 2, fold_change = fold_change
    )
    , Neuron = Mean_CellTypeEnriched_Expression(
      exM = exM, clusters = 0, fold_change = fold_change
    )
    , vRG = 0
    , oRG = 0
    , Endothelial = 0
    , Interneuron = 0
  )
  # Add metadata
  idx <- match(row.names(mnExDF), row.names(seuratO@meta.data))
  mnExDF$PHASE <- seuratO@meta.data$Phase[idx]
  mnExDF$CLUSTER <- seuratO@ident[idx]
  mnExDF$nUMI <- seuratO@meta.data$nUMI[idx]
  mnExDF$G2Mscore <- seuratO@meta.data$G2M.Score[idx]
  mnExDF$Sscore <- seuratO@meta.data$S.Score[idx]
  return(mnExDF)
}


expr_flag_DF <- Format_Mean_CellTypeEnriched_Expression(
  exM = noCentExM
  , seuratO = centSO
  , fold_change = 0.4
)
ggDF <- expr_flag_DF[ ,c(1:3, 9)]
ggDF <- melt(ggDF)
stdev_DF <- data.frame(
  Stdev = c(
    sd(ggDF$value[ggDF$CLUSTER %in% c(7,9)])
    , sd(ggDF$value[ggDF$CLUSTER %in% c(2)])
    , sd(ggDF$value[ggDF$CLUSTER %in% c(0)])
  )
  , Median = c(
    median(ggDF$value[ggDF$CLUSTER %in% c(7,9)])
    , median(ggDF$value[ggDF$CLUSTER %in% c(2)])
    , median(ggDF$value[ggDF$CLUSTER %in% c(0)])
  )
  , variable = c("RG", "IP", "Neuron")
)
stdev_DF$Median - stdev_DF$Stdev
ggplot(ggDF, aes(x = CLUSTER, y = value)) +
  facet_wrap(~variable) +
  geom_violin() +
  geom_jitter(size = 0.01, alpha = 0.1)
ggsave(paste0(outGraph, "CellTypeEnriched_violin.png"), width = 9, height = 6
 , dpi = 150)



 Positive_Negative_Expression_Flag <- function(
   exDF, highThreshold, lowThreshold) {

   # Expected input data frame format
   # (from Average_MarkersExp_Per_Cell)
   # vRG         oRG          RG          IP Endothelial
   # ACCAGCTAGCCT -0.12961098 -0.04274303 -0.16390677 -0.05747777  0.08723860
   # GGAAGGACTGCA  0.17994922 -0.08862593  0.47211654 -0.02872334  1.65468242
   # Neuron PHASE CLUSTER nUMI
   # ACCAGCTAGCCT  0.8625643    G1       3 6165
   # GGAAGGACTGCA -0.9443603   G2M      11 4166

   df <- exDF
   df$TYPE <- NA
   df$TYPE[df[ ,c("Neuron")] > 1] <- "Neuron+"
   df$TYPE[df[ ,c("IP")] > 1.3] <- "IP+"
   df$TYPE[df[ ,c("RG")] > 0.4] <- "RG+"
   df$TYPE[df[ ,c("IP")] > 1.3 & df[ ,c("RG")] > 0.4] <- "IP+ RG+"
   df$TYPE[df[ ,c("IP")] > 1.3 & df[ ,c("Neuron")] > 1] <- "Neuron+ IP+"
   df$TYPE[df[ ,c("Neuron")] > 1 & df[ ,c("RG")] > 0.4] <- "Neuron+ RG+"

   df$TYPE <- factor(df$TYPE, levels = c("Neuron+", "IP+", "RG+"
     , "IP+ RG+", "Neuron+ IP+", "Neuron+ RG+")
   )
   df$CLUSTER <- factor(as.character(df$CLUSTER)
     , levels = sort(unique(as.numeric(as.character(df$CLUSTER)))))

   return(df)
 }

expr_flag_DF <- Positive_Negative_Expression_Flag(
  exDF = expr_flag_DF
  , highThreshold = 0.5
  , lowThreshold = 0.25
)

df <- Mixed_CellTypeEnriched_By_Cluster_Percent(
  exM = noCentExM
  , seuratO = centSO
  , highThreshold = 0.5
  , lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , fold_change = 0.4
  , cluster_order = NULL)


# Paper - plot with reordered clusters
gg <- Mixed_CellTypeEnriched_By_Cluster_Percent_Barplot(
  exM = noCentExM
  , seuratO = centSO
  , fold_change = 0.4
  , highThreshold = 0.75, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  , title = paste0(graphCodeTitle
    , "\n\nPercent of cells passing combinations of cell type enriched gene expression filters"
    , "\nKeep CC"
    , "\n+ = > 0.5 normalized expression"
    , "\n- = < 0.25 normalized expression"
  )
)
gg + ggplot_set_theme_publication
# Save
ggsave(paste0(outGraph
  , "PercentMixedCluster_CellTypeEnriched_Barplot_paper_supp.pdf")
  , width = 5, height = 3)

# Paper - table of percents for making differentiation diagram
percentDF <- Mixed_CellTypeEnriched_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , fold_change = 0.4
  , highThreshold = 0.5, lowThreshold = 0.25
  , type_keep = c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+", "RG+", "IP+"
    , "Neuron+")
  , cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
percentDF <- percentDF[with(percentDF, CLUSTER == 10), ]
percentDF$Percent_NA_Removed <-
  c(NA, (percentDF$value[2:6] / sum(percentDF$value[2:6])) * 100, NA)
percentDF[ ,c(3:4)] <- round(percentDF[ ,c(3:4)], 1)
write.csv(percentDF, file = paste0(
    outTable, "PercentMixedCluster_CellTypeEnriched_paper_supp.csv")
  , quote = FALSE)

################################################################################

### Number of mixed marker cells by region

Calculate_Number_Mixed_Marker_Cells_By_Region <- function(){
  print("Calculate_Number_Mixed_Marker_Cells_By_Region")
  # Flag mixed marker cells
  mixed_marker_cells_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO
  )
  mixed_marker_cells_DF <- Positive_Negative_Expression_Flag(
    exDF = mixed_marker_cells_DF, highThreshold = 0.5, lowThreshold = 0.5
  )
  # Format
  idx <- mixed_marker_cells_DF$TYPE %in%
    c("IP+ RG+", "Neuron+ RG+", "Neuron+ IP+")
  mixed_marker_cells_DF <- mixed_marker_cells_DF[idx, ]
  mixed_marker_cells_DF$TYPE <- droplevels(mixed_marker_cells_DF$TYPE)
  # Add region
  mixed_marker_cells_DF$Region <- "GZ"
  mixed_marker_cells_DF$Region[
    row.names(mixed_marker_cells_DF) %in% centSO@meta.data$CELL[centSO@meta.data$REGION == "CP"]
    ] <- "CP"
  # Caculate number
  mixed_marker_cells_DF <- aggregate(TYPE~CLUSTER, mixed_marker_cells_DF, table)
  mixed_marker_cells_DF <- melt(data.frame(Region = mixed_marker_cells_DF[ ,1], mixed_marker_cells_DF[ ,2]))
  # Format
  colnames(mixed_marker_cells_DF) <- c("Region", "Markers", "Number")
  mixed_marker_cells_DF$Markers <-
    gsub("\\.+", "\\+", mixed_marker_cells_DF$Markers)
  return(mixed_marker_cells_DF)
}
mixed_marker_cells_DF <- Calculate_Number_Mixed_Marker_Cells_By_Region()
write.csv(mixed_marker_cells_DF
  , file = paste0(outTable, "Number_Mixex_marker_Cells_By_Region.csv")
  , quote = FALSE, row.names = FALSE
)

################################################################################

### Mixed color tSNE

Mixed_tSNE_Format_Data_For_GGplot <- function(green_genes, red_genes){

  genes1 <- green_genes
  genes2 <- red_genes

  # Collect tSNE and expression values
  tsneDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  exp1_DF <- Mean_Expression(tsneDF, genes = genes1, exM = noCentExM)
  tsneDF$Expression_1 <- exp1_DF$EXPRESSION
  exp2_DF <- Mean_Expression(tsneDF, genes = genes2, exM = noCentExM)
  tsneDF$Expression_2 <- exp2_DF$EXPRESSION

  # Expression limits
  tsneDF$Expression_1[tsneDF$Expression_1 > 1] <- 1
  tsneDF$Expression_1[tsneDF$Expression_1 < 0] <- 0
  tsneDF$Expression_2[tsneDF$Expression_2 > 1] <- 1
  tsneDF$Expression_2[tsneDF$Expression_2 < 0] <- 0

  # Flag cells with mixed markers (>0.5 expression)
  tsneDF$Mixed <- "No"
  tsneDF$Mixed[tsneDF$Expression_1 > 0.5 & tsneDF$Expression_2 > 0.5] <- "Yes"

  # Transform expression to 1-255 range for rgb function
  tsneDF$Expression_1 <- tsneDF$Expression_1/1
  tsneDF$Expression_2 <- tsneDF$Expression_2/1
  tsneDF$Expression_1 <- round(tsneDF$Expression_1 * 255, 0)
  tsneDF$Expression_2 <- round(tsneDF$Expression_2 * 255, 0)

  # scale <- rowSums(tsneDF[c("Expression_1", "Expression_2")])
  # scale <- scale[scale < 0.25] <- 0.25
  # scale <- round(scale * 1, 0)
  # tsneDF$Size <- scale

  # Convert to rgb
  tsneDF <- within(tsneDF
    , mix <- rgb(green = Expression_1, red = Expression_2, blue = 0
      , maxColorValue = 255
      # , alpha = scale
    )
  )
  return(tsneDF)
}

# RG IP
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "RG"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "IP"]
)
gg1 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  geom_point(data = tsneDF[tsneDF$Mixed == "Yes", ]
    , aes(x = tSNE_1, y = tSNE_2, color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = RG; Red = IP") +
  ggplot_set_theme_publication_nolabels
# RG Neuron
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "RG"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "Neuron"]
)
gg2 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  geom_point(data = tsneDF[tsneDF$Mixed == "Yes", ]
    , aes(x = tSNE_1, y = tSNE_2, color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = RG; Red = Neuron") +
  ggplot_set_theme_publication_nolabels
# IP Neuron
tsneDF <- Mixed_tSNE_Format_Data_For_GGplot(
  green_genes = kmDF$Gene.Symbol[kmDF$Grouping == "IP"]
  , red_genes = kmDF$Gene.Symbol[kmDF$Grouping == "Neuron"]
)
gg3 <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = mix), size = 0.2) +
  geom_point(data = tsneDF[tsneDF$Mixed == "Yes", ]
    , aes(x = tSNE_1, y = tSNE_2, color = mix), size = 0.2) +
  scale_color_identity() +
  ggtitle("Green = IP; Red = Neuron") +
  ggplot_set_theme_publication_nolabels
# Combine
Plot_Grid(list(gg1, gg2, gg3), ncol = 3, rel_height = 0.2
  , title = paste0(graphCodeTitle
    , "\n\nTwo channel mixed expression of marker genes"
  )
)
ggsave(paste0(outGraph, "tSNE_Mixed.png")
  , width = 12, height = 5)
################################################################################

### Number / Percent of cells in CC phase subset by markers
print("### Number / Percent of cells in CC phase subset by markers")

## Percent of cells in CC phase
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)

df <- rbind(
  data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$oRG > 0.5]), SUBSET = "vRG+ oRG+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$IP > 0.5]), SUBSET = "vRG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$oRG > 0.5 & df1$IP > 0.5]), SUBSET = "oRG+ IP+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$vRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "vRG+ Neuron+")
  , data.frame(Percent_Of_Table(df1$PHASE[df1$oRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "oRG+ Neuron+")
)
df$tableArg <- factor(df$tableArg, levels = c("G1", "S", "G2M"))
# Plot
ggplot(df, aes(x = SUBSET, y = Freq, fill = tableArg)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_discrete(name = "CC phase") +
  xlab("Cell subset") +
  ylab("Percent of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nPercent of cell subsets in cell cycle phases"
    , "\nSubset by mean expression levels of groups of marker genes"
    , "\n(> 0.5 log normalized mean expression)")
  )
ggsave(paste0(outGraph, "PercentCCphase_RG_IP_Subset_Barplot.pdf")
  , width = 7, height = 7)


## Number of cells in CC phase

df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)

df <- rbind(
  data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$IP > 0.5]), SUBSET = "RG+ IP+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5]), SUBSET = "IP+")
  , data.frame(table(df1$PHASE[df1$RG > 0.5 & df1$Neuron > 0.5]), SUBSET = "RG+ Neuron+")
  , data.frame(table(df1$PHASE[df1$IP > 0.5 & df1$Neuron > 0.5]), SUBSET = "IP+ Neuron+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$oRG > 0.5]), SUBSET = "vRG+ oRG+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$IP > 0.5]), SUBSET = "vRG+ IP+")
  , data.frame(table(df1$PHASE[df1$oRG > 0.5 & df1$IP > 0.5]), SUBSET = "oRG+ IP+")
  , data.frame(table(df1$PHASE[df1$vRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "vRG+ Neuron+")
  , data.frame(table(df1$PHASE[df1$oRG > 0.5 & df1$Neuron > 0.5]), SUBSET = "oRG+ Neuron+")
)
df$Var1 <- factor(df$Var1, levels = c("G1", "S", "G2M"))
# Plot
ggplot(df, aes(x = SUBSET, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_discrete(name = "CC phase") +
  xlab("Cell subset") +
  ylab("Number of cells") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cell subsets in cell cycle phases"
    , "\nSubset by mean expression levels of groups of marker genes"
    , "\n(> 0.5 log normalized mean expression)")
  )
ggsave(paste0(outGraph, "NumberCCphase_RG_IP_Subset_Barplot.pdf")
  , width = 7, height = 7)
################################################################################

me_DF <- ME_CellType_L[["RG_to_IP_8"]]
me_DF <- me_DF[me_DF$Cell_Subset_025 == "RG IP" & me_DF$Fold_Change_Cutoff == 0.25, ]
df1 <- me_DF[me_DF$variable == "IP eigengene", ]
df2 <- me_DF[me_DF$variable == "RG eigengene", ]
df3 <- merge(df1, df2, by = "CellID")
ggplot(df3, aes(x = value.x, y = value.y)) +
  geom_point() +
  xlab(df3$variable.x[1]) +
  ylab(df3$variable.y[1]) +
  ggtitle("S phase")
ggsave(paste0(outGraph, "Eigengene_Sphase_ScatterPlot.png"))

me_DF <- ME_CellType_L[["RG_to_IP_10"]]
me_DF <- me_DF[me_DF$Cell_Subset_025 == "RG IP" & me_DF$Fold_Change_Cutoff == 0.25, ]
df1 <- me_DF[me_DF$variable == "IP eigengene", ]
df2 <- me_DF[me_DF$variable == "RG eigengene", ]
df3 <- merge(df1, df2, by = "CellID")
ggplot(df3, aes(x = value.x, y = value.y)) +
  geom_point() +
  xlab(df3$variable.x[1]) +
  ylab(df3$variable.y[1]) +
  ggtitle("G2/M phase")
ggsave(paste0(outGraph, "Eigengene_G2Mphase_ScatterPlot.png"))
################################################################################

ssExM <- noCentExM[rownames(noCentExM) %in% c("STMN2", "PAX6", "PCNA"), ]
ssExM <- ssExM > 0.5
ssExM <- as.data.frame(t(ssExM))
sum(ssExM$PAX6 == TRUE & ssExM$STMN2 == TRUE & ssExM$PCNA == TRUE)
sum(ssExM$PAX6 == TRUE & ssExM$STMN2 == TRUE & ssExM$PCNA == FALSE)
sum(ssExM$PAX6 == TRUE & ssExM$STMN2 == FALSE & ssExM$PCNA == TRUE)
sum(ssExM$PAX6 == TRUE)
sum(ssExM$STMN2 == TRUE)
sum(ssExM$PCNA == TRUE)
