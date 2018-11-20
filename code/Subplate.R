# Damon Polioudakis
# 2018-06-14
# Plot known marker lists as heatmaps, violin plots, feature plots

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

options(stringsAsFactors = FALSE)

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(reshape2)
require(gridExtra)
require(ggplot2)
require(cowplot)
require(viridis)
require(WGCNA)
require(ggpubr)
source("Function_Library.R")
source("GGplot_Theme.R")

## Inputs

# Seurat cluster round 2
load("../analysis/analyzed_data/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_Cluster3_seuratO.Robj")

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

# biomaRt gene info
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Miller LCM from Allen
# Row annotations
millerLCM_row_annot_DF <- read.csv(
  "../allen_brain_data/miller_LCM/lmd_matrix_12840/rows_metadata.csv")
# Processed data from Luis
# [26] "layer.exp"
# [27] "layer.exp.sc"
# [28] "layer.max"
load("../allen_brain_data/miller_LCM/fetal_Ctx_layer_exp.Rdata")
idx <- match(rownames(layer.max), millerLCM_row_annot_DF$gene_symbol)
idx <- idx[! is.na(idx)]
millerLCM_row_annot_DF <- millerLCM_row_annot_DF[idx, ]

## Output Directories
out_graph <- "../analysis/graphs/Subplate/Subplate_"
out_table <- "../analysis/tables/Subplate/Subplate_"
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)

## Other variables
script_name <- "Subplate.R"
################################################################################

### Functions

Format_Miller_LCM_By_Region_For_ggPlot <- function(genes){
  print("Format_Miller_LCM_By_Region_For_ggPlot")
  # browser()
  # Subset to genes of interest
  genes <- genes[! duplicated(genes)]
  ss_row_annot_DF <- millerLCM_row_annot_DF[
    millerLCM_row_annot_DF$gene_symbol %in% genes, ]
  idx <- match(rownames(layer.max), millerLCM_row_annot_DF$gene_symbol)
  millerLCM_row_annot_DF <- millerLCM_row_annot_DF[idx, ]
  ex_M <- layer.max[
    rownames(layer.max) %in% ss_row_annot_DF$gene_symbol, ]
  # Select probe with highest expression for each gene
  ex_M <- t(ex_M)
  ex_M <- melt(ex_M)
  idx <- match(ex_M$Var2, millerLCM_row_annot_DF$gene_symbol)
  ex_M$gene_symbol <- millerLCM_row_annot_DF$gene_symbol[idx]
  ex_M$gene_symbol <- factor(ex_M$gene_symbol, levels = rev(genes))
  ex_M <- aggregate(value~gene_symbol+Var1, ex_M, max)
  ex_M <- ex_M[! duplicated(ex_M[ ,c("Var1", "gene_symbol")]), ]
  ex_M <- ex_M[! ex_M$Var1 == "", ]
  ex_M <- dcast(ex_M, Var1~gene_symbol, value.var = "value")
  rownames(ex_M) <- ex_M[ ,1]
  ex_M <- ex_M[ ,-1]
  # Center and scale
  ex_M <- scale(ex_M)
  ex_M <- melt(ex_M)
  # Set limits to -1.5 1.5
  ex_M$value[ex_M$value < -1.5] <- -1.5
  ex_M$value[ex_M$value > 1.5] <- 1.5
  # Add area
  ex_M$Area <- substr(ex_M$Var1, 1, 1)
  # Add zone
  ex_M$Zone <- substr(ex_M$Var1, 2, 3)
  ex_M$Zone <- factor(ex_M$Zone, levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ"))
  return(ex_M)
}
################################################################################

### Subplate marker expression

Plot_Subplate_Marker_Expression <- function(){
  print("Plot_Subplate_Marker_Expression")

  ## Miller LCM

  ## Full list
  # By area and zone
  gg_DF <- Format_Miller_LCM_By_Region_For_ggPlot(
    genes = c(as.character(kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate"])
      , c("NXPH4", "INPP4B", "CTGF", "HTR1D", "TPD52L1", "ADRA2A"))
    )
  # By zone
  gg_DF <- aggregate(value~Zone+Var2, gg_DF, mean)
  # Plot
  ggplot(gg_DF, aes(x = Var2, y = Zone, fill = value)) +
    geom_tile() +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
        , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
    ylab("Region") +
    xlab("Gene") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # theme(axis.text.x = element_blank()) +
    ggtitle(paste0(script_name
      , "\n\nSubplate markers (Molnar) expression in miller LCM"
      , "\nExpression mean centered and variance scaled across all samples"))
  ggsave(paste0(out_graph, "SPmolnarMarks_Miller_Heatmap.pdf")
    , width = 7, height = 5)

  ## Refined list
  # By area and zone
  gg_DF <- Format_Miller_LCM_By_Region_For_ggPlot(
    genes = c("CDH10", "CDH18", "GABRA5", "NR4A2", "SV2B", "CPLX3")
  )
  # By zone
  gg_DF <- aggregate(value~Zone+Var2, gg_DF, mean)
  # Plot
  ggplot(gg_DF, aes(x = Var2, y = Zone, fill = value)) +
    geom_tile() +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
        , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
    ylab("Region") +
    xlab("Gene") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # theme(axis.text.x = element_blank()) +
    ggtitle(paste0(script_name
      , "\n\nSubplate markers (Molnar) expression in miller LCM"
      , "\nExpression mean centered and variance scaled across all samples"))
  ggsave(paste0(out_graph, "SPmolnarMarksRefined_Miller_Heatmap.pdf")
    , width = 7, height = 5)

  ## Miller LCM derived subplate markers
  # By area and zone
  gg_DF <- Format_Miller_LCM_By_Region_For_ggPlot(
    genes = c("CDH18", "HS3ST3B1", "HAS3", "HCRTR2", "CTXN3", "EPHA8", "SERTM1"
      , "MGC12916", "DKK1", "PTGS2")
  )
  # By zone
  gg_DF <- aggregate(value~Zone+Var2, gg_DF, mean)
  # Plot
  ggplot(gg_DF, aes(x = Var2, y = Zone, fill = value)) +
    geom_tile() +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
        , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
    ylab("Region") +
    xlab("Gene") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # theme(axis.text.x = element_blank()) +
    ggtitle(paste0(script_name
      , "\n\nExpression of subplate markers derived from Miller LCM in Miller LCM"
      , "\nExpression mean centered and variance scaled across all samples"))
  ggsave(paste0(out_graph, "SPmillerLCMmarks_Miller_Heatmap.pdf")
    , width = 7, height = 5)


  ## Feature plots

  # Feature plot normalized, mean centered scaled
  # Mean expression
  ggL <- FeaturePlot(
    genes = c("CDH10", "CDH18", "GABRA5", "NR4A2", "SV2B", "CPLX3")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = ""
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.4)
    , title = paste0(script_name
      , "\n\nExpression of refined subplate marker list"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "SPmarksRefined_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 4)

  # Feature plot normalized, mean centered scaled
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("CDH10", "CDH18", "GABRA5", "NR4A2", "SV2B", "CPLX3")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.3)
    , title = paste0(script_name
      , "\n\nExpression of refined subplate marker list"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "SPmarksRefined_FeaturePlotIndividual_NormalizedCenteredScaled.png")
    , width = 12, height = 12)

  ## Subplate markers obtained from Miller LCM DE of SP vs CP, MZ, VZ, SVZ
  # Feature plot normalized, mean centered scaled
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("CDH18", "HS3ST3B1", "HAS3", "HCRTR2", "CTXN3", "EPHA8", "SERTM1"
      , "MGC12916", "DKK1", "PTGS2")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.1)
    , title = paste0(script_name
      , "\n\nExpression of subplate markers derived from Miller LCM"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "SPmarksMillerLCM_FeaturePlotIndividual_NormalizedCenteredScaled.png")
    , width = 12, height = 12)
  # Feature plot normalized
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("CDH18", "HS3ST3B1", "HAS3", "HCRTR2", "CTXN3", "EPHA8", "SERTM1"
      , "MGC12916", "DKK1", "PTGS2")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -0.5
    , limHigh = 2
    , centScale = FALSE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.1)
    , title = paste0(script_name
      , "\n\nExpression of subplate markers derived from Miller LCM"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "SPmarksMillerLCM_FeaturePlotIndividual_Normalized.png")
    , width = 12, height = 12)

  # ST18
  # Feature plot normalized, mean centered scaled
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("ST18")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.3)
    , title = paste0(script_name
      , "\n\nExpression of ST18"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "ST18_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 4)

  # Markers and ST18
  # Feature plot normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = c(as.character(kmDF$Gene.Symbol)[
      kmDF$Grouping %in% c("Subplate Miller LCM"
        , "Excitatory Deep Layer Cortical")]
      , "ST18")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
    , geneGrouping = c(as.character(kmDF$Grouping)[
      kmDF$Grouping %in% c("Subplate Miller LCM"
        , "Excitatory Deep Layer Cortical")]
      , "ST18")
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.3)
    , title = paste0(script_name
      , "\n\nExpression of markers and ST18"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "Marks_ST18_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 7)

  # Sub-clustering feature plot of markers and ST18
  # Feature plot normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = c(as.character(kmDF$Gene.Symbol)[
      kmDF$Grouping %in% c("Subplate Miller LCM"
        , "Excitatory Deep Layer Cortical")]
      , "ST18")
    , tsneDF = as.data.frame(so@dr$tsne@cell.embeddings)
    , seuratO = so
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
    , geneGrouping = c(as.character(kmDF$Grouping)[
      kmDF$Grouping %in% c("Subplate Miller LCM"
        , "Excitatory Deep Layer Cortical")]
      , "ST18")
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.3)
    , title = paste0(script_name
      , "\n\nExpression of markers and ST18"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph
    , "Marks_ST18_Subcluster_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 7)

  ## Temporal CP top 10 DE genes
  # Feature plot normalized, mean centered scaled
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("ERVMER34-1", "STK31", "MET", "STK31", "ANLN", "NR4A2", "LMO1"
    , "BHLHE40", "MUM1L1", "NR2F2", "GRIK1", "SPINT1")
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.1)
    , title = paste0(script_name
      , "\n\nExpression of temporal CP DE genes obtained from Miller LCM"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "TemporalCPtop10DE_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 14)

  ## Temporal CP top 20 DE genes
  # Feature plot normalized, mean centered scaled
  # Individual expression
  ggL <- FeaturePlot(
    genes = c("ERVMER34-1", "NPY", "STK31", "MET", "NTNG1", "ANLN", "NR4A2"
    , "LMO1", "BHLHE40", "SLC17A6", "FAM150B", "ASCL1", "MUM1L1", "CXCL12"
    , "EXD1", "ITGA2", "GRIK1", "FREM2", "NR2F2", "NTNG2", "SPINT1"
    , "CHRNA3"
    , "PLP1"
    , "MRAP2"
    , "IQCJ"
    , "IRX3"
    , "TMEM249"
    , "CLUL1"
    , "TMEM61"
    )
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r', rel_height = c(0.1)
    , title = paste0(script_name
      , "\n\nExpression of temporal CP DE genes obtained from Miller LCM"
      , "\nNormalized expression, mean centered and variance scaled"
      , "\n")
  )
  ggsave(paste0(out_graph, "TemporalCPtop30DE_FeaturePlot_NormalizedCenteredScaled.png")
    , width = 12, height = 34)
}
################################################################################

### Plot ST18

Plot_ST18_Expression <- function(){
  print("Plot_ST18_Expression")
  # Feature plot
  # Normalized, no mean centering scaling
  ggL <- FeaturePlot(
    genes = "ST18"
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = noCentExM
    , limLow = -0.5
    , limHigh = 2
    , geneGrouping = ""
    , centScale = FALSE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(script_name
    , "\n\nExpression of ST18"
    , "\nNormalized expression"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
  ggsave(paste0(out_graph, "ST18_FeaturePlot_Normalized.png")
    , width = 12, height = 4)

  # Feature plot
  # Normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = "ST18"
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = ""
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(script_name
    , "\n\nExpression of ST18"
    , "\nCluster round 2 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
  ggsave(paste0(out_graph
    , "ST18_FeaturePlot_Round2NormalizedCenteredScaled.png")
    , width = 12, height = 4)

  # Feature plot
  # Cluster round 1 normalized, mean centered scaled
  ggL <- FeaturePlot(
    genes = "ST18"
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = rd1CentExM
    , limLow = -1.5
    , limHigh = 1.5
    , geneGrouping = ""
    , centScale = TRUE
    , size = (400/nrow(centSO@scale.data))^(1/3)
  )
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
  # now add the title
  title <- ggdraw() + draw_label(paste0(script_name
    , "\n\nExpression of ST18"
    , "\nCluster round 1 normalized expression, mean centered and variance scaled"
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.4, 1))
  ggsave(paste0(out_graph, "ST18_FeaturePlot_Round1NormalizedCenteredScaled.png")
    , width = 12, height = 4)
}
################################################################################

### Genes that correlate to ST18 expression

# Correlation to ST18
Correlation_to_ST18 <- function(){
  print("Correlation_to_ST18")
  cor_ST18_DF <- cor(t(noCentExM)
    , noCentExM[rownames(noCentExM) %in% "ST18", ])
  cor_ST18_DF <- as.data.frame(cor_ST18_DF)
  colnames(cor_ST18_DF) <- "Correlation"
  cor_ST18_DF$Gene <- rownames(cor_ST18_DF)

  cor_ST18_DF <- cor_ST18_DF[order(cor_ST18_DF[,"Correlation"]), ]
  cor_ST18_DF <- cor_ST18_DF[! is.na(cor_ST18_DF[,"Correlation"]), ]

  sp_known_marks <- kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate"]
  print(mean(cor_ST18_DF$Correlation[cor_ST18_DF$Gene %in% sp_known_marks]))

  sp_miller_marks <- kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate Miller LCM"]
  print(mean(cor_ST18_DF$Correlation[cor_ST18_DF$Gene %in% sp_miller_marks]))

  return(cor_ST18_DF)
}
# Correlation_to_ST18 <- function(){
#   print("Correlation_to_ST18")
#   cor_ST18_DFL <- lapply(rownames(noCentExM), function(gene){
#     cor_O <- cor.test(noCentExM[rownames(noCentExM) %in% "ST18", ]
#     , noCentExM[rownames(noCentExM) %in% gene])
#     cor_ST18_DF <- data.frame(
#       Gene = gene
#       , Correlation = cor_O$estimate
#       , Pvalue = cor_O$p.value)
#     return(cor_ST18_DF)
#   })
#   cor_ST18_DF <- do.call("rbind", cor_ST18_DFL)
#   cor_ST18_DF <- cor_ST18_DF[order(cor_ST18_DF[,2]), ]
#   cor_ST18_DF <- cor_ST18_DF[! is.na(cor_ST18_DF[,2]), ]
#
#   sp_known_marks <- kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate"]
#   print(mean(cor_ST18_DF$Correlation[cor_ST18_DF$Gene %in% sp_known_marks]))
#
#   sp_miller_marks <- kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate Miller LCM"]
#   print(mean(cor_ST18_DF$Correlation[cor_ST18_DF$Gene %in% sp_miller_marks]))
#
#   return(cor_ST18_DF)
# }
################################################################################

### Cells expressing SP markers and ST18

## Violin plots to check expression levels

Plot_Expression_Levels_Violin <- function(){
  print("Plot_Expression_Levels_Violin")

  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate"]
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping[kmDF$Grouping %in% "Subplate"]
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(out_graph, "SPmarks_ViolinPlot_Normalized.png")
    , width = 9, height = 5)
  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = kmDF$Gene.Symbol[kmDF$Grouping %in% "Subplate"]
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = kmDF$Grouping[kmDF$Grouping %in% "Subplate"]
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nMarker gene expression by cluster"
      , "\nNormalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(out_graph, "SPmarks_ViolinPlot_NormalizedCenteredScaled.png")
    , width = 9, height = 5)

  # Subplate Miller LCM obtained markers
  subplate_marks <- c("CDH18", "HS3ST3B1", "HAS3", "HCRTR2", "CTXN3", "EPHA8"
    , "SERTM1", "MGC12916", "DKK1", "PTGS2")
  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = subplate_marks
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = rep("Subplate", length(subplate_marks))
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nSubplate marker genes obtained from Miller LCM DE expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(out_graph, "SPmillerLCMmarks_ViolinPlot_Normalized.png")
    , width = 9, height = 5)
  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = subplate_marks
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = rep("Subplate", length(subplate_marks))
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nSubplate marker genes obtained from Miller LCM DE expression by cluster"
      , "\nNormalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(out_graph
    , "SPmillerLCMmarks_ViolinPlot_NormalizedCenteredScaled.png")
    , width = 9, height = 5)

  # ST18
  # Normalized
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = "ST18"
    , exprM = noCentExM
    , clusterIDs = centSO@ident
    , grouping = ""
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nST18 expression by cluster"
      , "\nNormalized expression"
      , "\n"))
  ggsave(paste0(out_graph, "ST18_ViolinPlot_Normalized.png")
    , width = 9, height = 5)
  # Normalized centered scaled
  Gene_Expression_By_Cluster_ViolinPlot(
    genes = "ST18"
    , exprM = as.matrix(centSO@scale.data)
    , clusterIDs = centSO@ident
    , grouping = ""
  ) +
    ggtitle(paste0(script_name
      , "\n"
      , "\nST18 expression by cluster"
      , "\nNormalized mean centered scaled expression"
      , "\n"))
  ggsave(paste0(out_graph, "ST18_ViolinPlot_NormalizedCenteredScaled.png")
    , width = 9, height = 5)
}

Plot_Intersection_ST18_SPmarks <- function(){
  print("Plot_Intersection_ST18_SPmarks")
  ## Cells expressing both ST18 and SP markers

  sp_known_marks <- c("CDH18", "HS3ST3B1", "HAS3", "HCRTR2", "CTXN3", "EPHA8"
    , "SERTM1", "MGC12916", "DKK1", "PTGS2")

  # Cells expressing ST18
  ex_sp_M <- noCentExM[rownames(noCentExM) %in% sp_known_marks, ]
  table(colMeans(ex_sp_M) > 0.1)
  sp_marks_cell_IDs <- colnames(ex_sp_M)[colMeans(ex_sp_M) > 0.1]
  ex_ST18 <- noCentExM[rownames(noCentExM) %in% "ST18", ]
  table(ex_ST18 > 0.5)
  st18_cell_IDs <- names(ex_ST18)[ex_ST18 > 0.5]

  cell_cat_DF <- data.frame(Cell = names(centSO@ident), Cluster = centSO@ident)
  cell_cat_DF$Subplate <- cell_cat_DF$Cell %in% sp_marks_cell_IDs
  cell_cat_DF$ST18 <- cell_cat_DF$Cell %in% st18_cell_IDs
  cell_cat_DF$Subplate_Cluster <-
    cell_cat_DF$Cell %in% names(so@ident)[so@ident == 2]
  cell_cat_DF$Subplate_ST18 <- ifelse(
    cell_cat_DF$Subplate == TRUE & cell_cat_DF$ST18 == TRUE, "ST18+ & SP marker +"
      , ifelse(cell_cat_DF$ST18 == TRUE, "ST18+"
        , ifelse(cell_cat_DF$Subplate == TRUE, "SP marker +", NA)))

  gg_DF <- cell_cat_DF[! is.na(cell_cat_DF$Subplate_ST18), ]
  gg_DF$Cluster <- factor(gg_DF$Cluster
    , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
  ggplot(gg_DF, aes(x = Subplate_ST18, fill = Subplate_Cluster)) +
    geom_bar() +
    ggplot_set_theme_publication +
    ggtitle(paste0(script_name
      , "\n\nNumber of cells expressing ST18, subplate markers, or both"))
  ggsave(paste0(out_graph, "express_ST18_SPmarks_barplot.pdf")
    , width = 5, height = 4)

  gg_DF <- cell_cat_DF[! is.na(cell_cat_DF$Subplate_ST18), ]
  gg_DF$Cluster <- factor(gg_DF$Cluster
    , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
  ggplot(gg_DF, aes(x = Cluster, fill = Subplate_ST18)) +
    geom_bar() +
    ggplot_set_theme_publication +
    ggtitle(paste0(script_name
      , "\n\nNumber of cells expressing ST18, subplate markers, or both"))
  ggsave(paste0(out_graph, "express_ST18_SPmarks_by_cluster_barplot.pdf")
    , width = 5, height = 4)

  gg_DF <- data.frame(
    Number = c(
      sum(cell_cat_DF$ST18 == TRUE)
      , sum(cell_cat_DF$Subplate == TRUE))
    , Both = c(sum(cell_cat_DF$ST18 == TRUE & cell_cat_DF$Subplate == TRUE)
        , sum(cell_cat_DF$ST18 == TRUE & cell_cat_DF$Subplate == TRUE))
    , Genes = c("ST18", "Subplate"))
  gg_DF <- melt(gg_DF)

  ggplot(gg_DF, aes(x = Genes, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    ggplot_set_theme_publication +
    ggtitle(paste0(script_name
      , "\n\nNumber of cells expressing ST18, subplate markers, or both"))
  ggsave(paste0(out_graph, "express_ST18_SPmarks_barplot2.pdf")
    , width = 5, height = 4)
}
################################################################################

### DE

# DE subplate cluster versus all cells

DE_By_CellIDs <- function(
  cell_IDs_1, cell_IDs_2, positive_DE_label = NA, negative_DE_label = NA) {
  ids <- c(cell_IDs_1, cell_IDs_2)
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% ids]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% ids
    , c("nUMI", "librarylab", "individual", "CELL")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$groups <- "cell_IDs_1"
  termsDF$groups[termsDF$CELL %in% cell_IDs_2] <- "cell_IDs_2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ groups+nUMI+librarylab+individual")
  # Format LM DE
  deDF <- data.frame(
    Log2_FC_Group1_vs_Group2 = deLM$coefmat[ ,"groupscell_IDs_2"]
    , Pvalue = deLM$pvalmat[ ,"groupscell_IDs_2"]
  )
  deDF$Gene = row.names(deDF)
  # Make DE cell_IDs_1 positive fold change
  deDF$Log2_FC_Group1_vs_Group2 <- deDF$Log2_FC_Group1_vs_Group2 * -1
  # Label DE group
  deDF$DE_Group <- NA
  deDF$DE_Group[deDF$Log2_FC_Group1_vs_Group2 > 0] <- positive_DE_label
  deDF$DE_Group[deDF$Log2_FC_Group1_vs_Group2 < 0] <- negative_DE_label
  # Order by fold change
  deDF <- deDF[order(deDF$Log2_FC_Group1_vs_Group2), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method = "BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))
  return(deDF)
}

DE_SPcluster_Versus_Allcells <- function(){
  print("DE_SPcluster_Versus_Allcells")
  # SP cluster cells versus all cells
  sp_cluster_ids <- names(so@ident)[so@ident %in% 2]
  all_other_ids <- names(centSO@ident)[! names(centSO@ident) %in% sp_cluster_ids]
  de_sp_vs_all_DF <- DE_By_CellIDs(
    cell_IDs_1 = sp_cluster_ids, cell_IDs_2 = all_other_ids)
  return(de_sp_vs_all_DF)
}

DE_SP_Cluster_Versus_Other_Deep_Layer_Cluster_Cells <- function(){
  print("DE_SP_Cluster_Versus_Other_Deep_Layer_Cluster_Cells")
  # SP cluster cells versus deep layer cluster cells
  # SP cells: subcluster 3:2 and cluster 13
  # sp_cluster_ids <- names(so@ident)[so@ident %in% 2]
  sp_cluster_ids <- c(
    names(so@ident)[so@ident %in% c(2)]
    , names(centSO@ident)[centSO@ident %in% 13]
  )
  all_other_cluster3_ids <- names(so@ident)[! names(so@ident) %in% sp_cluster_ids]
  de_sp_vs_cluster3_DF <- DE_By_CellIDs(
    cell_IDs_1 = sp_cluster_ids, cell_IDs_2 = all_other_cluster3_ids)
  return(de_sp_vs_cluster3_DF)
}

# Compile gene lists
Compile_Gene_Lists <- function(
  de_sp_vs_all_DF, de_sp_vs_cluster3_DF, cor_ST18_DF){
  print("Compile_Gene_Lists")
  gene_group_DF <- rbind(
    data.frame(Gene = de_sp_vs_all_DF$Gene[
      de_sp_vs_all_DF$Log2_FC_Group1_vs_Group2 > 0.4 &
      de_sp_vs_all_DF$FDR < 0.05]
      , Group = "SP_vs_all_DE04")
    , data.frame(Gene = tail(de_sp_vs_all_DF$Gene[
      de_sp_vs_all_DF$Log2_FC_Group1_vs_Group2 > 0.4 &
      de_sp_vs_all_DF$FDR < 0.05], 100)
      , Group = "SP_vs_all_DE04_top100")
    , data.frame(Gene = tail(de_sp_vs_all_DF$Gene[
      de_sp_vs_all_DF$Log2_FC_Group1_vs_Group2 > 0.4 &
      de_sp_vs_all_DF$FDR < 0.05], 80)
      , Group = "SP_vs_all_DE04_top80")
    , data.frame(Gene = de_sp_vs_cluster3_DF$Gene[
      de_sp_vs_cluster3_DF$Log2_FC_Group1_vs_Group2 > 0.4 &
      de_sp_vs_cluster3_DF$FDR < 0.5]
      , Group = "SP_vs_cluster3_DE04")
    # , data.frame(Gene = cor_ST18_DF$Gene[cor_ST18_DF$Pvalue < 0.05]
    #   , Group = "ST18_correlate_pval05")
    # , data.frame(Gene = cor_ST18_DF$Gene[cor_ST18_DF$Pvalue < 0.01]
    #   , Group = "ST18_correlate_pval01")
    # , data.frame(Gene = cor_ST18_DF$Gene[cor_ST18_DF$Pvalue < 0.001]
    #   , Group = "ST18_correlate_pval001")
    , data.frame(Gene = tail(cor_ST18_DF$Gene, 100)
      , Group = "ST18_correlate_top100")
    , data.frame(Gene = tail(cor_ST18_DF$Gene, 80)
      , Group = "ST18_correlate_top80")
    , data.frame(Gene = tail(cor_ST18_DF$Gene, 20)
      , Group = "ST18_correlate_top20")
  )
  genes <- intersect(
    gene_group_DF$Gene[gene_group_DF$Group == "SP_vs_all_DE04_top80"]
    , gene_group_DF$Gene[gene_group_DF$Group == "ST18_correlate_top80"])
  gene_group_DF <- rbind(gene_group_DF,
     data.frame(Gene = genes, Group = "Intersect_SPDE_ST18cor"))
  return(gene_group_DF)
}

Output_Compiled_Gene_Lists <- function(){
  write.csv(gene_group_DF, file = paste0(out_table, "Marker_Lists.csv")
  , quote = FALSE, header = TRUE, row.names = FALSE)
}

Plot_DE_genes_Feature_Plot <- function(gene_group_DF){
  print("Plot_DE_genes_Feature_Plot")
  ## Feature plots of DE genes
  # Plot
  ggL <- FeaturePlot(
    genes = gene_group_DF$Gene
    , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
    , seuratO = centSO
    , exM = centSO@scale.data
    , limLow = -1.5, limHigh = 1.5
    , centScale = TRUE
    , geneGrouping = gene_group_DF$Group)
  ggL <- lapply(ggL, function(gg){
    gg <- gg + ggplot_set_theme_publication
    return(gg)
  })
  ggL[[1]] <- ggL[[1]] + theme(legend.position = "none")
  ggL[[2]] <- ggL[[2]] + theme(legend.position = "none")
  Plot_Grid(ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(script_name
      , "\n\nMean expression of subplate cluster DE genes")
    )
  ggsave(paste0(out_graph, "SPcluster_DE_genes_FeaturePlot.png"), height = 18)
}

Plot_SP_DE_Genes_In_Miller_LCM <- function(gene_group_DF){
  print("Plot_SP_DE_Genes_In_Miller_LCM")
  # Allen LCM
  # By area and zone
  gg_DFL <- lapply(split(gene_group_DF, gene_group_DF$Group)
    , function(ss_gene_group_DF){
    genes <- ss_gene_group_DF$Gene
    gg_DF <- Format_Miller_LCM_By_Region_For_ggPlot(genes = genes)
    return(gg_DF)
  })
  # Plot
  gg_L <- lapply(names(gg_DFL), function(name){
    gg_DF <- gg_DFL[[name]]
    ggplot(gg_DF, aes(x = Var2, y = Zone, fill = value)) +
      facet_wrap(~Area, ncol = 6) +
      geom_tile() +
      scale_fill_distiller(name = "Normalized\nexpression", type = "div"
          , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
      ylab("Region") +
      xlab("Gene") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # theme(axis.text.x = element_blank()) +
      ggtitle(name)
  })
  Plot_Grid(gg_L, ncol = 1, title = paste0(script_name
    , "\n\nMiller LCM expression"
    , "\nExpression mean centered and variance scaled across all samples"))
  ggsave(paste0(out_graph, "Miller_By_Area_Heatmap.pdf")
    , width = 16, height = 26, dpi = 150)

  # By zone
  gg_DFL <- lapply(gg_DFL, function(gg_DF){
    aggregate(value~Zone+Var2, gg_DF, mean)
  })
  # Plot
  gg_L <- lapply(names(gg_DFL), function(name){
    gg_DF <- gg_DFL[[name]]
    ggplot(gg_DF, aes(x = Var2, y = Zone, fill = value)) +
      geom_tile() +
      scale_fill_distiller(name = "Normalized\nexpression", type = "div"
          , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
      ylab("Region") +
      xlab("Gene") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # theme(axis.text.x = element_blank()) +
      ggtitle(name)
  })
  Plot_Grid(gg_L, ncol = 1, title = paste0(script_name
    , "\n\nMiller LCM expression"
    , "\nExpression mean centered and variance scaled across all samples"))
  ggsave(paste0(out_graph, "Miller_Heatmap.pdf")
    , width = 16, height = 26)
}

Calculate_Eigengene_SP_DE_Gene_Lists_In_MillerLCM <- function(gene_group_DF){
  print("Calculate_Eigengene_SP_DE_Gene_Lists_In_MillerLCM")
  # By area and zone
  me_DFL <- lapply(split(gene_group_DF, gene_group_DF$Group)
    , function(ss_gene_group_DF){
    genes <- ss_gene_group_DF$Gene
    me_DF <- Calculate_Eigengene_In_MillerLCM(genes = genes)
    return(me_DF)
  })
  return(me_DFL)
}
Calculate_Eigengene_In_MillerLCM <- function(genes){
  print("Calculate_Eigengene_In_MillerLCM")
  idx <- match(rownames(layer.max), millerLCM_row_annot_DF$gene_symbol)
  genes_DF <- millerLCM_row_annot_DF[idx, ]
  genes_DF$Class <- ifelse(genes_DF$gene_symbol %in% genes, "SP", "Non-SP")
  me_DF <- moduleEigengenes(t(layer.max), genes_DF$Class)$eigengenes
  me_DF$Sample <- colnames(layer.max)
  # Add area
  me_DF$Area <- substr(me_DF$Sample, 1, 1)
  # Add zone
  me_DF$Zone <- substr(me_DF$Sample, 2, 3)
  me_DF$Zone <- factor(me_DF$Zone
    , levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ"))
  return(me_DF)
}
Plot_Eigengene_SP_DE_Gene_Lists_In_MillerLCM <- function(me_DFL){
  gg_L <- lapply(names(me_DFL), function(name){
    gg_DF <- me_DFL[[name]]
    ggplot(gg_DF, aes(x = Zone, y = MESP, fill = Zone)) +
      geom_boxplot() +
      stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "SP") +
      ggplot_set_theme_publication +
      ggtitle(name)
  })
  Plot_Grid(gg_L, ncol = 2, title = paste0(script_name
    , "\n\nMiller LCM expression"
    , "\nEigengene of SP markers"))
  ggsave(paste0(out_graph, "EG_Miller_Boxplot.pdf")
    , width = 7, height = 15)
}

Run_Analysis_SP_cluster_DE_genes <- function(cor_ST18_DF){
  print("Run_Analysis_SP_cluster_DE_genes")
  # DE subplate cluster
  de_sp_vs_all_DF <- DE_SPcluster_Versus_Allcells()
  de_sp_vs_cluster3_DF <- DE_SP_Cluster_Versus_Other_Deep_Layer_Cluster_Cells()
  # Compile gene lists
  gene_group_DF <- Compile_Gene_Lists(
    de_sp_vs_all_DF = de_sp_vs_all_DF
    , de_sp_vs_cluster3_DF = de_sp_vs_cluster3_DF
    , cor_ST18_DF = cor_ST18_DF)
  Output_Compiled_Gene_Lists()
  # Feature plot
  Plot_DE_genes_Feature_Plot(gene_group_DF = gene_group_DF)
  # Plot in Miller LCM
  Plot_SP_DE_Genes_In_Miller_LCM(gene_group_DF = gene_group_DF)
  # Calculate eigengene of DE gene lists
  me_DFL <- Calculate_Eigengene_SP_DE_Gene_Lists_In_MillerLCM(
    gene_group_DF = gene_group_DF
  )
  # Plot eigengene in Miller LCM
  Plot_Eigengene_SP_DE_Gene_Lists_In_MillerLCM(me_DFL = me_DFL)
}
################################################################################

Main_Function <- function(){
  Plot_Subplate_Marker_Expression()
  Plot_ST18_Expression()
  cor_ST18_DF <- Correlation_to_ST18()
  Plot_Expression_Levels_Violin()
  Plot_Intersection_ST18_SPmarks()
  Run_Analysis_SP_cluster_DE_genes(cor_ST18_DF = cor_ST18_DF)
}

Main_Function()
################################################################################


#
#
#
# intersect(
#   tail(de_sp_vs_all_DF$Gene, 50)
#   , tail(de_sp_vs_cluster3_DF$Gene, 50)
# )
#
#
#
#
#
#
#
#
# # Miller LCM from Allen
# in_millerLCM <- list.files("../allen_brain_data/miller_LCM", full.names = TRUE)
# in_millerLCM <- in_millerLCM[grep("lmd", in_millerLCM)]
# # Column annotations
# millerLCM_col_annot_DFL <- lapply(in_millerLCM, function(in_millerLCM){
#   millerLCM_col_annot_DF <- read.csv(paste0(in_millerLCM, "/columns_metadata.csv"))
#   millerLCM_col_annot_DF$Match_Key <- with(millerLCM_col_annot_DF
#     , paste(gsub(".*_", "", in_millerLCM), structure_id))
#   return(millerLCM_col_annot_DF)
# })
# millerLCM_col_annot_DF <- do.call("rbind", millerLCM_col_annot_DFL)
# millerLCM_col_annot_DF$Index <- c(1:nrow(millerLCM_col_annot_DF))
# # Expression matrices
# ex_DFL <- lapply(in_millerLCM, function(in_millerLCM){
#   read.csv(paste0(in_millerLCM, "/expression_matrix.csv"), header = FALSE)
# })
# millerLCM_ex_DF <- do.call("cbind", ex_DFL)
# # Row annotations
# millerLCM_row_annot_DF <- read.csv(
#   "../allen_brain_data/miller_LCM/lmd_matrix_12840/rows_metadata.csv")
#
# webtool_col_order <- read.csv("../allen_brain_data/miller_LCM/webtool_columns.csv")
# webtool_col_order$Match_Key <- with(webtool_col_order
#   , paste(donor_id, structure_id))
#
# # Reorder to webtool sample order
# webtool_col_order <- webtool_col_order[webtool_col_order$Match_Key %in% millerLCM_col_annot_DF$Match_Key, ]
# idx <- match(webtool_col_order$Match_Key, millerLCM_col_annot_DF$Match_Key)
# millerLCM_col_annot_DF <- millerLCM_col_annot_DF[idx, ]
# millerLCM_ex_DF <- millerLCM_ex_DF[ ,idx]
#
# # Get maximum expression probe
# genes <- unique(millerLCM_row_annot_DF$entrez_id)
# keepind <- matrix(nrow = 0, ncol = 0);
# for (ii in 1:length(genes)) {
#   genematchind <- which(millerLCM_row_annot_DF$entrez_id == genes[ii]);
#   if (length(genematchind) > 1) {
#     themeans <- rowMeans(millerLCM_ex_DF[genematchind, ]);
#     maxind <- which(themeans == max(themeans))[1];
#     keepind <- c(keepind, genematchind[maxind]);
#   } else {
#     keepind = c(keepind, genematchind);
#   }
# }
# millerLCM_ex_DF = millerLCM_ex_DF[keepind, ];
# rownames(millerLCM_ex_DF) <- millerLCM_row_annot_DF$entrez_id[
#   match(rownames(millerLCM_ex_DF), rownames(millerLCM_row_annot_DF))]
# millerLCM_row_annot_DF = millerLCM_row_annot_DF[keepind, ];
#
# # Convert entrez id to gene symbols
# # Remove genes with no hgnc symbol
# millerLCM_ex_DF <- millerLCM_ex_DF[
#   rownames(millerLCM_ex_DF) %in% bmDF$entrezgene, ]
# hgnc_symbol <- bmDF$hgnc_symbol[
#   match(rownames(millerLCM_ex_DF), bmDF$entrezgene)]
# millerLCM_ex_DF <- millerLCM_ex_DF[! duplicated(hgnc_symbol), ]
# rownames(millerLCM_ex_DF) <- hgnc_symbol[! duplicated(hgnc_symbol)]
#
#
# Format_Miller_LCM_For_GGplot <- function(
#   millerLCM_ex_DF
#   , genes
#   , millerLCM_col_annot_DF
#   , millerLCM_row_annot_DF){
#   # browser()
#   print("Format_Miller_LCM_For_GGplot")
#   colnames(millerLCM_ex_DF) <-
#     paste(as.character(millerLCM_col_annot_DF$structure_name)
#     , c(1:nrow(millerLCM_col_annot_DF)))
#   millerLCM_ex_DF <- millerLCM_ex_DF[rownames(millerLCM_ex_DF) %in% genes, ]
#   # Change outliers to NA
#   idx_M <- apply(millerLCM_ex_DF, 1, function(expr){
#     mn <- mean(expr)
#     stdev <- sd(expr)
#     idx <- expr > mn + stdev*2 | expr < mn - stdev*2
#     return(idx)
#   })
#   idx_M <- t(idx_M)
#   millerLCM_ex_DF[idx_M] <- NA
#   millerLCM_ex_DF <- as.data.frame(t(scale(t(millerLCM_ex_DF))))
#   millerLCM_ex_DF$Gene <- rownames(millerLCM_ex_DF)
#   millerLCM_ex_DF <- melt(millerLCM_ex_DF)
#   millerLCM_ex_DF$value[millerLCM_ex_DF$value > 1.5] <- 1.5
#   millerLCM_ex_DF$value[millerLCM_ex_DF$value < -1.5] <- -1.5
#   millerLCM_ex_DF$Gene <- factor(millerLCM_ex_DF$Gene
#     , levels = genes
#   )
#   millerLCM_ex_DF$variable <- factor(
#     millerLCM_ex_DF$variable
#     , levels = unique(millerLCM_ex_DF$variable))
#   return(millerLCM_ex_DF)
# }
# ggDF <- Format_Miller_LCM_For_GGplot(
#   millerLCM_ex_DF = millerLCM_ex_DF
#   , millerLCM_col_annot_DF = millerLCM_col_annot_DF
#   , millerLCM_row_annot_DF = millerLCM_row_annot_DF
#   , genes = tail(gene_group_DF$Gene))
#   # , genes = gene_group_DF$Gene[! duplicated(gene_group_DF$Gene)])
#
# ggplot(ggDF, aes(x = variable, y = Gene, fill = value)) +
#   geom_tile() +
#   scale_fill_distiller(name = "Normalized\nexpression", type = "div"
#       , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
#   xlab("Region") +
#   ylab("Gene") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle(paste0(script_name
#     , "\nAllen LCM expression"))
# ggsave(paste0(out_graph, "Miller_Heatmap.png")
#   , width = 8, height = 8, dpi = 150)
