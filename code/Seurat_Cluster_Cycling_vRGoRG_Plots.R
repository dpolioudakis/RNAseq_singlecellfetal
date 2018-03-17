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
source("Function_Library.R")
source("Seurat_Cluster_Cycling_vRGoRG_Functions.R")

options(stringsAsFactors = FALSE)

## Inputs

# Keep CC genes from variable gene list used for clustering
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TESTcluster0278910_seuratO.Robj")
# centSO <- ssCentSO; rm(ssCentSO)
# noCentExM <- ssNoCentExM; rm(ssNoCentExM)

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

# DE cell types
deDF <- read.csv(paste0(outTable, "DE_CellTypes.csv"))

# DE mixed marker cells
transition_state_DE_DF <- read.csv(
  paste0(outTable, "DE_CellTransitionStates.csv"))

# Eigengenes of DE genes
# ME_CellType_L
load(file = paste0(outData, "ME_CellType.rdata"))
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

  print("Plot_ME_CellType_Genes")

  me_markerFlag_DFL <- split(me_markerFlag_DF
    , list(me_markerFlag_DF$Fold_Change_Cutoff
    , me_markerFlag_DF$variable)
  )

  ggL <- lapply(names(me_markerFlag_DFL), function(name){

    print(name)

    me_markerFlag_DF <- me_markerFlag_DFL[[name]]

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
      scale_fill_manual(values = color_1) +
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
      #
      # gg <- plot_grid(gg, ncol = 1, rel_heights = c(1,0.5))
      return(gg)
  })
  ggL[4:6] <- lapply(ggL[4:6], function(gg){
    gg + scale_fill_manual(values = color_2)
  })
  return(ggL)
}

# Plotting

# RG IP cluster 8
gg1L <- Plot_ME_CellType_Genes(
  me_markerFlag_DF = ME_CellType_L[["RG_to_IP_8"]]
  , color_1 = "#8dd3c7", color_2 = "#bebada"
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
  , color_1 = "#8dd3c7", color_2 = "#bebada"
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
  , color_1 = "#8dd3c7", color_2 = "#fb8072"
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
  , color_1 = "#8dd3c7", color_2 = "#fb8072"
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
  , color_1 = "#bebada", color_2 = "#fb8072"
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
  , color_1 = "#bebada", color_2 = "#fb8072"
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
  # RG IP
  gg1L[2]
  , gg2L[2]
  , gg1L[5]
  , gg2L[5]
  # RG neuron
  , gg3L[2]
  , gg4L[2]
  , gg3L[5]
  , gg4L[5]
  # IP neuron
  , gg5L[2]
  , gg6L[2]
  , gg5L[5]
  , gg6L[5]
))
ggL <- lapply(ggL, function(gg){
  gg <- gg + ggtitle("") + xlab("")
  return(gg)
})
Plot_Grid(ggL, ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nME of DE genes")
)
ggsave(paste0(outGraph, "DE_ME_boxplot_paper.pdf")
  , height = 16, width = 6)

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
  percent_signature_DF <- data.frame(
    Percent = c(length(intersect(genes1, genes2))/length(genes1) * 100
      , length(intersect(genes3, genes4))/length(genes3) * 100)
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
  S_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_Sphase"
    , fold_change = 0.01
    , fdr = 1.1
    , transition_state = "RG vs RGIP"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_IP"
    , comparison_2 = "RG_vs_RGIP_G2Mphase"
    , fold_change = 0.01
    , fdr = 1.1
    , transition_state = "RG vs RGIP"
    , phase = "G2M phase"
  )

  , S_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_Sphase"
    , fold_change = 0.01
    , fdr = 1.1
    , transition_state = "RG vs RGNeuron"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
    comparison_1 = "RG_vs_Neuron"
    , comparison_2 = "RG_vs_RGNeuron_G2Mphase"
    , fold_change = 0.01
    , fdr = 1.1
    , transition_state = "RG vs RGNeuron"
    , phase = "G2M phase"
  )

  , S_phase = Percent_Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_Sphase"
    , fold_change = 0.01
    , fdr = 1.1
    , transition_state = "IP vs IPNeuron"
    , phase = "S phase"
  )
  , G2M_phase = Percent_Signature_DE_Genes(
      comparison_1 = "IP_vs_Neuron"
    , comparison_2 = "IP_vs_IPNeuron_G2Mphase"
    , fold_change = 0.01
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
ggplot(ggDF, aes(x = Phase, y = Percent, fill = Signature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(~Transition_State, scales = "free") +
  ylim(c(0,100)) +
  scale_fill_manual(values = c("#8dd3c7", "#bebada", "#fb8072")) +
  ggplot_set_theme_publication +
  ggtitle(paste0(graphCodeTitle
    , "\n\nPercent of DE signature genes DE in expected direction in mixed marker cells"))
ggsave(paste0(outGraph, "DE_PercentMixed_FDR1_FC0.pdf")
  , width = 8, height = 4)


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

  subset_transition_state_DE_DF$DE_Group_2 <- as.factor(DE_group_label)
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
      , "\nDE Signature: ", ggDF$DE_Group_2
    ))
})
Plot_Grid(ggL, ncol = 4, title = "DE Mixed")
ggsave(paste0(outGraph, "DE_Mixed_barplot.png"), width = 11, height = 9)
################################################################################

### Percent of cells passing expression filters per cluster

Mixed_Marker_By_Cluster_Percent_Barplot <- function(
  exM, seuratO, highThreshold, lowThreshold, title, cluster_order = NULL) {

  df <- Average_MarkersExp_Per_Cell(exM = exM, seuratO = seuratO)
  df <- Positive_Negative_Expression_Flag(
    exDF = df, highThreshold = 0.5, lowThreshold = 0.5)
  df$TYPE[df$TYPE %in% c("Neuron+", "IP+", "RG+", "RG+ vRG- oRG-")] <- "NA"
  df$TYPE <- droplevels(df$TYPE)
  df <- aggregate(df$TYPE, list(df$CLUSTER), Percent_Of_Table)

  df2 <- as.data.frame(df[ ,2])
  df2$CLUSTER <- df[ ,1]
  df <- melt(df2)

  # Set cluster order
  if (! is.null(cluster_order)){
    df$CLUSTER <- factor(df$CLUSTER, levels = cluster_order)
  }
  df <- df[! is.na(df$CLUSTER), ]

  gg <- ggplot(df, aes(x = CLUSTER, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
    # scale_fill_brewer(type = "qual", palette = "Set2", direction = 1) +
    xlab("Cluster") +
    ylab("Percent of cells") +
    ggtitle(title)

  return(gg)
}

# Plot
gg1 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.5
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.5 normalized expression"
)
gg2 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.5 normalized expression\n- = < 0.25 normalized expression"
)
gg3 <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.75, lowThreshold = 0.25
  , title = "Keep CC\n+ = > 0.75 normalized expression\n- = < 0.25 normalized expression"
)
# Plot grid
pg <- plot_grid(gg1, gg2, gg3, ncol = 2)
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nPercent of cells passing combinations of marker expression filters"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "PercentMixedCluster_Barplot.pdf")
  , width = 12, height = 14)

# Paper - plot with reordered clusters
gg <- Mixed_Marker_By_Cluster_Percent_Barplot(
  exM = noCentExM, seuratO = centSO
  , highThreshold = 0.5, lowThreshold = 0.25
  , cluster_order = c(9,7,8,10,2,0,1,12,4,3,14,5,6,11,13,15,16)
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
  ggplot_set_theme_publication
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
  ggplot_set_theme_publication
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
  ggplot_set_theme_publication
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
