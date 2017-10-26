# Damon Polioudakis
# 2017-06-15
# Interneuron marker expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3 +
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
require(cowplot)
source("Function_Library.R")

## Inputs

# Log normalized, regressed nUMI and percent mito
# seuratO
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")

# biomaRt gene info
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Interneuron_Marker_Expression.R"
outGraph <- "../analysis/graphs/Interneuron_Marker_Expression_"
outTable <- "../analysis/tables/Interneuron_Marker_Expression_"

## Output Directories
outGraphDir <- dirname(outGraph)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Plots

# Genes to check
genes <- c("DLX1", "DLX2", "DLX5", "DLX6", "SOX6", "TTF1", "LHX6", "OLIG1"
  , "OLIG2", "MKI67", "BUB1")

## Violin plots of expression by cluster
ggL <- lapply(genes, function(gene) {
  # Expression per cluster
  ggDF <- noCentExM
  ggDF <- data.frame(EXPRESSION = ggDF[row.names(ggDF) == gene, ])
  ggDF$CLUSTER <- centSO@ident
  ggplot(ggDF, aes(x = CLUSTER, y = EXPRESSION)) +
    geom_violin(aes(fill = CLUSTER)) +
    geom_jitter(size = 0.05, height = 0, alpha = 0.3) +
    theme(legend.position = "none") +
    ylab("Normalized expression") +
    xlab("Clusters") +
    ggtitle(gene)
})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'l')
# now add the title
title <- paste0(graphCodeTitle
  , "\n"
  , "\nInterneuron marker expression by cluster"
  , "\nNormalized expression"
  , "\n")
# rel_heights values control title margins
title <- ggdraw() + draw_label(title)
plot_grid(title, pg, ncol = 1
  , rel_heights = c(length(ggL)*0.1, length(ggL)))
ggsave(paste0(outGraph, "violinPlots.png"), width = 14, height = 4+length(ggL))

## Feature plots

# Collect tSNE values for ggplot
ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)

# Feature plot
# Normalized, mean centered scaled
# Loop through and plot each group of genes
ggL <- lapply(genes, function(gene) {
  print(gene)
  ggDF <- Mean_Expression(ggDF, gene, centSO@scale.data)
  ggDF <- Set_Limits(ggDF, limLow = -1.5, limHigh = 1.5)
  ggFp <- Feature_Plot_CentScale(ggDF, limLow = -1.5, limHigh = 1.5
    , title = paste0("\n", gene)
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
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'l')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE plot, each point is a cell"
  , "\nColor indicates normalized expression, mean centered, variance scaled"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1
  , rel_heights = c(length(ggL)*0.05, length(ggL)))
ggsave(paste0(
  outGraph, "FeaturePlot_NormCentScale.png")
  , width = 14, height = 5+length(ggL), limitsize = FALSE)

# Feature plot
# Normalized, not mean centered scaled
# Loop through and plot each group of genes
ggL <- lapply(genes, function(gene) {
  print(gene)
  ggDF <- Mean_Expression(ggDF, gene, noCentExM)
  ggDF <- Set_Limits(ggDF, limLow = -1, limHigh = 3)
  ggFp <- Feature_Plot(ggDF, limLow = -1, limHigh = 3
    , title = paste0("\n", gene)
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
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'l')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(3, 0.3))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\ntSNE plot, each point is a cell"
  , "\nColor indicates normalized expression"
  , "\n"))
# rel_heights values control title margins
plot_grid(title, pg, ncol = 1
  , rel_heights = c(length(ggL)*0.05, length(ggL)))
ggsave(paste0(
  outGraph, "FeaturePlot_Norm.png")
  , width = 14, height = 5+length(ggL), limitsize = FALSE)

## Heatmaps

# Normalized, mean centering scaling
geneGroupDF <- data.frame(GENE = genes, GROUP = "")
ggL <- Heatmaps_By_Cluster_Combined(geneGroupDF, exprM = centSO@scale.data
  , seuratO = centSO, lowerLimit = -1.5, upperLimit = 1.5
  , clusters1 = c(0:1), clusters2 = c(2:10), clusters3 = c(11:17)
)
# plot_grid combine
pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of genes of interest"
  , "\nx-axis: Genes"
  , "\ny-axis: Cells ordered by cluster"
  , "\nNormalized expression, mean centered, variance scaled"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.3, 1))
ggsave(paste0(outGraph, "ExprHeatmap_NormalizedCenteredScaled.png")
  , width = 12, height = 4+nrow(geneGroupDF)/3, limitsize = FALSE)
################################################################################
# ## Assign annotated cluster names to clusters
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
# new.cluster.ids <- c(
#   "Excitatory Upper Layer Neuron 1"
#   , "Excitatory Neuron"
#   , "Excitatory Upper Layer Neuron 2"
#   , "Excitatory Deep Layer Neuron"
#   , "Intermediate Progenitors"
#   , "Interneuron"
#   , "Mitotic Progenitors"
#   , "oRG"
#   , "Oligodendrocyte Precursor"
#   , "Endothelial")
# seuratO@ident <- plyr::mapvalues(seuratO@ident, from = current.cluster.ids
#   , to = new.cluster.ids)
# # Stash numerical cluster identities if want to use later
# seuratO <- StashIdent(seuratO, save.name = "Cluster_Numbers")


inExDF <- seuratO@scale.data[ ,seuratO@ident == "Interneuron"]

genes <- c("DLX1", "DLX2", "SOX6", "TTF1", "LHX6", "OLIG1", "OLIG2", "MKI67", "BUB1")

bmDF[bmDF$ensembl_gene_id == "ENSG00000136352", ]

row.names(inExDF)[row.names(inExDF) %in% c("NKX2A", "TITF1", "TTF-1", "TTF1"
  , "NMTC1", "NK-2", "TEBP", "BCH", "BHC", "NKX2-1")]

ggDF <- data.frame(t(inExDF[row.names(inExDF) %in% genes, ])
  , CELL_ID = colnames(inExDF))
ggDF$REGION <- "GZ"
ggDF$REGION[grep("CP", ggDF$CELL_ID)] <- "CP"
ggDF <- melt(ggDF)
ggDF$variable <- as.character(ggDF$variable)
ggDF$variable[ggDF$variable == "TTF1"] <- "NKX2.1"
unique(ggDF$variable)

# Boxplot
ggplot(ggDF, aes(x = variable, y = value, fill = REGION)) +
  geom_boxplot() +
  xlab("Genes") +
  ylab("Log2 normalized expression") +
  ggtitle("Expression in Interneuron Cluster")
ggsave(paste0(outGraph, "Boxplot.png"))

# Boxplot for only cells expressing gene of interest
gg2DF <- ggDF[ggDF$value > 0.5, ]
ggplot(gg2DF, aes(x = variable, y = value, fill = REGION)) +
  geom_boxplot() +
  xlab("Genes") +
  ylab("Log2 normalized expression") +
  ggtitle(paste0("Expression in Interneuron Cluster"
    , "\nSubset to cells in which gene is expressed"))
ggsave(paste0(outGraph, "Expressing_Boxplot.png"))

gg3DF <- data.frame(LHX6 = ggDF$value[ggDF$variable == "LHX6"]
  , NKX2.1 = ggDF$value[ggDF$variable == "NKX2.1"])
ggplot(gg3DF, aes(x = LHX6, y = NKX2.1)) +
  geom_point(shape = 1) +
  ggtitle(paste0("Log2 normalized expression for cells in interneuron cluster"))
ggsave(paste0(outGraph, "LHX6_vs_NKX21.png"))

# Heatmap
# Set limits
ggDF$value[ggDF$value < 0] <- 0
ggDF$value[ggDF$value > 3] <- 3
ggplot(ggDF, aes(x = CELL_ID, y = variable, fill = value)) +
  geom_tile() +
  facet_grid(~REGION, space = "free", scales = "free") +
  scale_fill_distiller(name = "Log2 normalized\nexpression", type = "div"
    , palette = 5, direction = -1, limits = c(0, 3)) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90)) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(text = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab("Genes") +
  xlab("Cells") +
  ggtitle("Expression in Interneuron Cluster")
ggsave(paste0(outGraph, "Heatmap.png"))

# Table of percent expressing
with(gg2DF, by(gg2DF, variable, function(x) t.test(value ~ REGION, data = x)))

df <- tapply(ggDF$value, list(ggDF$variable, ggDF$REGION)
  , function(v1) sum(v1 > 0.5)/length(v1))
df <- round(df, 3) * 100
write.csv(df, paste0(outTable, "Expressing_Percent.csv"), quote = FALSE)