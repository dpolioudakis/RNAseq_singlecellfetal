# Damon Polioudakis
# 2017-12-09
# Analyize Monocle output using Telley gene lists
################################################################################

rm(list = ls())

require(monocle)
require(Seurat)
require(cowplot)
require(ggplot2)
require(viridis)
require(reshape2)
require(biomaRt)
require(ggdendro)
source("Function_Library.R")

# qsub task ID for selecting Seurat cluster ID to run through Monocle
args <- commandArgs(trailingOnly = TRUE)

## Inputs

# Monocle round 2

# How to reorder list of monocle objects for plotting
toOrder <- c("0-1-4-12", "0-1-2", "0-1", "3-14", "5-6", "7-9", "0", "1", "2", "3"
  , "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")

# Monocle object
inMo <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inMo <- inMo[grep("Monocle_Round2_monocleO", inMo)]
moL <- lapply(inMo, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(mo_filtered)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_monocleO_cluster", "", inMo)
names <- gsub(".Robj", "", names)
names(moL) <- names
# Reorder list of monocle objects for plotting
moL <- moL[toOrder]

# Pseudotime DE
inPtDE <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inPtDE <- inPtDE[grep("Monocle_Round2_DEpseudotime_cluster", inPtDE)]
ptDeL <- lapply(inPtDE, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(ptDE)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEpseudotime_cluster", "", inPtDE)
names <- gsub(".Robj", "", names)
names(ptDeL) <- names
# Reorder list of DE for plotting
ptDeL <- ptDeL[toOrder]

# State DE
inStDE <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inStDE <- inStDE[grep("DEstate_cluster", inStDE)]
# inStDE <- inStDE[8:15]
stDeLDF <- lapply(inStDE, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  # load("../analysis/Monocle_Round2/Monocle_Round2_Comp1-20_DEstate_cluster5-6.Robj")
  return(stateDE)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEstate_cluster", "", inStDE)
names <- gsub(".Robj", "", names)
names(stDeLDF) <- names
# Reorder list of DE for plotting
stDeLDF <- stDeLDF[toOrder]

# Branch DE
inBeam <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inBeam <- inBeam[grep("DEbranch_cluster", inBeam)]
beamLLDF <- lapply(inBeam, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(beamLDF)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEbranch_cluster", "", inBeam)
names <- gsub(".Robj", "", names)
names(beamLLDF) <- names
# Reorder list of DE for plotting
beamLLDF <- beamLLDF[toOrder]

# Seurat clustering
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# BrainSpan developmental transcriptome
bsDF <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/expression_matrix.csv", header = FALSE)
rnames <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/rows_metadata.csv")
row.names(bsDF) <- rnames$gene_symbol
bsMtDF <- read.csv("../source/BrainSpan_DevTranscriptome/genes_matrix_csv/columns_metadata.csv")

## Variables
graphCodeTitle <- "Monocle_Round2_Analysis_PseudoMaturation.R"
outGraph <- "../analysis/graphs/Monocle_Round2_Analysis_PseudoMaturation/Comp1-10/Monocle_Round2_Analysis_PseudoMaturation_"
outTable <- "../analysis/tables/Monocle_Round2_Analysis_PseudoMaturation/Comp1-10/Monocle_Round2_Analysis_PseudoMaturation_"
outRdat <- "../analysis/Monocle_Round2_Analysis_PseudoMaturation/Comp1-10/Monocle_Round2_Analysis_PseudoMaturation_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outRdat), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

Seurat_Heatmap_By_Cluster_Hclust_Genes <- function(genes, exM) {
  
  # Subset expression matrix
  exM <- exM[row.names(exM) %in% genes, ]
  
  # Center scale gene expression by cell
  exM <- t(scale(t(exM)))
  
  # Obtain the dendrogram
  dend <- as.dendrogram(hclust(d = dist(exM), method = "ward.D2"))
  dend_data <- dendro_data(dend)
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  # Use the dendrogram label data to position the gene labels
  gene_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, gene = as.character(label), height = 1))
  # Limits for the vertical axes
  gene_axis_limits <- with(
    gene_pos_table,
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1
  # Facet dendrogram so it lines up with faceted heatmaps
  segment_data$Facet <- ""
  
  # Dendrogram plot
  plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_x_reverse(expand = c(0, 0.5)) + 
    scale_y_continuous(breaks = gene_pos_table$y_center, 
      labels = gene_pos_table$gene, 
      limits = gene_axis_limits,
      expand = c(0, 0)) + 
    facet_wrap(~Facet) +
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) +
    theme(strip.background = element_blank())
  
  ggDF <- melt(exM)
  ggDF$Var1 <- factor(ggDF$Var1, levels = gene_pos_table$gene)
  
  idx <- match(ggDF$Var2, row.names(mo_filtered@phenoData@data))
  ggDF$Pseudotime <- mo_filtered@phenoData@data$Pseudotime[idx]
  ggDF <- ggDF[order(ggDF$Pseudotime), ]
  ggDF$Var2 <- factor(ggDF$Var2, levels = unique(ggDF$Var2))

  ggDF$value[ggDF$value > limHigh] <- limHigh
  ggDF$value[ggDF$value < limLow] <- limLow
  # ggplot
  gg <- ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    # facet_grid(GROUP~SEURAT_CLUSTERS, space = "free", scales = "free") +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10))
  
  
   
  
  # Heatmap plot
  geneGroupDF <- data.frame(GENE = row.names(exM), GROUP = "")
  ggL <- lapply(c(0:17), function(cluster){
    tryCatch(
      Heatmap_By_Cluster(
        geneGroupDF = geneGroupDF
        , exprM = as.matrix(centSO@scale.data)
        , seuratO = centSO
        , clusters = cluster
        , lowerLimit = -1.5
        , upperLimit = 1.5
        , geneOrder = gene_pos_table$gene
      )
      , error = function(e) NULL)
  })
  # Remove nulls from ggplot list
  ggL <- ggL[! sapply(ggL, is.null)]
  # Extract legend
  legend <- get_legend(ggL[[1]])
  # Format - remove axis labels
  # ggL[[1]] <- ggL[[1]] + theme(
  #   axis.title.x = element_blank()
  #   , legend.position = "none"
  #   , strip.text.y = element_blank())
  ggL[1:length(ggL)] <- lapply(ggL[1:length(ggL)], function(gg) {
    gg + theme(
      strip.text.y = element_blank()
      , legend.position = "none"
      , axis.title.y = element_blank()
      , axis.text.y = element_blank()
      , axis.ticks.y = element_blank()
      , axis.title.x = element_blank()
      # margin: top, right, bottom, and left
      , plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")
    )
  })
  
  # Combine individual heatmaps and dendrogram
  rel_widths <- as.vector(log((table(centSO@ident) + 1), 10)) + 1
  rel_widths <- c(20, rel_widths, 1)
  # Combine
  pg <- plot_grid(plotlist = list(plt_dendr, gg), ncol = 2
    , rel_widths = c(1,2), align = 'h', axis = 'b')
  
  return(pg)
}
################################################################################

exM <- noCentExM
genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.01]



pg <- plot_grid(plotlist = list(plt_dendr, gg), ncol = 2
  , rel_widths = c(1,2), align = 'h', axis = 'b')
ggsave(paste0(outGraph, "Cluster", paste0(cluster, collapse = "-"), "heatmap.png")
  , width = 7, height = 10)

png(paste0(outGraph, "Cluster", paste0(cluster, collapse = "-"), "heatmapMonocle.png")
  , width = 7, height = 10, units = "in", res = 300)
plot_pseudotime_heatmap(mo_filtered[as.character(genes),]
  , cores = 1
  , show_rownames = FALSE)
  # , ~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs)
dev.off()
################################################################################

### Brainspan developmental transcriptome

unique(bsMtDF$structure_name)
# Structures with cortex in name:
# [1] "occipital neocortex"
# [2] "primary motor-sensory cortex (samples)"
# [3] "posterior (caudal) superior temporal cortex (area 22c)"
# [4] "anterior (rostral) cingulate (medial prefrontal) cortex"
# [5] "dorsolateral prefrontal cortex"
# [6] "orbital frontal cortex"
# [7] "inferolateral temporal cortex (area TEv, area 20)"
# [8] "ventrolateral prefrontal cortex"
# [9] "parietal neocortex"
# [10] "temporal neocortex"
# [11] "primary auditory cortex (core)"
# [12] "primary visual cortex (striate cortex, area V1/17)"
# [13] "primary motor cortex (area M1, area 4)"
# [14] "posteroventral (inferior) parietal cortex"
# [15] "primary somatosensory cortex (area S1, areas 3,1,2)"
# [16] "cerebellar cortex"

# Subset to cortex, remove visual and cerebellar
df <- bsMtDF[grep("cortex", bsMtDF$structure_name), ]
df[df$structure_name != "cerebellar cortex" &
    df$structure_name != "primary visual cortex (striate cortex, area V1/17)", ]
ssBsDF <- bsDF[ ,df$column_num]
# Subset to genes of interest
idx <- match(c(
  # RG
  "VIM", "HES1", "HOPX", "ITGB5", "CARHSP1", "ZFHX4", "LITAF", "MAFF"
  # Excitatory - not deep layer
  , "TUBB3", "STMN2", "SATB2", "CUX1", "CSRP2"
  # Excitatory deep layer
  , "BCL11B", "SOX5", "TBR1", "LCORL", "ST18", "KAT6B"
  # Excitatory - upper or deep layer
  , "NEUROD6", "MAP2", "YWHAB"
  # Interneuron
  , "DLX1", "DLX5", "CXCR4", "CALB2", "CITED2"), rnames$gene_symbol)
df1 <- rnames[idx, ]
ssBsDF <- ssBsDF[df1$row_num, ]
row.names(ssBsDF) <- df1$gene_symbol
# Genes as columns
ssBsDF <- as.data.frame(t(ssBsDF))
# Remove outliers by stdev
ssBsDF <- as.data.frame(apply(ssBsDF, 2, function(x) {
  Remove_Outliers_By_SD(x, nStdev = 2.5)
}))
# Add age
ssBsDF$AGE <- factor(df$age, levels = unique(bsMtDF$age))

# mean expression at each time point
ssBsDF <- melt(ssBsDF)
mnDF <- aggregate(value~variable+AGE, data = ssBsDF, mean, na.rm = TRUE)

# Split by gene and plot
ldf <- split(ssBsDF, ssBsDF$variable)
# Duplicate some genes to plot next to genes of interest
idx <- match(c(
  # RG
  "VIM", "HES1", "HOPX", "ITGB5", "CARHSP1", "ZFHX4", "LITAF", "MAFF"
  # Excitatory - not deep layer
  , "TUBB3", "STMN2", "SATB2", "CUX1", "CSRP2"
  # Excitatory deep layer
  , "BCL11B", "SOX5", "TBR1", "LCORL", "ST18", "KAT6B"
  # Excitatory - upper or deep layer
  , "TUBB3", "STMN2", "NEUROD6", "MAP2", "YWHAB"
  # Interneuron
  , "DLX1", "DLX5", "CXCR4", "CALB2", "CITED2"), names(ldf))
ldf <- ldf[idx]
names(ldf)
# Loop through and plot
ggL <- lapply(ldf, function(df) {
  # ggplot
  ggplot(df, aes(x = AGE, y = value, group = 1)) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    geom_jitter(size = 0.1, width = 0.2) +
    stat_summary(fun.y = "mean", color = "red", geom = "line") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # ggtitle(df$variable[1])
    xlab("Age") +
    ylab("Normalized expression")
})
# plot grid
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n"
  , "\nExpression across Brainspan developmental transcriptome"
  , "\nOutliers >2.5 SD from mean expression removed"
  , "\n"
)
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.075, 1))
# Save
ggsave(paste0(outGraph, "BrainSpan.png"), width = 14, height = 30)
################################################################################