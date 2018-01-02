# Damon Polioudakis
# 2017-08-15
# Run monocle
################################################################################

rm(list = ls())

require(monocle)
require(Seurat)
require(cowplot)
require(ggplot2)
require(viridis)
require(reshape2)
require(biomaRt)
source("Function_Library.R")

# qsub task ID for selecting Seurat cluster ID to run through Monocle
args <- commandArgs(trailingOnly = TRUE)

## Inputs

# Monocle round 2

# How to reorder list of monocle objects for plotting
toOrder <- c("0-1-2-4-12", "0-1-4-12", "0-1-2", "0-1", "3-14", "5-6", "7-9", "0"
  , "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"
  , "15")

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

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# Molyneaux markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

# Telley transcriptional waves
tyDF <- read.csv("../source/Telley_2015_ST3_TranscriptionalWaves.csv", header = TRUE)

# Human TFs
tfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_TF_EnsemblID.txt")
# Human chromatin remodeling factors
crDF <- read.table(
  "../source/AnimalTFDB_Homo_sapiens_chr_remodeling_factor_EnsemblID.txt")
# Human co-factors
cfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_cofactor_EnsemblID.txt")

# biomaRt gene info
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Monocle_Round2_Analysis.R"
outGraph <- "../analysis/graphs/Monocle_Round2_Analysis/Comp1-10/Monocle_Round2_"
outTable <- "../analysis/tables/Monocle_Round2_Analysis/Comp1-10/Monocle_Round2_"
outRdat <- "../analysis/Monocle_Round2_Analysis/Comp1-10/Monocle_Round2_"

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

### Plots of states and clusters
print("### Plots of states and clusters")

## Trajectory plots colored by state or covariate - component 1 and 2
pgL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Visualize trajectory in reduced dimensional space
  ggL <- list(
    plot_cell_trajectory(mo_filtered, 1, 2, color_by = "State", cell_size = 0.01) +
      theme(legend.position = "none") + ggtitle("Colored by state")
    , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "REGION", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "BRAIN", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "LIBRARY", cell_size = 0.01)
  )
  # Change legend size
  ggL <- lapply(ggL, function(gg) {
    gg + guides(colour = guide_legend(override.aes = list(size = 5)))})
  
  # Plot grid
  pg <- plot_grid(plotlist = ggL, ncol = 4, axis = 'b', align = "h")
  # now add the title
  title <- paste("Cluster: ", cluster)
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  return(pg)
})
# Plot grid - 4 columns
pg <- plot_grid(plotlist = pgL, ncol = 1)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by state and covariates"
  , "\nComponent 1 and 2"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_Covariates_Comp12.png")
  , width = 13, height = 95, limitsize = FALSE)

## Trajectory plots colored by state or covariate - component 1 and 3
pgL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Visualize trajectory in reduced dimensional space
  ggL <- list(
    plot_cell_trajectory(mo_filtered, 1, 3, color_by = "State", cell_size = 0.01) +
      theme(legend.position = "none") + ggtitle("Colored by state")
    , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "REGION", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "BRAIN", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "LIBRARY", cell_size = 0.01)
  )
  # Change legend size
  ggL <- lapply(ggL, function(gg) {
    gg + guides(colour = guide_legend(override.aes = list(size = 5)))})
  
  # Plot grid
  pg <- plot_grid(plotlist = ggL, ncol = 4, axis = 'b', align = "h")
  # now add the title
  title <- paste("Cluster: ", cluster)
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  return(pg)
})
# Plot grid
pg <- plot_grid(plotlist = pgL, ncol = 1)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by state and covariates"
  , "\nComponent 1 and 3"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_Covariates_Comp13.png")
  , width = 13, height = 95, limitsize = FALSE)

## Trajectory plots faceted by state - component 1 and 2
pgL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Plot faceted by state
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "State", cell_size = 0.05) +
    facet_wrap(~State, ncol = 4, scales = "free") +
    theme(legend.position = "none") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nMonocle trajectory faceted by state"
      , "\nComponent 1 and 2"
      , "\nCluster ", cluster
      , "\n"))
  ggsave(paste0(outGraph, "Trajectory_Facet_Comp12_Cluster", cluster, ".png")
    , width = 13, height = 2+0.75*length(unique(pData(mo_filtered))))
})
## Trajectory plots faceted by state - component 1 and 3
pgL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Plot faceted by state
  plot_cell_trajectory(mo_filtered, 1, 3, color_by = "State", cell_size = 0.05) +
    facet_wrap(~State, ncol = 4, scales = "free") +
    theme(legend.position = "none") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nMonocle trajectory faceted by state"
      , "\nComponent 1 and 3"
      , "\nCluster ", cluster
      , "\n"))
  ggsave(paste0(outGraph, "Trajectory_Facet_Comp13_Cluster", cluster, ".png")
    , width = 13, height = 2+0.75*length(unique(pData(mo_filtered))))
})

## Trajectory colored by Seurat clusters

pgL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Set Seurat cluster factor levels
  mo_filtered@phenoData$cluster <- factor(mo_filtered@phenoData$cluster
    , levels = sort(as.numeric(as.character(unique(mo_filtered@phenoData$cluster)))))
  
  # Plot Seurat clusters on trajectory
  l <- list(
    plot_cell_trajectory(mo_filtered, 1, 2, color_by = "cluster", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "cluster", cell_size = 0.01)
    , plot_cell_trajectory(mo_filtered, 1, 4, color_by = "cluster", cell_size = 0.01)
  )
  # Change legend size
  l <- lapply(l, function(gg) {
    gg + guides(colour = guide_legend(override.aes = list(size = 7)))
  })
  # extract the legend from one of the plots
  legend <- get_legend(l[[1]])
  # Remove legends from plots
  l <- lapply(l, function(gg) {gg + theme(legend.position = "none")})
  # plot_grid combine tSNE graphs
  pg <- plot_grid(plotlist = l, ncol = 3, align = 'h', axis = 't')
  # add the legend to the row we made earlier. Give it one-third of the width
  # of one plot (via rel_widths).
  pg <- plot_grid(pg, legend, rel_widths = c(3, 2))
  # Cluster title
  title <- ggdraw() + draw_label(paste0("Cluster: ", cluster))
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
  return(pg)
})
# Save
# Plot grid
pg <- plot_grid(plotlist = pgL, ncol = 1)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by Seurat clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_SeuratCluster.png")
  , width = 13, height = 70, limitsize = FALSE)


## Trajectory colored by pseudotime

# Comp 1 and 2
ggL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "Pseudotime"
    , cell_size = 0.01) +
    scale_color_viridis() +
    theme(legend.position = "right") +
    ggtitle(paste0("Cluster: ", cluster))
})
# Save
# Plot grid
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by pseudotime"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_Pseudotime_Comp12.png")
  , width = 13, height = 26, limitsize = FALSE)

# Comp 1 and 3
ggL <- lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  plot_cell_trajectory(mo_filtered, 1, 3, color_by = "Pseudotime"
    , cell_size = 0.01) +
    scale_color_viridis() +
    theme(legend.position = "right") +
    ggtitle(paste0("Cluster: ", cluster))
})
# Save
# Plot grid
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by pseudotime"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_Pseudotime_Comp13.png")
  , width = 13, height = 26, limitsize = FALSE)
################################################################################

### Plots of trajectory colored by expression of genes

## Luis markers

# Subset to marker genes of interest for Luis' excel file
# Cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
kmDFL <- split(kmDF, kmDF$Grouping)

# Luis markers
# Output Directories
outDir <- paste0(outGraph, "Trajectory_LuisMarkers")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_Trajectory_LuisMarkers")
# Plot
lapply(names(moL), function(cluster) {
  mo_filtered <- moL[[cluster]]
  ggL <- lapply(kmDFL, function(kmDF) {
    genes <- kmDF$Gene.Symbol
    print(genes)
    gg <- Plot_Trajectory_Gene_Expression(
      monocleO = mo_filtered
      , genes = genes
      # , exprM = as.matrix(exprs(mo_filtered))
      , exprM = as.matrix(centSO@scale.data)
      , limHigh = 1.5, limLow = -1.5, title = kmDF$Grouping[1])
    gg <- gg + theme(text = element_text(size = 12))
    return(gg)
  })
  pg <- plot_grid(plotlist = ggL, ncol = 3)
  # now add the title
  title <- paste0(graphCodeTitle
    , "\n\nMean expression of groups of published marker genes for cluster: ", cluster
    , "\nNormalized mean centered variance scaled"
    , "\n")
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
  ggsave(paste0(outSubGraph, "_Cluster", cluster, ".png")
    , width = 12, height = 24)
})


## Molyneaux markers
# Output Directories
outDir <- paste0(outGraph, "Trajectory_Molyneaux_NormCentScale")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_Trajectory_Molyneaux_NormCentScale")
# Plot
pgL <- lapply(names(moL)[c(1,2,3,6,7)], function(cluster) {
  mo_filtered <- moL[[cluster]]
  # Molyneaux
  # Seurat normalized centered scaled expression
  ggL <- apply(mmDF, 1, function(v1) {
    genes <- v1[["hgnc_symbol"]]
    print(genes)
    Plot_Trajectory_Gene_Expression(monocleO = mo_filtered
      , genes = genes
      # , exprM = as.matrix(exprs(mo_filtered))
      , exprM = as.matrix(centSO@scale.data)
      , limHigh = 1.5, limLow = -1.5, title = paste(genes, v1[["mgi_symbol"]]))
  })
  plot_grid(plotlist = ggL, ncol = 2)
  pg <- plot_grid(plotlist = ggL, ncol = 3)
  # now add the title
  title <- paste0(graphCodeTitle
    , "\n\nExpression of Molyneaux layer genes for cluster: ", cluster
    , "\nNormalized mean centered variance scaled"
    , "\n")
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.025, 1))
  # Save
  ggsave(paste0(
    outSubGraph, "_Cluster", cluster, ".png")
    , width = 12, height = 90, limitsize = FALSE)
})

# 
# pgL <- lapply(names(moL)[c(1,2,3,6,7)], function(cluster) {
#   mo_filtered <- moL[[cluster]]
#   # Monocle normalized expression
#   # Tmp expression matrix from monocle object
#   m <- as.matrix(exprs(mo_filtered))
#   # Loop through genes
#   ggL <- apply(mmDF, 1, function(v1) {
#     genes <- v1[["hgnc_symbol"]]
#     print(genes)
#     Plot_Trajectory_Gene_Expression(monocleO = mo_filtered
#       , genes = genes
#       # , exprM = as.matrix(exprs(mo_filtered))
#       , exprM = m
#       , limHigh = 3, limLow = -1, title = paste(genes, v1[["mgi_symbol"]]))
#   })
#   pg <- plot_grid(plotlist = ggL, ncol = 3)
#   # now add the title
#   title <- paste0(graphCodeTitle
#     , "\n\nExpression of Molyneaux layer genes for cluster: ", cluster
#     , "\nNormalized mean centered variance scaled"
#     , "\n")
#   title <- ggdraw() + draw_label(title)
#   # rel_heights values control title margins
#   pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.025, 1))
#   ggsave(paste0(
#     outGraph, "Trajectory_Molyneaux_MonocleNormExpr_Cluster", cluster, ".png")
#     , width = 12, height = 90, limitsize = FALSE)
#   # Delete tmp expression matrix
#   rm(m)
# })

# # Jitter expression plot
# blast_genes <- row.names(subset(fData(mo_filtered),
#   gene_short_name %in% mmDF$hgnc_symbol[1:5]))
# plot_genes_jitter(mo_filtered[blast_genes,], grouping = "State", min_expr = 0.1)
# ggsave(paste0(
#   outGraph, "State_Molyneaux_MonocleNormExpr_Cluster", cluster, ".png")
#   , width = 12, height = 12)

# # Pseudotime gene module heatmaps
# diff_test_res <- differentialGeneTest(mo_filtered
#   , fulmoLdelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
#   , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")
################################################################################

### Plotting genes that change as a function of pseudotime

## Histogram of p-values
ggL <- lapply(names(ptDeL), function(cluster){
  
  diff_test_res <- ptDeL[[cluster]]
  
  ggplot(diff_test_res, aes(x = pval)) +
    geom_histogram() +
    ggtitle(paste0("Cluster: ", cluster))
})
pg <- plot_grid(plotlist = ggL, ncol = 4)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nHistogram of p-values testing each gene for DE as function of pseudotime"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "Hist_Pvalue.png"), width = 13, height = 24)

## Histogram of q-values
ggL <- lapply(names(ptDeL), function(cluster){
  
  diff_test_res <- ptDeL[[cluster]]
  
  ggplot(diff_test_res, aes(x = qval)) +
    geom_histogram() +
    ggtitle(paste0("Cluster: ", cluster))
})
pg <- plot_grid(plotlist = ggL, ncol = 4)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nHistogram of q-values testing each gene for DE as function of pseudotime"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
ggsave(paste0(outGraph, "Hist_Qvalue.png"), width = 13, height = 24)

## Heatmap of top 100 genes
# Output Directories
outDir <- paste0(outGraph, "Pseudotime_Heatmap")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_Pseudotime_Heatmap")
# Plot
lapply(names(ptDeL), function(cluster){
  
  diff_test_res <- ptDeL[[cluster]]
  mo_filtered <- moL[[cluster]]
  
  # Order by qval
  diff_test_res <- diff_test_res[order(diff_test_res$qval), ]
  
  # Subset to top 100 by qval
  sig_gene_names <- row.names(diff_test_res)[
    diff_test_res$use_for_ordering == TRUE][1:100]
  
  # Plot
  png(paste0(outSubGraph, "_Cluster", cluster, ".png")
    , width = 12, height = 16, units = "in", res = 300)
  tryCatch(plot_pseudotime_heatmap(mo_filtered[sig_gene_names,]
    , num_clusters = 3
    , cores = 1
    , show_rownames = T)
    , error = function(cond) {
      message(paste0("Error for cluster: ", cluster))
      message("Error message:")
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    }
  )
  dev.off()
})


## Heatmap of TFs, co-factors, chromatin remodelers

# Data frame of TFs, co-factors, chromatin remodelers
tccDF <- rbind(tfDF, crDF, cfDF)
# 1970
nrow(tccDF)
# Add gene symbol
tccDF$GENE <- bmDF$hgnc_symbol[match(tccDF$V1, bmDF$ensembl_gene_id)]
# Check match using ensembl
tccDF$ENSEMBL <- bmDF$ensembl_gene_id[match(tccDF$V1, bmDF$ensembl_gene_id)]

# Heatmap
# Output Directories
outDir <- paste0(outGraph, "Pseudotime_Heatmap_TfCfCr")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_Pseudotime_Heatmap_TfCfCr")
# Plot
lapply(names(ptDeL), function(cluster){
  
  diff_test_res <- ptDeL[[cluster]]
  mo_filtered <- moL[[cluster]]
  
  # Subset to TFs, co-factors, chromatin remodelers and by qval
  diff_test_res <- diff_test_res[diff_test_res$qval < 0.05 & diff_test_res$gene_short_name %in% tccDF$GENE, ]
  
  # Order by qval
  diff_test_res <- diff_test_res[order(diff_test_res$qval), ]
  
  genes <- row.names(diff_test_res)
  
  png(paste0(outSubGraph, "_Cluster", cluster, ".png")
    , width = 12, height = 16, units = "in", res = 300)
  tryCatch(plot_pseudotime_heatmap(mo_filtered[genes,]
    , num_clusters = 3
    , cores = 1
    , show_rownames = T)
    , error = function(cond) {
      message(paste0("Error for cluster: ", cluster))
      message("Error message:")
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    }
  )
  dev.off()
})

# Heatmap of state DE genes
# Output Directories
outDir <- paste0(outGraph, "StateDE_Heatmap")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_StateDE_Heatmap")
# Plot
lapply(names(stDeLDF), function(cluster) {
  stateDeDF <- stDeLDF[[cluster]]
  mo <- moL[[cluster]]
  nCellsDF <- data.frame(table(pData(mo)$State))
  print(nCellsDF)
  apply(nCellsDF, 1, function(y){
    print(y)
    state <- as.numeric(as.character(y[["Var1"]]))
    nCells <- as.numeric(as.character(y[["Freq"]]))
    if (nCells > 50) {
      # stateDeDF <- stateDeDF[order(stateDeDF$qval), ]
      stateDeDF <- stateDeDF[stateDeDF$LOG_FC > 0.7 & stateDeDF$STATE == state, ]
      if (nrow(stateDeDF) > 0) {
        print(paste0("Graphing state ", state))
        # stateDeDF <- stateDeDF[stateDeDF$qval < 0.05]
        geneGroupDF <- data.frame(GENE = stateDeDF$GENE)
        geneGroupDF$GROUP <- ""
        
        Heatmap_By_State(geneGroupDF, exprM = noCentExM, monocleO = mo
          , clusters = 1:max(as.numeric(as.character(pData(mo)$State)))
          , lowerLimit = -1.5, upperLimit = 1.5, geneOrder = geneGroupDF$GENE
          , centScale = TRUE)
        ggsave(paste0(
          outSubGraph, "_Cluster", cluster, "_State", state, ".png")
          , height = 3+0.2*nrow(geneGroupDF))
        # # plot_grid combine
        # pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
        # # now add the title
        # title <- ggdraw() + draw_label(title)
        # # rel_heights values control title margins
        # pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
      }
    }
  })
})

# # Jitter of DE genes
# lapply(names(stDeLDF), function(cluster) {
#   stateDeDF <- stDeLDF[[cluster]]
#   mo <- moL[[cluster]]
#   plot_genes_jitter(mo[row.names(stateDeDF)[1], ], grouping = "State", color_by = "State"
#     , nrow = 1, ncol = NULL, plot_trend = TRUE)
#   ggsave(paste0(outGraph, "StateDE_Jitter_Cluster", cluster, ".png")
#     , width = 13, height = 6)
# })

beamDF <- beamLLDF[["4"]][[1]]
mo_filtered <- moL[["4"]]
png(paste0(outGraph, "BranchDE_Heatmap_Cluster_Branch.png")
  , width = 12, height = 16, units = "in", res = 300)
plot_genes_branched_heatmap(mo_filtered[row.names(subset(beamDF, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4,
  cores = 1,
  use_gene_short_name = T,
  show_rownames = T)
dev.off()
################################################################################
  
  
  
  
  