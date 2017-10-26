# Damon Polioudakis
# 2017-08-15
# Run monocle
################################################################################

rm(list = ls())

require(monocle)
require(Seurat)
require(cowplot)
require(ggplot2)
source("Function_Library.R")

# load("../analysis/Monocle/Monocle_monocleO.Robj")
load("../analysis/Monocle/Monocle_PC1-40_RG_IPC_Neuron_monocleO.Robj")

## Inputs

# Seurat object
load("../analysis/Seurat_Cluster_DS2-11/Seurat_Cluster_DS2-11_FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_PC1to40_seuratO.Robj")

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# Molyneaux markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

## Variables
graphCodeTitle <- "Monocle.R"
outGraph <- "../analysis/graphs/Monocle/Monocle_PC1-40_RG_IPC_Neuron_"
outTable <- "../analysis/tables/Monocle/Monocle_PC1-40_RG_IPC_Neuron_"
outRdat <- "../analysis/Monocle/Monocle_PC1-40_RG_IPC_Neuron_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdat)
dir.create(outRdatDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Run Monocle

# Uses size factor normalization:
# The size factor is the median ratio of the sample over a "pseudosample": for
# each gene, the geometric mean of all samples.

# Subset Seurat object to specifc columns
ids <- c(0, 1, 2, 3, 4, 7, 8, 9, 10, 13, 14)
cellIDs <- names(centSO@ident)[centSO@ident %in% ids]
centSO <- FilterCells(object = centSO, subset.names = NULL, cells.use = cellIDs)

# Add clustering to metadata
v1 <- centSO@meta.data$res.0.6
names(v1) <- row.names(centSO@meta.data)
centSO <- AddMetaData(centSO, v1, "cluster")
# row.names(metDF) <- metDF$CELL
# metDF <- metDF[metDF$CELL %in% row.names(centSO@meta.data), ]
# centSO <- AddMetaData(centSO, metDF)

# Subset cells to Seurat filtered cells
exDF <- centSO@raw.data
exDF <- exDF[colnames(exDF) %in% colnames(centSO@scale.data)]
metDF <- centSO@meta.data
metDF$BRAIN <- as.factor(metDF$BRAIN)
# Order metadata to match expression matrix colnames
idx <- match(colnames(exDF), row.names(metDF))
metDF <- metDF[idx, ]

# # Subsetting for testing
# idx <- sample(1:ncol(exDF), 2000)
# exDF <- exDF[ ,idx]
# metDF <- metDF[idx, ]

feature_data = data.frame(gene_short_name = rownames(exDF))
rownames(feature_data) = feature_data$gene_short_name

pd <- new("AnnotatedDataFrame", data = metDF)
fd <- new("AnnotatedDataFrame", data = feature_data) 

mo <- newCellDataSet(cellData = as(as.matrix(exDF), "sparseMatrix"), 
  phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, 
  expressionFamily = negbinomial.size())

# Filter genes and cells
# retains all 
mo <- detectGenes(mo, min_expr = 0.1)
print(head(fData(mo)))
# genes expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(mo), num_cells_expressed >= 10))
print(head(pData(mo)))

# If you are using RPC values to measure expression, as we are in this vignette,
# it's also good to look at the distribution of mRNA totals across the cells:
pData(mo)$Total_mRNAs <- Matrix::colSums(exprs(mo))

mo <- mo[,pData(mo)$Total_mRNAs < 1e6]

pdf(paste0(outGraph, "nUMI_density.pdf"))
qplot(Total_mRNAs, data = pData(mo), color = as.factor(BRAIN), geom = "density")
qplot(Total_mRNAs, data = pData(mo), color = REGION, geom = "density")
dev.off()

print("Size factors and dispersion")
# Size factors and dispersion
mo <- estimateSizeFactors(mo) 
mo <- estimateDispersions(mo)

print("Dispersed genes to use for pseudotime ordering")
# Dispersed genes to use for pseudotime ordering
disp_table <- dispersionTable(mo) 
ordering_genes <- subset(disp_table
  , mean_expression >= 0.01 & dispersion_empirical >= 0.25 * dispersion_fit)$gene_id
mo_filtered <- setOrderingFilter(mo, ordering_genes) 
# Plot
png(paste0(outGraph, "OrderingGenesDispersion.png"))
plot_ordering_genes(mo_filtered)
dev.off()

print("Variance explained by each PC")
# Variance explained by each PC
png(paste0(outGraph, "PCA_VarianceExplained.png"))
plot_pc_variance_explained(mo_filtered, verbose = TRUE
  , use_existing_pc_variance = TRUE, return_all = FALSE) 
dev.off()

print("Reduce data dimensionality")
# Reduce data dimensionality
# Use number of genes expressed or total mRNAs?
mo_filtered <- reduceDimension(mo_filtered, max_components = 40,
  residualModelFormulaStr = "~individual + librarylab + Total_mRNAs"
  , verbose = TRUE)

print("Order cells along trajectory")
# Order cells along trajectory
mo_filtered <- orderCells(mo_filtered)

save(mo, mo_filtered, file = paste0(outRdat, "monocleO.Robj"))
################################################################################

### Pseudotime gene module correlations

# Pseudotime gene module heatmaps
diff_test_res <- differentialGeneTest(mo_filtered
  , fullModelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
  , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")

save(mo, mo_filtered, diff_test_res, file = paste0(outRdat, "monocleO.Robj"))

# Order by qval
diff_test_res <- diff_test_res[order(diff_test_res$qval), ]

sig_gene_names <- row.names(diff_test_res)[
  diff_test_res$use_for_ordering == TRUE][1:100]
png(paste0(outGraph, "Pseudotime_Heatmap_Cluster", cluster, ".png")
  , width = 12, height = 16, units = "in", res = 300)
plot_pseudotime_heatmap(mo_filtered[sig_gene_names,],
  num_clusters = 3,
  cores = 1,
  show_rownames = T)
dev.off()
################################################################################

### Plots

## Trajectory colored by state or covariates

# Visualize trajectory in reduced dimensional space
ggL <- list(
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "REGION", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "REGION", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "State", cell_size = 0.01) +
    theme(legend.position = "none") + ggtitle("Colored by state")
  , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "State", cell_size = 0.01) +
    theme(legend.position = "none") + ggtitle("Colored by state")
  , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "BRAIN", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "BRAIN", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 2, color_by = "LIBRARY", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "LIBRARY", cell_size = 0.01)
)
# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 5)))})

# Plot grid - 4 columns
pg <- plot_grid(plotlist = ggL, ncol = 4)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by state and covariates"
  , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "Trajectory_ncol4.png")
  , width = 13, height = 8)

# Plot grid - 2 columns
pg <- plot_grid(plotlist = ggL, ncol = 2)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by state and covariates"
  , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "Trajectory_ncol2.png")
  , width = 13, height = 26)

# Plot faceted by state
plot_cell_trajectory(mo_filtered, 1, 2, color_by = "State", cell_size = 0.01) +
  facet_wrap(~State, ncol = 4) +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMonocle trajectory faceted by state"
    , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
    , "\n"))
ggsave(paste0(outGraph, "Trajectory_Facet_Comp12.png")
  , width = 13, height = 26)
# Plot faceted by state
plot_cell_trajectory(mo_filtered, 1, 3, color_by = "State", cell_size = 0.01) +
  facet_wrap(~State, ncol = 4) +
  theme(legend.position = "none") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMonocle trajectory faceted by state"
    , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
    , "\n"))
ggsave(paste0(outGraph, "Trajectory_Facet_Comp13.png")
  , width = 13, height = 26)


## Trajectory colored by Seurat clusters

# Set Seurat cluster factor levels
mo_filtered@phenoData$cluster <- factor(mo_filtered@phenoData$cluster
  , levels = sort(as.numeric(as.character(unique(mo_filtered@phenotype$cluster)))))

# Plot Seurat clusters on trajectory
ggL <- list(
  plot_cell_trajectory(mo_filtered, 1, 2, color_by = "cluster", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 3, color_by = "cluster", cell_size = 0.01)
  , plot_cell_trajectory(mo_filtered, 1, 4, color_by = "cluster", cell_size = 0.01)
)
# gg <- l[[1]]
# gg <- ggplot(gg$data, aes(x = data_dim_1, y = data_dim_2
#   , shape = cluster, color = cluster)) +
#   scale_shape_manual(values = 1:nlevels(gg$data$cluster)) +
#   geom_point() + 
#   theme(legend.position = "right")
# gg$layer[[2]] <- l[[1]]$layer[[1]]

# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 7)))
  })
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 1, align = 'h', axis = 't')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(3, 2))
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by Seurat clusters"
  , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# Save
ggsave(paste0(outGraph, "Trajectory_SeuratCluster.png"), width = 13, height = 20)

# Plot Seurat clusters on trajectory split as grid component 1 vs 2
ggL <- lapply(sort(as.numeric(as.character(unique(mo_filtered@phenoData$cluster))))
  , function(cluster) {
    
    mo_filtered@phenoData$clusterX <- FALSE
    mo_filtered@phenoData$clusterX[mo_filtered@phenoData$cluster %in% cluster] <- TRUE
    
    # Plot Seurat clusters on trajectory
    gg <- plot_cell_trajectory(mo_filtered, 1, 2, color_by = "clusterX", cell_size = 0.001) +
      scale_colour_manual(name = "Cluster:"
        , values = c("#a6cee3", "red")) +
      ggtitle(paste0("Cluster: ", cluster)) +
      theme(legend.position = "right")
    return(gg)
})
# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 7)))
})
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 't')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.2))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by Seurat clusters"
  , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
# Save
ggsave(paste0(outGraph, "Trajectory_SeuratClusterGrid_Comp12.png"), width = 16
  , height = 14)

# Plot Seurat clusters on trajectory split as grid component 1 vs 3
ggL <- lapply(sort(as.numeric(as.character(unique(mo_filtered@phenoData$cluster))))
  , function(cluster) {
    
    mo_filtered@phenoData$clusterX <- FALSE
    mo_filtered@phenoData$clusterX[mo_filtered@phenoData$cluster %in% cluster] <- TRUE
    
    # Plot Seurat clusters on trajectory
    gg <- plot_cell_trajectory(mo_filtered, 1, 3, color_by = "clusterX", cell_size = 0.001) +
      scale_colour_manual(name = "Cluster:"
        , values = c("#a6cee3", "red")) +
      ggtitle(paste0("Cluster: ", cluster)) +
      theme(legend.position = "right")
    return(gg)
  })
# Change legend size
ggL <- lapply(ggL, function(gg) {
  gg + guides(colour = guide_legend(override.aes = list(size = 7)))
})
# extract the legend from one of the plots
legend <- get_legend(ggL[[1]])
# Remove legends from plots
ggL <- lapply(ggL, function(gg) {gg + theme(legend.position = "none")})
# plot_grid combine tSNE graphs
pg <- plot_grid(plotlist = ggL, ncol = 4, align = 'h', axis = 't')
# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pg <- plot_grid(pg, legend, rel_widths = c(1, 0.2))
# now add the title
title <- ggdraw() + draw_label(paste0(graphCodeTitle
  , "\n\nMonocle trajectory colored by Seurat clusters"
  , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"))
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.15, 1))
# Save
ggsave(paste0(outGraph, "Trajectory_SeuratClusterGrid_Comp13.png"), width = 16
  , height = 14)

# ## Trajectory colored by pseudotime
# ggL <- list(
#   plot_cell_trajectory(mo_filtered, 1, 2, color_by = "Pseudotime"
#     , cell_size = 0.01) +
#     scale_color_viridis() +
#     theme(legend.position = "right")
#   plot_cell_trajectory(mo_filtered, 1, 3, color_by = "Pseudotime"
#     , cell_size = 0.01) +
#     scale_color_viridis() +
#     theme(legend.position = "right")
#   plot_cell_trajectory(mo_filtered, 1, 4, color_by = "Pseudotime"
#     , cell_size = 0.01) +
#     scale_color_viridis() +
#     theme(legend.position = "right")
# )
# # plot_grid combine tSNE graphs
# pg <- plot_grid(plotlist = ggL, ncol = 1, align = 'h', axis = 't')
# # now add the title
# title <- paste0(graphCodeTitle
#   , "\n\nMonocle trajectory colored by pseudotime"
#   , "\nCells from progenitor, IPC, and excitatory neuron Seurat clusters"
#   , "\n")
# title <- ggdraw() + draw_label(title)
# # rel_heights values control title margins
# pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
# # Save
# ggsave(paste0(outGraph, "Trajectory_Pseudotime.png"), width = 13, height = 20)
################################################################################

## Luis markers

# Subset to marker genes of interest for Luis' excel file
# Cleanup marker data frame
kmDF <- kmDF[! kmDF$Gene.Symbol == "", ]
kmDF <- kmDF[! is.na(kmDF$Grouping), ]
kmDF$Grouping <- factor(kmDF$Grouping, levels = unique(kmDF$Grouping))
kmDFL <- split(kmDF, kmDF$Grouping)

# Plot Luis markers - Seurat norm center scaled expression
ggL <- lapply(kmDFL, function(kmDF) {
  genes <- kmDF$Gene.Symbol
  print(genes)
  gg <- Plot_Trajectory_Gene_Expression(monocleO = mo_filtered
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
  , "\n\nMean expression of groups of published marker genes"
  , "\nNormalized mean centered variance scaled"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.2, 1))
ggsave(paste0(outGraph, "Trajectory_LuisMarkers_NormCentScale.png")
  , width = 12, height = 6+length(ggL)/2)


## Molyneaux markers

# Molyneaux
# Seurat normalized centered scaled expression
ggL <- apply(mmDF, 1, function(v1) {
  genes <- v1[["hgnc_symbol"]]
  print(genes)
  Plot_Trajectory_Gene_Expression(monocleO = mo_filtered
    , exprColumn = "EXPRESSION", genes = genes
    # , exprM = as.matrix(exprs(mo_filtered))
    , exprM = as.matrix(centSO@scale.data)
    , limHigh = 1.5, limLow = -1.5, title = genes)
})
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of Molyneaux layer genes"
  , "\nNormalized mean centered variance scaled"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(outGraph, "Trajectory_Molyneaux_NormCentScale.png")
  , width = 12, height = length(ggL), limitsize = FALSE)

# Monocle normalized expression
# Tmp expression matrix from monocle object
m <- as.matrix(exprs(mo_filtered))
# Loop through genes
ggL <- apply(mmDF, 1, function(v1) {
  genes <- v1[["hgnc_symbol"]]
  print(genes)
  Plot_Trajectory_Gene_Expression(monocleO = mo_filtered
    , exprColumn = "EXPRESSION", genes = genes
    # , exprM = as.matrix(exprs(mo_filtered))
    , exprM = m
    , limHigh = 1.5, limLow = -1.5, title = genes)
})
pg <- plot_grid(plotlist = ggL, ncol = 3)
# now add the title
title <- paste0(graphCodeTitle
  , "\n\nExpression of Molyneaux layer genes"
  , "\nNormalized mean centered variance scaled"
  , "\n")
title <- ggdraw() + draw_label(title)
# rel_heights values control title margins
pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.05, 1))
ggsave(paste0(
  outGraph, "Trajectory_Molyneaux_MonocleNormExpr.png")
  , width = 12, height = length(ggL)*1.25, limitsize = FALSE)
# Delete tmp expression matrix
rm(m)
################################################################################


# # modules of genes that co-vary across pseudotime
# diff_test_res <- differentialGeneTest(mo_filtered,
#   fullModelFormulaStr = "~sm.ns(Pseudotime)")
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
# plot_pseudotime_heatmap(mo_filtered,
#   num_clusters = 3,
#   cores = 1,
#   show_rownames = T)