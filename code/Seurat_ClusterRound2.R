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
require(irlba)
require(gridExtra)
require(cowplot)
require(viridis)
require(tidyverse)
require(reshape2)
source("Function_Library.R")
source("GGplot_Theme.R")
# require(xlsx)

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
# args <- 6
print(args)

## Inputs

# Seurat clustering object
# PC 1-40
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM
clustersL <- sort(as.numeric(as.character(unique(centSO@ident))))
clusterID <- clustersL[[as.numeric(args[1])]]
print(paste0("Cluster ID: ", clusterID))

## Variables
cluster_annot_tb <- tribble(
    ~cluster_number, ~cluster_annot
    , "9",  "vRG"
    , "7",  "oRG"
    , "8",  "PgS"
    , "10", "PgG2M"
    , "2",  "IP"
    , "0",  "ExN"
    , "1",  "ExM"
    , "4",  "ExCal"
    , "3",  "ExDp1"
    , "13", "ExDp2"
    , "5",  "InSST"
    , "6",  "InCALB2"
    , "11", "OPC"
    , "12", "End"
    , "14", "Per"
    , "15", "Mic"
  )
script_name <- "Seurat_ClusterRound2.R"
date <- format(Sys.Date(), "%Y%m%d")
# date <- "20181222"
outGraph <- paste0("../analysis/graphs/Seurat_ClusterRound2/ClusterRound2/"
  , date, "/Seurat_ClusterRound2_")
outData <- paste0(
  "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/"
  , date, "/Seurat_ClusterRound2_")
out_table <- paste0(
  "../analysis/tables/Seurat_ClusterRound2/ClusterRound2/"
  , date, "/Seurat_ClusterRound2_")

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
################################################################################

### Plotting function

main_function <- function(){
  so_and_reg_exm_l <- run_seurat()
  so <- so_and_reg_exm_l[[1]]
  noCentExM <- so_and_reg_exm_l[[2]]
  rd1CentExM <- so_and_reg_exm_l[[3]]
  rm(so_and_reg_exm_l)
  so <- seurat_cluster_and_tsne(
    so = so, pcs = list(c(1:5), c(1:6), c(1:8), c(1:10), c(1:40))
    , resolutions = c(0.4, 0.5, 0.54, 0.6, 0.7, 0.8))
  save(so, noCentExM, rd1CentExM
    , file = paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))
}

main_plotting_function <- function(){
  load(paste0("../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/"
  , "20190103", "/Seurat_ClusterRound2_", "Cluster", clusterID, "_seuratO.Robj"))
  # load(paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))
  plot_pcs_variance(so = so)
  plot_genes_with_highest_pc_loadings(so = so)
  plot_seurat_pcs_tsne_cluster_tests(so = so)
  plot_tsnes_colored_by_covariates(so = so)
  plot_40k_cell_tsne_colored_by_subclustering(so = so)
  plot_40k_cell_tsne_subclustered(so = so)
  plot_genes_with_highest_pc_loadings_expression_heatmap(
    so = so, seurat_40k_obj = centSO, cluster_col = "clusters_pc1to5_res_0.7"
    , cluster_id = clusterID)
  plot_genes_with_highest_pc_loadings_expression_heatmap(
    so = so, seurat_40k_obj = centSO, cluster_col = "clusters_pc1to10_res_0.5"
    , cluster_id = clusterID)
  plot_tsne_colored_by_pc_scores(so = so, tsne_slot = "tsne_pc1to5")
  plot_tsne_colored_by_pc_scores(so = so, tsne_slot = "tsne_pc1to10")
}
################################################################################

### Run Seurat

run_seurat <- function(){
  print("run_seurat")

  ## Split by cluster

  ## Select Seurat cluster to re-cluster
  # List of cluster IDs
  # clustersL <- append(list(c(0, 1, 4), c(0, 1), c(3, 13), c(5, 6), c(7, 9))
  #   , sort(as.numeric(as.character(unique(centSO@ident)))))
  # clustersL <- sort(as.numeric(as.character(unique(centSO@ident))))
  # clusterID <- clustersL[[as.numeric(args[1])]]
  # print(paste0("Cluster ID: ", clusterID))

  Subset_Seurat_Raw_Counts_By_Cluster <- function(clusterID, seuratO){
    print("Subset_Seurat_Raw_Counts_By_Cluster")
    exDF <- centSO@raw.data
    clusterIDs_cellIDs <- centSO@ident
    clusterIDs_cellIDs <- clusterIDs_cellIDs[clusterIDs_cellIDs %in% clusterID]
    exDF <- exDF[, colnames(exDF) %in% names(clusterIDs_cellIDs)]
    return(exDF)
  }

  Create_Subset_Seurat_Object <- function(seuratO, clusterID){
    print("Create_Subset_Seurat_Object")
    exDF <- Subset_Seurat_Raw_Counts_By_Cluster(
      clusterID = clusterID, seuratO = seuratO
    )
    so <- CreateSeuratObject(raw.data = exDF
      , min.cells = 0, min.genes = 0
      # (ln (transcripts-per-10,000 + 1))
      , normalization.method = "LogNormalize", scale.factor = 10000
      , project = paste0("Cluster ", clusterID) , do.scale = FALSE, do.center = FALSE)
    # Add metadata
    metDF <- centSO@meta.data
    ssMetDF <- metDF[row.names(metDF) %in% colnames(so@raw.data), ]
    so <- AddMetaData(so, metadata = ssMetDF)
    return(so)
  }

  # Initialize Seurat object for each cluster
  so <- Create_Subset_Seurat_Object(seuratO = centSO, clusterID = clusterID)

  # Save centered scaled expression from Seurat round 1
  rd1CentExM <- centSO@scale.data[
    , colnames(centSO@scale.data) %in% colnames(so@raw.data)]
  idx <- match(colnames(so@raw.data), colnames(rd1CentExM))
  rd1CentExM <- rd1CentExM[ ,idx]

  ## Detection of variable genes across the single cells

  print("### Detection of variable genes across the single cells")

  so <- FindVariableGenes(so, mean.function = ExpMean
      , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3
      , y.cutoff = 0.5
      )

  # Check length of var.genes
  print(paste0("Cluster: ", clusterID, " length var.genes: "
   , length(so@var.genes)
  ))

  ## Regress out unwanted sources of variation

  print("### Regress out unwanted sources of variation")

  # note that this overwrites so@scale.data. Therefore, if you intend to use
  # ScaleData, you can set do.scale=F and do.center=F in the original object to
  # save some time.

  ## Regress covariates and center scale

  so <- ScaleData(so, vars.to.regress = c("nUMI", "librarylab", "individual"))

  ## No center or scale

  # Make dataframe of covariates to regress out
  covDF <- data.frame(nUMI = so@meta.data$nUMI
    , librarylab = so@meta.data$librarylab
    , individual = so@meta.data$individual)
  # Regress out confounding variables
  RegressCovariates <- function (exM, covDF) {
    exRegCovM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
    rownames(exRegCovM) <- rownames(exM)
    colnames(exRegCovM) <- colnames(exM)
    # ncol(covDF)+1 when condition has 2 levels
    coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
    for (i in 1:nrow(exM)) {
      if (i%%1000 == 0) {print(paste("Done for", i, "genes..."))}
      mod <- lm(as.numeric(exM[i, ]) ~ ., data = covDF)
      # The full data - the undesired covariates
      exRegCovM[i,] <- coef(mod)[1] + mod$residuals
      # lmmod1 <- lm(as.numeric(exM[i, ]) ~ condition + age + sex + pmi, data = covDF)
    }
    return(exRegCovM)
  }

  noCentExM <- RegressCovariates(so@data, covDF)

  ## Perform linear dimensional reduction

  print("### Perform linear dimensional reduction")

  # Perform PCA on the scaled data. By default, the genes in object\@var.genes are
  # used as input, but can be defined using pc.genes. We have typically found that
  # running dimensionality reduction on genes with high-dispersion can improve
  # performance. However, with UMI data - particularly after using ScaleData, we
  # often see that PCA returns similar (albeit slower) results when run on much
  # larger subsets of genes, including the whole transcriptome.
  # r2SO <- PCA(r2SO, pc.genes = r2SO@var.genes, do.print = TRUE
  #   , pcs.print = 5, genes.print = 5)

  # Run PCA with the IRLBA package (iteratively computes the top dimensions,
  # dramatic increase in speed since we only use a fraction of the PCs anyways) if
  # you see the warning "did not converge–results might be invalid!; try
  # increasing maxit or fastpath=FALSE", try increasing maxit
  so <- RunPCA(object = so, pc.genes = so@var.genes, pcs.compute = 50
    , pcs.print = 1:2, genes.print = 5, maxit = 500, weight.by.var = FALSE)

  # ProjectPCA scores each gene in the dataset (including genes not included in
  # the PCA) based on their correlation with the calculated components. Though we
  # don't use this further here, it can be used to identify markers that are
  # strongly correlated with cellular heterogeneity, but may not have passed
  # through variable gene selection. The results of the projected PCA can be
  # explored by setting use.full=T in the functions below.
  # Saves in slot @pca.x.full
  so <- ProjectPCA(so, pcs.print = 0)

  save(so, noCentExM, rd1CentExM
    , file = paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))

  ## Run tSNE and Seurat clustering

  print("### Run tSNE and Seurat clustering")

  ## Clustering parameters to save
  so <- RunTSNE(so, dims.use = 1:10, do.fast = TRUE)
  so <- FindClusters(so, dims.use = 1:10, resolution = 0.5
    , print.output = 0, save.SNN = TRUE)

  save(so, noCentExM, rd1CentExM
    , file = paste0(outData, "Cluster", clusterID, "_seuratO.Robj"))

  return(list(so, noCentExM, rd1CentExM))
}
################################################################################

### Add clustering and tsne with other parameters

seurat_cluster <- function(
  so, resolutions, pcs, cluster_col_name_pfx = "cluster_"){
  print("seurat_cluster")
  for (i in resolutions){
    print(paste0("resolution ", i))
    cluster_col_name <- paste0(cluster_col_name_pfx, "res_", i)
    default_cluster_col_name <- paste0("res.", i)
    #  cluster
    so <- FindClusters(object = so, dims.use = pcs, resolution = i
      , print.output = 0, save.SNN = TRUE)
    # name cluster column in metadata
    so@meta.data[[cluster_col_name]] <-
      so@meta.data[[default_cluster_col_name]]
  }
  return(so)
}

seurat_cluster_and_tsne <- function(
  so, pcs = list(c(1:5), c(1:10)), resolutions = 0.5){

  print("seurat_cluster_and_tsne")

  for (i in pcs){
    print(paste0("pcs: ", paste0(i, collapse = ",")))
    # name for tsne field in seurat object
    reduction_name <- paste0("tsne_pc", min(i), "to", max(i))
    # name for cluster column in metadata
    cluster_col_name <- paste0("clusters_pc", min(i), "to", max(i), "_")
    so <- RunTSNE(
      object = so, dims.use = i, do.fast = TRUE
      , reduction.name = reduction_name)
    so <- seurat_cluster(
      so = so, resolutions = resolutions, pcs = i
      , cluster_col_name_pfx = cluster_col_name)
  }
  return(so)
}
################################################################################

### PCs variance plot

plot_pcs_variance <- function(so){
  print("plot_pcs_variance")
  # PC variance plot
  gg <- DimElbowPlot(so, reduction.type = "pca", dims.plot = 50)
  gg <- gg + ggtitle(
      paste0(script_name
        , "\n\nPC variance for each cluster"
        , "\nCluster: ", clusterID)
      ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(0, max(so@dr$pca@sdev)))
  ggsave(paste0(outGraph, "Cluster", clusterID, "_PCElbowPlot.pdf")
    , width = 5, height = 5)
}
################################################################################

### PCs tests tSNE plots

plot_seurat_pcs_tsne_cluster_tests <- function(so){

  print("plot_seurat_pcs_tsne_cluster_tests")

  clusterings_df <- data.frame(clusterings = colnames(so@meta.data)[
    grep("clusters_", colnames(so@meta.data))])
  clusterings_df$tsnes <- clusterings_df$clusterings
  clusterings_df$tsnes <- gsub("clusters_", "tsne_", clusterings_df$tsnes)
  clusterings_df$tsnes <- gsub("_res_.*", "", clusterings_df$tsnes)
  clusterings_df <- clusterings_df

  gg_l <- apply(clusterings_df, 1, function(x){
    clustering <- x[["clusterings"]]
    tsne <- x[["tsnes"]]
    tryCatch(
      {
        so@dr$tsne <- so@dr[[tsne]]
        cluster_ids <- factor(so@meta.data[[clustering]])
        names(cluster_ids) <- rownames(so@meta.data)
        so@ident <- cluster_ids
        gg <- TSNEPlot(so, pt.size = 0.01, do.return = TRUE, do.label = TRUE)
        gg <- gg + ggplot_set_theme_publication
        gg <- gg + ggtitle(clustering)
        return(gg)
      }
      , error = function(cond) {
        print(paste0("Error for dimensions: ", tsne, " ", clustering))
        message(cond)
        # Choose a return value in case of error
        return(NA)
      }
    )
  })

  gg_l <- gg_l[! is.na(gg_l)]
  # Combine tSNE graphs
  Plot_Grid(gg_l, ncol = 3, rel_height = 0.1
    , title = paste0(script_name
      , "\n\nSeurat cluster and tSNE with different PCs used"
      , "\nCluster: ", clusterID)
  )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_ClusterTests.png")
    , width = 12, height = 4+1*nrow(clusterings_df), limitsize = FALSE)
}
################################################################################
#
# ### Clustering tests tsne plots
#
# plot_seurat_resolution_tests_tsnes <- function(
#   resolutions = c(0.3, 0.4, 0.5, 0.54, 0.6, 0.7)){
#
#   print("plot_seurat_resolution_tests_tsnes")
#
#   ## Resolution test
#   Seurat_Resolution_Test <- function(seuratO, resolutions, dims_use){
#     print("Seurat_Resolution_Test")
#     so <- tryCatch({
#       RunTSNE(seuratO, dims.use = dims_use, do.fast = TRUE)
#       } ,
#       error = function(cond) {
#         print(paste0("Error for dimensions: ", dims_use))
#         message(cond)
#         # Choose a return value in case of error
#         return(NA)
#       }
#     )
#     ggL <- lapply(resolutions, function(resolution) {
#       tryCatch({
#         print(paste0("Resolution: ", resolution))
#         # so <- FindClusters(so, dims.use = dims_use, resolution = resolution
#         #   , print.output = 0, save.SNN = TRUE)
#         res_col_name <- paste0("res.", resolution)
#         cluster_ids <- as.factor(so@meta.data[[res_col_name]])
#         names(cluster_ids) <- rownames(so@meta.data)
#         so@ident <- cluster_ids
#         gg <- TSNEPlot(so, pt.size = 0.01, do.return = TRUE, do.label = TRUE)
#         gg <- gg + ggtitle(paste0("Resolution: ", resolution))
#         gg <- gg + ggplot_set_theme_publication
#         return(gg)
#       } ,
#       error = function(cond) {
#         print(paste0("Error for resolution: ", resolution))
#         message(cond)
#         # Choose a return value in case of error
#         return(NA)
#       }
#     )
#     })
#     return(ggL)
#   }
#
#   # Plot using PC1-40
#   ggL <- Seurat_Resolution_Test(
#     seuratO = so
#     , dims_use = 1:40
#     , resolutions = resolutions
#   )
#   ggL <- ggL[! is.na(ggL)]
#   # Combine tSNE graphs
#   pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
#     , title = paste0(script_name
#       , "\n\nSeurat cluster and tSNE with different resolutions"
#       , "\nCluster: ", clusterID
#       , "\nPC 1-40"
#     )
#   )
#   ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-40_ResolutionTest.png")
#     , width = 12, height = 7)
#
#   # Plot using PC1-15
#   ggL <- Seurat_Resolution_Test(
#     seuratO = so
#     , dims_use = 1:15
#     , resolutions = resolutions
#   )
#   ggL <- ggL[! is.na(ggL)]
#   # Combine tSNE graphs
#   pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
#     , title = paste0(script_name
#       , "\n\nSeurat cluster and tSNE with different resolutions"
#       , "\nCluster: ", clusterID
#       , "\nPC 1-10"
#     )
#   )
#   ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-15_ResolutionTest.png")
#     , width = 12, height = 7)
#
#   # Plot using PC1-10
#   ggL <- Seurat_Resolution_Test(
#     seuratO = so
#     , dims_use = 1:10
#     , resolutions = resolutions
#   )
#   ggL <- ggL[! is.na(ggL)]
#   # Combine tSNE graphs
#   pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
#     , title = paste0(script_name
#       , "\n\nSeurat cluster and tSNE with different resolutions"
#       , "\nCluster: ", clusterID
#       , "\nPC 1-10"
#     )
#   )
#   ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-10_ResolutionTest.png")
#     , width = 12, height = 7)
#
#   # Plot using PC1-8
#   ggL <- Seurat_Resolution_Test(
#     seuratO = so
#     , dims_use = 1:8
#     , resolutions = resolutions
#   )
#   ggL <- ggL[! is.na(ggL)]
#   # Combine tSNE graphs
#   pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
#     , title = paste0(script_name
#       , "\n\nSeurat cluster and tSNE with different resolutions"
#       , "\nCluster: ", clusterID
#       , "\nPC 1-10"
#     )
#   )
#   ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-8_ResolutionTest.png")
#     , width = 12, height = 7)
#
#   # Plot using PC1-5
#   ggL <- Seurat_Resolution_Test(
#     seuratO = so
#     , dims_use = 1:5
#     , resolutions = resolutions
#   )
#   ggL <- ggL[! is.na(ggL)]
#   # Combine tSNE graphs
#   pg <- Plot_Grid(ggL, ncol = 3, rel_height = 0.3
#     , title = paste0(script_name
#       , "\n\nSeurat cluster and tSNE with different resolutions"
#       , "\nCluster: ", clusterID
#       , "\nPC 1-5"
#     )
#   )
#   ggsave(paste0(outGraph, "Cluster", clusterID, "_tSNE_PC1-5_ResolutionTest.png")
#     , width = 12, height = 7)
#
# }
################################################################################

### Feature plot of covariates

plot_tsnes_colored_by_covariates <- function(so){

  print("plot_tsnes_colored_by_covariates")

  Feature_Plot_Of_Covariate <- function(seuratO, covariate){
    ("Feature_Plot_Of_Covariate")
    # Collect tSNE values
    ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
    # Add cluster identity
    ggDF$CLUSTER <- so@ident
    # Add metadata
    ggDF <- data.frame(ggDF, so@meta.data[match(row.names(ggDF), so@meta.data$CELL), ])
    ggDF$BRAIN <- as.factor(ggDF$BRAIN)
    # Plot
    gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF[[covariate]])) +
      geom_point(size = 0.1, alpha = 0.5) +
      guides(colour = guide_legend(override.aes = list(size = 7))) +
      ggplot_set_theme_publication +
      ggtitle(covariate)
    if(class(ggDF[[covariate]]) == "numeric"){gg <- gg + scale_color_viridis()}
    return(gg)
  }

  ## PC 1-10
  # tSNE colored by clustering
  ggTsne <- TSNE_Plot(so) +
    theme(legend.position = "none") +
    ggplot_set_theme_publication
  # tSNE colored by covariate
  gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
  gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
  gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
  gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
  gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
  # plot grid
  pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
    , align = "v", axis = "r", ncol = 3, rel_height = 0.3
    , title = paste0(script_name
      , "\n\nSeurat clustering round 2 and covariates"
      , "\nCluster: ", clusterID
      , "\nSeurat variable genes used for clustering"
      , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-10 resolution 0.5"
      , "\n")
    )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-10_covariates.png")
    , width = 14, height = 7, limitsize = FALSE)

  ## PC 1-15
  # Re do tSNE and clustering
  so <- RunTSNE(so, dims.use = 1:15, do.fast = TRUE)
  so <- FindClusters(so, dims.use = 1:15, resolution = 0.5
  , print.output = 0, save.SNN = TRUE)
  # tSNE colored by clustering
  ggTsne <- TSNE_Plot(so) +
    theme(legend.position = "none") +
    ggplot_set_theme_publication
  # tSNE colored by covariate
  gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
  gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
  gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
  gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
  gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
  # plot grid
  pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
    , align = "v", axis = "r", ncol = 3, rel_height = 0.3
    , title = paste0(script_name
      , "\n\nSeurat clustering round 2 and covariates"
      , "\nCluster: ", clusterID
      , "\nSeurat variable genes used for clustering"
      , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-15 resolution 0.5"
      , "\n")
    )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-15_covariates.png")
    , width = 14, height = 7, limitsize = FALSE)

  ## PC 1-5
  # Re do tSNE and clustering
  so <- RunTSNE(so, dims.use = 1:5, do.fast = TRUE)
  so <- FindClusters(so, dims.use = 1:5, resolution = 0.5
  , print.output = 0, save.SNN = TRUE)
  # tSNE colored by clustering
  ggTsne <- TSNE_Plot(so) +
    theme(legend.position = "none") +
    ggplot_set_theme_publication
  # tSNE colored by covariate
  gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
  gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
  gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
  gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
  gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
  # plot grid
  pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
    , align = "v", axis = "r", ncol = 3, rel_height = 0.3
    , title = paste0(script_name
      , "\n\nSeurat clustering round 2 and covariates"
      , "\nCluster: ", clusterID
      , "\nSeurat variable genes used for clustering"
      , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-5 resolution 0.5"
      , "\n")
    )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-5_covariates.png")
    , width = 14, height = 7, limitsize = FALSE)

  ## PC 1-40
  # Re do tSNE and clustering
  so <- RunTSNE(so, dims.use = 1:40, do.fast = TRUE)
  so <- FindClusters(so, dims.use = 1:40, resolution = 0.5
  , print.output = 0, save.SNN = TRUE)
  # tSNE colored by clustering
  ggTsne <- TSNE_Plot(so) +
    theme(legend.position = "none") +
    ggplot_set_theme_publication
  # tSNE colored by covariate
  gg1 <- Feature_Plot_Of_Covariate(so, "BRAIN")
  gg2 <- Feature_Plot_Of_Covariate(so, "LIBRARY")
  gg3 <- Feature_Plot_Of_Covariate(so, "REGION")
  gg4 <- Feature_Plot_Of_Covariate(so, "nGene")
  gg5 <- Feature_Plot_Of_Covariate(so, "nUMI")
  # plot grid
  pg <- Plot_Grid(list(ggTsne, gg1, gg2, gg3, gg4, gg5)
    , align = "v", axis = "r", ncol = 3, rel_height = 0.3
    , title = paste0(script_name
      , "\n\nSeurat clustering round 2 and covariates"
      , "\nCluster: ", clusterID
      , "\nSeurat variable genes used for clustering"
      , "\ntSNE PC 1-40, cluster round 2 tSNE PC 1-40 resolution 0.5"
      , "\n")
    )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_PC1-40_covariates.png")
    , width = 14, height = 7, limitsize = FALSE)
}
################################################################################

### Plot genes with highest PC loadings

plot_genes_with_highest_pc_loadings <- function(so){

  print("plot_genes_with_highest_pc_loadings")

  # Plot genes with highest PC loadings
  ggL <- lapply(1:10, function(pc) {
    df <- rbind(data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc])[1:10])
    , data.frame(PC = sort(so@dr$pca@gene.loadings[ ,pc], decreasing = TRUE)[10:1])
      )
    df$GENE <- factor(row.names(df), levels = row.names(df))
    gg <- ggplot(df, aes(x = PC, y = GENE)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      ggtitle(paste0("Cluster: ", clusterID, "\nPC ", pc))
    return(gg)
  })
  Plot_Grid(ggL, ncol = 5, align = 'v', axis = 'r', rel_height = 0.2
    , title = paste0(script_name
      , "\n\nGenes with highest PC loadings"
      , "\nCluster: ", clusterID)
  )
  ggsave(paste0(outGraph, "cluster", clusterID, "_PCAplots.pdf"), width = 15
    , height = 12, limitsize = FALSE)
}
################################################################################

### tSNE plot colored by top PC scores

plot_tsne_colored_by_pc_scores <- function(so = so, tsne_slot){

  print("plot_tsne_colored_by_pc_scores")

  so@dr$tsne <- so@dr[[tsne_slot]]

  # Feature plots of top PC scores
  ggL <- lapply(c(1:10), function(i) {
    # Collect tSNE values
    ggDF <- as.data.frame(so@dr$tsne@cell.embeddings)
    # Add PC score
    ggDF$PC <- so@dr$pca@cell.embeddings[ ,i]
    gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = PC)) +
      geom_point(size = 1, alpha = 0.5) +
      # guides(colour = guide_legend(override.aes = list(size = 7))) +
      scale_color_distiller(name = "PC score", type = "div"
        , palette = 5, direction = -1) +
      ggtitle(paste0("PC: ", i)) +
      ggplot_set_theme_publication
    return(gg)
  })
  # plot grid
  pg <- plot_grid(plotlist = ggL, align = "v", axis = "l", ncol = 3)
  # now add the title
  title <- ggdraw() + draw_label(paste0(script_name
    , "\n\nSeurat cluster and tSNE colored by PC scores"
    , "\nCluster: ", clusterID
    , "\n", tsne_slot
    , "\n"))
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  # save
  ggsave(paste0(outGraph, "cluster", clusterID, "_PCscore_", tsne_slot, ".png")
  , width = 15, height = 15)
}
################################################################################

### Color 40k cell tSNE by sub-clustering

plot_40k_cell_tsne_colored_by_subclustering <- function(so){

  print("plot_40k_cell_tsne_colored_by_subclustering")

  # Collect tSNE values
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Add round 1 cluster identity
  ggDF$Cluster <- centSO@ident

  # plot round 2 clusters
  clusterings <- names(so@meta.data)[grep("clusters_pc", names(so@meta.data))]
  gg_l <- lapply(clusterings, function(clustering){
    # Add round 2 cluster identity
    idx <- match(rownames(ggDF), rownames(so@meta.data))
    ggDF$Cluster_Round2 <- so@meta.data[[clustering]][idx]
    # Plot
    ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster_Round2)) +
      geom_point(size = 0.02, alpha = 0.5) +
      guides(colour = guide_legend(override.aes = list(size = 3)
        , title = "Cluster")) +
      ggplot_set_theme_publication +
      ggtitle(clustering)
  })

  # plot round 1 clusters
  tnse_r1_gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster)) +
    geom_point(size = 0.02, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2
      , title = "Cluster")) +
    ggplot_set_theme_publication +
    ggtitle("Clustering round 1")

  # add tsne colored by round 1 to round 2 clusterings
  Plot_Grid(append(list(tnse_r1_gg), gg_l)
    , ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(script_name
      , "\n\nSeurat tSNE round 1 colored by round 2 clustering"
      , "\nCluster: ", clusterID
      , "\n")
  )
  ggsave(paste0(outGraph, "cluster", clusterID, "_tSNE_R1_Cluster_R1.png")
    , width = 14, height = 4+2.5*length(clusterings), limitsize = FALSE)
}
################################################################################

### Plot 40k cell tSNE sub-clustered

plot_40k_cell_tsne_subclustered <- function(so){

  print("plot_40k_cell_tsne_subclustered")

  # Collect tSNE values
  ggDF <- as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Add round 1 cluster identity
  ggDF$Cluster <- centSO@ident

  # plot round 2 clusters
  clusterings <- names(so@meta.data)[grep("clusters_pc", names(so@meta.data))]
  gg_l <- lapply(clusterings, function(clustering){
    # Add round 2 cluster identity
    idx <- match(rownames(ggDF), rownames(so@meta.data))
    ggDF$Cluster_Round2 <- so@meta.data[[clustering]][idx]
    # Remove cells not in cluster
    ggDF <- ggDF[! is.na(ggDF$Cluster_Round2), ]
    # Plot
    ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster_Round2)) +
      geom_point(size = 0.02, alpha = 0.5) +
      guides(colour = guide_legend(override.aes = list(size = 3)
        , title = "Cluster")) +
      ggplot_set_theme_publication +
      ggtitle(clustering)
  })

  # plot round 1 clusters
  tnse_r1_gg <- ggplot(ggDF, aes(x = tSNE_1, y = tSNE_2, col = ggDF$Cluster)) +
    geom_point(size = 0.02, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2
      , title = "Cluster")) +
    ggplot_set_theme_publication +
    ggtitle("Clustering round 1")

  # add tsne colored by round 1 to round 2 clusterings
  Plot_Grid(append(list(tnse_r1_gg), gg_l)
    , ncol = 2, rel_height = 0.1, align = 'v', axis = 'r'
    , title = paste0(script_name
      , "\n\nSeurat tSNE round 1 colored by round 2 clustering"
      , "\nCluster: ", clusterID
      , "\n")
  )
  ggsave(paste0(
      outGraph, "cluster", clusterID, "_tSNE_R1_subcluster_subset.png")
    , width = 10, height = 4+1.75*length(clusterings), limitsize = FALSE)
}
################################################################################

### Plot expression heatmap of genes loading highly on PCs of cluster

plot_genes_with_highest_pc_loadings_expression_heatmap <- function(
  so, seurat_40k_obj, cluster_col, cluster_id){

  print("plot_genes_with_highest_pc_loadings_expression_heatmap")

  # genes with high pc loadings
  gene_group_df <- lapply(1:10, function(pc){
    rbind(
      data.frame(loading = sort(so@dr$pca@gene.loadings[ ,pc])[1:10])
      , data.frame(loading = sort(so@dr$pca@gene.loadings[ ,pc], decreasing = TRUE)[1:10])[10:1, ,drop = FALSE]
    ) %>%
    # data.frame(loading = sort(so@dr$pca@gene.loadings[ ,pc]
    #     , decreasing = TRUE)[1:20] %>% rev) %>%
      rownames_to_column(var = "Gene") %>%
      mutate(Group = pc)
    }) %>% do.call("rbind", .)

  # set clustering to use
  cluster_ids <- factor(so@meta.data[[cluster_col]])
  names(cluster_ids) <- rownames(so@meta.data)
  so@ident <- cluster_ids

  # Gather cell IDs and cluster IDs for heatmap
  cellid_clusterid <- seurat_40k_obj@ident %>% enframe(
    name = "cell_id", value = "cluster_ids") %>%
    left_join(
      (so@ident %>% enframe(name = "cell_id", value = "sub_cluster_id"))
      , by = "cell_id") %>%
    mutate_if(is.factor, as.character) %>%
    mutate(cluster_ids = if_else(cluster_ids == cluster_id, paste0(.$cluster_ids, "_", .$sub_cluster_id), cluster_ids))  %>%
    left_join(cluster_annot_tb
      , by = c("cluster_ids" = "cluster_number")) %>%
      select(cell_id, cluster_ids) %>% deframe

  expr_m <- as.matrix(seurat_40k_obj@scale.data)

  cluster_annot <- cluster_annot_tb %>%
    filter(cluster_number == cluster_id) %>%
    pull(cluster_annot)

  cluster_order <- cluster_annot_tb %>%
    left_join(.
        , so@ident %>%
          unique %>%
          sort %>%
          as_tibble %>%
          mutate(cluster = cluster_id %>% as.character) %>%
          mutate(cluster_subcluster = paste0(cluster, "_", value))
      , by = c("cluster_number" = "cluster")) %>%
    mutate(cluster_number = if_else(
      cluster_number == cluster_id, cluster_subcluster, cluster_number)) %>%
    pull(cluster_number)

  graph_height = 2 + nrow(gene_group_df)*0.15

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = expr_m
    , cellID_clusterID = cellid_clusterid
    , clusters = cellid_clusterid %>% unique
    , clusterOrder = cluster_order
    , lowerLimit = -1.5
    , upperLimit = 1.5
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of by sub-cluster of genes loading highly on PCs of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered variance scaled by gene"
  ))
  ggsave(
    paste0(outGraph, "cluster", clusterID, "_", cluster_col
      , "_pc_loading_genes_expression_heatmap_zscore.png")
    , height = graph_height, width = 13, limitsize = FALSE)
}
################################################################################

### make table of cell metadata

make_cell_metadata_table <- function(){

  print("make_cell_metadata_table")

  # compile cell ids and subclusters
  in_subclust_so_l <- list.files(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20190103/"
    , full.names = TRUE)
  subclust_cellid_tb <- in_subclust_so_l %>%
    map(., .f = function(in_subclust_so){
      load(in_subclust_so)
      cluster_number <- so@project.name %>% gsub("Cluster ", "", .) %>% as.numeric
      if(cluster_number %in% c(1:4)){
        tibble(
          Cell = so@meta.data$CELL
          , Subcluster = so@meta.data[["clusters_pc1to10_res_0.5"]]
          , Cluster = cluster_number
          )
      }
      if(cluster_number %in% c(5:15)){
        tibble(
          Cell = so@meta.data$CELL
          , Subcluster = so@meta.data[["clusters_pc1to5_res_0.7"]]
          , Cluster = cluster_number
          )
      }
    }) %>%
    bind_rows() %>%
    mutate(Cluster = as.character(Cluster)) %>%
    left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
    mutate(Subcluster = paste0(cluster_annot, "_", Subcluster)) %>%
    mutate(Cluster = cluster_annot) %>%
    select(-cluster_annot, -Cluster)

  ## Cell metadata
  mdat_paper_DF <- centSO@meta.data
  mdat_paper_DF <- mdat_paper_DF[ ,c("CELL", "res.0.54", "BRAIN", "REGION"
    , "NEXTERA", "LIBRARY", "nGene", "nUMI", "percent.mito"
    , "S.Score", "G2M.Score", "Phase")]
  mdat_paper_DF$Cluster_number <- mdat_paper_DF$"res.0.54"
  # Cluster annotations
  cluster_annot <- c(
    "9" = "vRG"
    , "7" = "oRG"
    , "8" = "PgS"
    , "10" = "PgG2M"
    , "2" = "IP"
    , "0" = "ExN"
    , "1" = "ExM"
    , "4" = "ExCal"
    , "3" = "ExDp1"
    , "13" = "ExDp2"
    , "5" = "InSST"
    , "6" = "InCALB2"
    , "11" = "OPC"
    , "12" = "End"
    , "14" = "Per"
    , "15" = "Mic"
    , "16" = "NA"
  )
  idx <- match(mdat_paper_DF$Cluster_number, names(cluster_annot))
  mdat_paper_DF$Cluster <- cluster_annot[idx]
  # Donor IDs
  donor_annot <- c("2" = "368", "3" = "370", "4" = "371", "5" = "372")
  idx <- match(mdat_paper_DF$BRAIN, names(donor_annot))
  mdat_paper_DF$BRAIN <- donor_annot[idx]
  mdat_paper_DF <- mdat_paper_DF[ ,c(1, 14:13, 3:12)]
  # Gestation week
  gw_annot <- c("368" = "17", "370" = "18", "371" = "17", "372" = "18")
  idx <- match(mdat_paper_DF$BRAIN, names(gw_annot))
  mdat_paper_DF$GW <- gw_annot[idx]
  # Convert to percentage
  mdat_paper_DF$percent.mito <- round(mdat_paper_DF$percent.mito * 100, 2)
  # Order by cluster
  mdat_paper_DF$Cluster_number <- factor(mdat_paper_DF$Cluster_number
    , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15,16))
  mdat_paper_DF <- mdat_paper_DF[order(mdat_paper_DF$Cluster_number), ]
  # Remove cluster 16 cells
  mdat_paper_DF <- mdat_paper_DF[! mdat_paper_DF$Cluster_number == 16, ]
  colnames(mdat_paper_DF) <- c("Cell", "Cluster", "Cluster_number", "Donor"
    , "Layer", "Index", "Library", "Number_genes_detected", "Number_UMI"
    , "Percentage_mitochondrial", "S_phase_score", "G2M_phase_score", "Phase"
    , "Gestation_week")
  # Round
  mdat_paper_DF$S_phase_score <- signif(mdat_paper_DF$S_phase_score, 2)
  mdat_paper_DF$G2M_phase_score <- signif(mdat_paper_DF$G2M_phase_score, 2)
  mdat_paper_DF <- mdat_paper_DF[ ,c(1:2,4:5,14,6:13)]

  # add subcluster ids
  mdat_paper_tb <- left_join(mdat_paper_DF, subclust_cellid_tb
      , by = c("Cell" = "Cell")) %>%
    as_tibble() %>%
    select(Cell, Cluster, Subcluster, everything())

  # output table
  write.csv(mdat_paper_tb, file = paste0(out_table, "cell_metadata.csv")
    , quote = FALSE, row.names = FALSE)
}
################################################################################

### Run

# main_function()
main_plotting_function()
################################################################################
