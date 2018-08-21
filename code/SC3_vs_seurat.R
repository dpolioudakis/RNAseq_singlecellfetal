# Will Connell
# 2017-08-31
# Compare Clustering Tools Full 40k Cells

# module load gcc/4.9.3
# module load R/3.3.0

###########################################################################
# Seurat
require(reshape2)
require(Seurat)
require(Matrix)
require(gridExtra)
# SC3
library(SC3)
library(scater)
library(dplyr)

# Functions ---------------------------------------------------------------

# Jaccard Index: Seurat vs. SC3
Jaccard_Index <- function(v1, v2) {
  sum(v1 %in% v2) / (length(v1) + length(v2) - sum(v1 %in% v2))
}

graphCodeTitle <- c("SC3_vs_seurat")
dir.create(graphCodeTitle)

# Load data ---------------------------------------------------------------

# load Seurat object
load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM
seuratO <- centSO
rm(centSO)
rm(noCentExM)
#First lets stash our identities for later
seuratO <- StashIdent(seuratO, save.name = "SeuratCluster")


# Perform basic SC3 clustering --------------------------------------------

# take raw data of only filtered cells and genes
raw.data <- seuratO@raw.data %>%
  select(colnames(seuratO@scale.data)) %>%
  filter(rownames(.) %in% rownames(seuratO@scale.data))
rownames(raw.data) <- rownames(seuratO@scale.data)

# subsample
# raw.data <- raw.data[ , sample(ncol(raw.data), 5000)]
raw.data <- as.matrix(raw.data)

# cell annotation
metDF <- dplyr::filter(metDF, metDF$CELL %in% colnames(raw.data))
# pd <- new("AnnotatedDataFrame", data = metDF)
# rownames(pd) <- metDF$CELL
# SCESEt object
sceset <- SingleCellExperiment(assays = list(counts = raw.data
  , logcounts = log2(raw.data + 1)),colData = metDF)
rm(raw.data)
# define feature names in feature_symbol column
rowData(sceset)$feature_symbol <- rownames(sceset)
# calculate QC metrics
sceset <- scater::calculateQCMetrics(sceset, use_spikes = FALSE)

# regress out nUMI, donor (brain) and lab (library)
design <- model.matrix(~sceset$nUMI +
                         sceset$BRAIN +
                         sceset$LIBRARY)
# normalize and store expression values (default store in "exprs" which is found in sceset@assayData$exprs
# and used for sc3 clustering method)
#set_exprs(sceset, "exprs_orig") <- exprs(sceset)
metadata(sceset)$log.exprs.offset <- 1
sceset <- scater::normaliseExprs(sceset, design = design,
                 method = "none", exprs_values = "counts",
                 return_norm_as_exprs = TRUE)
# estimate k
#sceset <- sc3_estimate_k(sceset)

# run SC3
sceset <- sc3(sceset, ks = 16:16, biology = FALSE, gene_filter = TRUE
  , svm_num_cells = 15000, svm_max = 5000, rand_seed = 88)
save(sceset, file = sprintf("%s/SC3_vs_seurat.RData", graphCodeTitle))
print(sceset)
print(str(sceset))


# Generate comparison plots -----------------------------------------------

# rename SC3 object
sc3.res <- colData(sceset)
print(str(sc3.res))
meta_data <- seuratO@meta.data
cluster_numb <- 16
for(i in 16){

  clust.df <- data.frame(cell=sc3.res@rownames
    , clust.id=as.numeric(sc3.res[[paste0("sc3_",cluster_numb,"_clusters")]]))
  clust.df <- clust.df[! is.na(clust.df$clust.id), ]
  clust.df$cell <- as.character(clust.df$cell)
  cluster_numb <- cluster_numb + 1

  # check same cells are shared between each object
  print(sprintf("Cell numbers of Seurat object: %d", nrow(seuratO@meta.data)))
  print(sprintf("Cell numbers of SC3 object: %d", nrow(clust.df)))
  print(sprintf('Shared cells total: %d', sum(seuratO@meta.data$CELL %in% clust.df$cell)))

  shared_cells <- inner_join(meta_data, clust.df, by = c('CELL' = 'cell'))
  print(str(shared_cells))
  shared_cells$SC3Cluster <- as.character(shared_cells$clust.id)

  seuratO <- SubsetData(seuratO, cells.use = shared_cells$CELL)
  seuratO@meta.data <- shared_cells
  rownames(seuratO@meta.data) <- shared_cells$CELL
  seuratO@meta.data$SeuratCluster <- as.character(seuratO@meta.data$SeuratCluster)




  # tSNE plot comparisons ---------------------------------------------------
  ## Plot tSNE - NUMBERING
  # tSNE plot - Seurat
  plot1 <- TSNEPlot(seuratO, do.label = TRUE, group.by = "SeuratCluster"
                    , pt.size = 0.2, do.return = TRUE, no.legend = TRUE)
  plot1 <- plot1 + ggtitle(paste0(graphCodeTitle
                                  , "\n"
                                  , "\ntSNE plot, each point is a cell"
                                  , "\nColor indicates cluster assignment"
                                  , "\nSeurat clustering"
                                  , "\n"))
  # tSNE plot - SC3
  plot2 <- TSNEPlot(seuratO, do.label = TRUE, group.by = "SC3Cluster", pt.size = 0.2,
                    do.return = TRUE, no.legend = TRUE)
  plot2 <- plot2 + ggtitle(paste0(graphCodeTitle
                                  , "\n"
                                  , "\ntSNE plot, each point is a cell"
                                  , "\nColor indicates cluster assignment"
                                  , "\nSC3 clustering"
                                  , "\n"))
  # Save
  png(sprintf("%s/clust_%d_tSNE_numbering.png", graphCodeTitle, i), width = 25, height = 12, units = "in", res = 300)
  grid.arrange(plot1, plot2, ncol = 2)
  dev.off()



  ## Plot  tSNE - NO NUMBERING
  # tSNE plot - Seurat
  plot1 <- TSNEPlot(seuratO, do.label = FALSE, group.by = "SeuratCluster"
                    , pt.size = 0.2, do.return = TRUE, no.legend = TRUE)
  plot1 <- plot1 + ggtitle(paste0(graphCodeTitle
                                  , "\n"
                                  , "\ntSNE plot, each point is a cell"
                                  , "\nColor indicates cluster assignment"
                                  , "\nSeurat clustering"
                                  , "\n"))
  # tSNE plot - SC3
  plot2 <- TSNEPlot(seuratO, do.label = FALSE, group.by = "SC3Cluster", pt.size = 0.2,
                    do.return = TRUE, no.legend = TRUE)
  plot2 <- plot2 + ggtitle(paste0(graphCodeTitle
                                  , "\n"
                                  , "\ntSNE plot, each point is a cell"
                                  , "\nColor indicates cluster assignment"
                                  , "\nSC3 clustering"
                                  , "\n"))
  # Save
  png(sprintf("%s/clust_%d_tSNE_no_numbering.png", graphCodeTitle, i), width = 25, height = 12, units = "in", res = 300)
  grid.arrange(plot1, plot2, ncol = 2)
  dev.off()


  # Jaccard Index Heatmaps --------------------------------------------------
  # coerce character columns of clusters to numeric for Jaccard Index function
  seuratO@meta.data$SeuratCluster <- as.numeric(seuratO@meta.data$SeuratCluster)

  # Empty matrix
  # we only want 16 total slots for Seurat (not 17), because Damon would like to
  # eliminate in cluster labels 0-16, the 16th
  jiM <- matrix(NA
                , length(unique(seuratO@meta.data$SeuratCluster))-1
                , length(unique(seuratO@meta.data$SC3Cluster)))

  # Fill with Jaccard index
  # we want order of clusters that matches the biology
  for (p in 1:(length(unique(seuratO@meta.data$SeuratCluster))-1)){
    for (j in 1:length(unique(seuratO@meta.data$SC3Cluster))){
      v1 <- row.names(seuratO@meta.data)[seuratO@meta.data$SeuratCluster == p]
      v2 <- row.names(seuratO@meta.data)[seuratO@meta.data$SC3Cluster == j]
      jiM[p,j] <- Jaccard_Index(v1, v2)
    }
  }
  # Plot Jaccard Index
  colnames(jiM) <- c(1:length(unique(seuratO@meta.data$SC3Cluster)))
  jiM <- as.data.frame(jiM)
  jiM$Seurat <- c(0:15)
  jiM <- jiM %>%
    select(Seurat, everything())
  jiM <- jiM[match(c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15), jiM$Seurat),]
  jiM$Seurat <- as.factor(jiM$Seurat)
  jiM$Seurat <- reorder(jiM$Seurat, 1:16)
  # heatmap gradient
  jiM_melt3 <- melt(jiM, id.vars='Seurat')

  # plot
  pdf(sprintf("%s/clust_%d_jaccard.pdf", graphCodeTitle, i), height = 10, width = 10)
  print(
    ggplot(jiM_melt3, aes(Seurat, variable)) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(jiM_melt3$value, 2))) + # write the values +
      scale_fill_gradient(low = "white",
                          high = "red") +
      theme(plot.title = element_text(hjust = 0.5, size = 25),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20)) +
      ylab(label = "SC3 ") +
      ggtitle('SC3 Clustering Method Comparison', subtitle = paste0("Jaccard Indices between ", nrow(seuratO@meta.data), " cells"))
  )
  dev.off()

}
