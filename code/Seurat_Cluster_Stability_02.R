# Damon Polioudakis
# 2017-03-29
# Clustering stability using bootstrapping + Jaccard method from Christian
# Hennig “Cluster-wise assessment of cluster stability”

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(Seurat)
require(Matrix)
require(fpc)
require(methods)
require(dplyr)
source("Function_Library.R")

## Input
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
rm(noCentExM)
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_centSO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Bootsrapped clustering results
inTables <- list.files(
  "../analysis/tables/Seurat_Cluster_Stability/DS2-11/Iterations"
  , full.names = TRUE
)

## Variables
graphCodeTitle <- "Seurat_Cluster_Stability.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"
outTable <- "../analysis/tables/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"
outTable_iterations <- "../analysis/tables/Seurat_Cluster_Stability/DS2-11/Iterations/Seurat_Cluster_Stability_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_Stability/DS2-11/Seurat_Cluster_Stability_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outTable_iterations), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 10)))
theme_update(plot.title = element_text(size = 10))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Jaccard index to test cluster stability

# Number of bootstraps
nBt <- length(inTables)
# Number of clusters in original
nOrCl <- length(unique(centSO@ident))
# Make matrix to fill in bootstrap jaccard results
jcrdM <- matrix(data = NA, nrow = nOrCl, ncol = nBt)
row.names(jcrdM) <- sort(unique(centSO@ident))
colnames(jcrdM) <- paste("Boot", seq(1, nBt))

# Bootstrap, recluster, + jaccard to test cluster stability
for (i in 1:nBt){
  print(paste("Bootstrap:", i))
  bootstrapped_clusters <- read.csv(inTables[i])
  # Intersection of original and bootstrap samples
  int <- intersect(
    names(centSO@ident)
    # , gsub("*_[0-9]+", "", names(btSO@ident)))
    , as.character(bootstrapped_clusters$Cell_ID)
  )
  # Subset original to intersection samples
  clsOr <- centSO@ident[names(centSO@ident) %in% int]
  # Remove unique integer IDs added to bootstrap samples
  # names(clsSb) <- gsub("*_[0-9]+", "", names(btSO@ident))
  clsSb <- bootstrapped_clusters$Cluster
  names(clsSb) <- bootstrapped_clusters$Cell_ID
  # Subset bootstrap to intersection samples
  clsSb <- clsSb[names(clsSb) %in% int]
  # Remove duplicates
  clsSb <- clsSb[unique(names(clsSb))]
  # Merge original and bootstrap to have same order
  clDF <- merge(data.frame(clsOr), data.frame(clsSb), by = "row.names")
  # Make TRUE FALSE list of cluster membership for clujaccard() (fpc package)
  # Original
  origCl <- lapply(unique(clDF$clsOr), function(cl) {clDF$clsOr == cl})
  names(origCl) <- unique(clDF$clsOr)
  # Bootstrap
  btClL <- lapply(unique(clDF$clsSb), function(cl) {clDF$clsSb == cl})
  # Max jaccard index for each original cluster and bootstrap clusters
  jcrd_L <- lapply(names(origCl), function(name) {
    oCl <- origCl[[name]]
    jaccard <- max(sapply(btClL, function(btCl) {clujaccard(oCl, btCl)}))
    data.frame(name, jaccard)
  })
  jcrd_DF <- do.call("rbind", jcrd_L)
  # Add as column of matrix (rows are clusters)
  for(j in 1:nrow(jcrd_DF)){
    cluster <- as.numeric(jcrd_DF$name[j])
    print(cluster)
    jcrdM[cluster,i] <- jcrd_DF$jaccard[j]
  }
}

# Mean jaccard index for each cluster
cluster_reorder <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15,16)
mean_jcrd_DF <- data.frame(
  Cluster = factor(names(rowMeans(jcrdM)), levels = names(rowMeans(jcrdM)))
  , Cluster_Reorder = factor(cluster_reorder, levels = cluster_reorder)
  , Jaccard_Mean = rowMeans(jcrdM)
)
# Remove cluster 16
mean_jcrd_DF <- mean_jcrd_DF[1:16, ]
jcrdM <- jcrdM[1:16, ]

# Output as csv
write.csv(mean_jcrd_DF, paste0(outTable, "JaccardMean.csv"), quote = FALSE
  , row.names = FALSE)
# Jaccard index for each boostrap
write.csv(jcrdM, paste0(outTable, "JaccardBoot.csv"), quote = FALSE
  , row.names = TRUE)

# Barplot of jaccard index
ggplot(mean_jcrd_DF, aes(y = Jaccard_Mean, x = Cluster_Reorder
  , fill = Cluster_Reorder)) +
  geom_bar(stat = "identity") +
  ylab("Jaccard index") +
  ggtitle(paste0(graphCodeTitle, "\n\nSeurat cluster stability: Jaccard Index"))
ggsave(paste0(outGraph, "Jaccard_Index.pdf"))
# Paper
ggplot(mean_jcrd_DF, aes(y = Jaccard_Mean, x = Cluster
  , fill = Cluster_Reorder)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  # ylim(c(0,1)) +
  ylab("Jaccard index") +
  ggplot_set_theme_publication +
  ggtitle(paste0(graphCodeTitle, "\n\nSeurat cluster stability: Jaccard Index"))
ggsave(paste0(outGraph, "Jaccard_Index_paper.pdf"))
################################################################################
