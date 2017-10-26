
# Damon Polioudakis
# 2016-11-15
# Determine number of marker genes output by Seurat
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)

## Input marker lists
inMrk <- c("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample50p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample25p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample10p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_DsCl50p_Marker_Genes_Clusters_Vs_All.txt"
  , "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_DsCl10p_Marker_Genes_Clusters_Vs_All.txt")

# Load marker lists derived from different alignments as list of dataframes
mrkLDF <- lapply(inMrk, function(inMrk) read.table(inMrk, header = TRUE))
names(mrkLDF) <- c("Exon", "Exon ds50p", "Exon ds25p", "Exon ds10p"
  , "Exon dsCell50p", "Exon dsCell10p")

# Make output graphs directory
dir.create("../analysis/graphs", recursive = TRUE)

## Variables
graphCodeTitle <- "Seurat_Number_Of_Markers.R"
outGraphPfx <- "../analysis/graphs/Seurat_Number_Of_Markers_"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 10))
################################################################################

## Total number of marker genes for all clusters for each alignment
df <- data.frame(lapply(mrkLDF, function(mrkDF) {length(mrkDF$gene)}))

# Format for ggplot
ggDF <- data.frame(NUMBER = t(df))
ggDF$ALIGNMENT <- row.names(ggDF)
ggDF$ALIGNMENT <- c("Exon", "Exon downsample 50%", "Exon downsample 25%", "Exon downsample 10%"
  , "Exon downsample cells 50%", "Exon downsample cells 10%")

# Bar plot
ggplot(ggDF, aes(x = ALIGNMENT, y = NUMBER, fill = ALIGNMENT)) +
  geom_bar(stat = "identity") +
  geom_text(data = ggDF
    , aes(label = paste(NUMBER), y = (NUMBER + 100))
    , col = "black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.2,1,0.2,0.2), "in")) +
  xlab("Alignment type") +
  ylab("Number of marker genes") +
  ggtitle(paste0(graphCodeTitle
    , "\nTotal number of marker genes for all clusters output by Seurat"
    , "\n"))
ggsave(paste0(outGraphPfx, "Total.pdf"), width = 6, height = 6)

## Number of marker genes per cluster for each alignment
df <- data.frame(lapply(mrkLDF, function(mrkDF) {
  length(mrkDF$gene) / length(unique(mrkDF$cluster))}))
df <- round(df, 1)

# Format for ggplot
ggDF <- data.frame(NUMBER = t(df))
ggDF$ALIGNMENT <- row.names(ggDF)
ggDF$ALIGNMENT <- c("Exon", "Exon downsample 50%", "Exon downsample 25%", "Exon downsample 10%"
  , "Exon downsample cells 50%", "Exon downsample cells 10%")

# Bar plot
ggplot(ggDF, aes(x = ALIGNMENT, y = NUMBER, fill = ALIGNMENT)) +
  geom_bar(stat = "identity") +
  geom_text(data = ggDF
    , aes(label = paste(NUMBER), y = (NUMBER + 10))
    , col = "black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.2,1,0.2,0.2), "in")) +
  xlab("Alignment type") +
  ylab("Number of marker genes") +
  ggtitle(paste0(graphCodeTitle
    , "\nNumber of marker genes per cluster output by Seurat"
    , "\n"))
ggsave(paste0(outGraphPfx, "Per_Cluster.pdf"), width = 6, height = 6)
