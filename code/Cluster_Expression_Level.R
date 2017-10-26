# Damon Polioudakis
# 2017-05-09
# Pooled expression level for each gene by cluster
################################################################################

rm(list = ls())

require(Seurat)
require(reshape2)
require(dplyr)

# Drop-seq
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 11)))
theme_update(plot.title = element_text(size = 11))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)

pdf("../analysis/graphs/Cluster_Expression_Level.pdf", width = 12, height = 9)

ldf <- split(as.data.frame(t(fetb@raw.data)), fetb@ident)
ll <- lapply(ldf, function(df) {rowSums(t(df))})
ll <- lapply(ll, function(l) sort(l, decreasing = TRUE))
df <- data.frame(ll)
df$RANK <- c(1:nrow(df))
ggDF <- melt(df, id.vars = "RANK")

## Assign annotated cluster names to clusters
current.cluster.ids <- c("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9")
new.cluster.ids <- c(
  "Excitatory Upper Layer Neuron 1"
  , "Excitatory Neuron"
  , "Excitatory Upper Layer Neuron 2"
  , "Excitatory Deep Layer Neuron"
  , "Intermediate Progenitors"
  , "Interneuron"
  , "Mitotic Progenitors"
  , "oRG"
  , "Oligodendrocyte Precursor"
  , "Endothelial")
ggDF$variable <- as.character(ggDF$variable)
ggDF$variable <- plyr::mapvalues(ggDF$variable, from = current.cluster.ids
  , to = new.cluster.ids)

ggDF$variable <- factor(ggDF$variable, levels = new.cluster.ids)

ggplot(ggDF, aes(y = value, x = RANK, group = 1)) +
  geom_line() +
  facet_wrap("variable", scales = "free_x") +
  coord_cartesian(ylim = c(0, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Raw counts") +
  xlab("Genes ordered by expression") +
  ggtitle(paste0(
    "Cluster_Expression_Level.R"
    , "\n"
    , "\nSummed expression level for each gene for each cluster"
    , "\nRaw counts"
    , "\nY-axis limit set to 100"
    , "\n"
  ))

ldf <- split(as.data.frame(t(fetb@scale.data)), fetb@ident)
ll <- lapply(ldf, function(df) {rowSums(t(df))})
ll <- lapply(ll, function(l) sort(l, decreasing = TRUE))
df <- data.frame(ll)
df$RANK <- c(1:nrow(df))
ggDF <- melt(df, id.vars = "RANK")

ggplot(ggDF, aes(y = value, x = RANK, group = 1)) +
  geom_line() +
  facet_wrap("variable", scales = "free_x") +
  coord_cartesian(ylim = c(0, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Mean centered scaled log2 normalized expression") +
  xlab("Genes ordered by expression") +
  ggtitle(paste0(
    "Cluster_Expression_Level.R"
    , "\n"
    , "\nSummed expression level for each gene for each cluster"
    , "\nMean centered scaled log2 normalized expression"
    , "\nY-axis limit set to 100"
    , "\n"
  ))

dev.off()
