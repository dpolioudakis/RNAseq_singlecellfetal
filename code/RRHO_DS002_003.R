# Damon Polioudakis
# 2017-02-15
# RRHO of Seurat output marker genes
################################################################################

rm(list=ls())

# setwd("~/Desktop/Single_Cell_RNAseq/code/dev")
# load(file = "PollenCellTypesGenes.R")
# mkDF <- read.table("Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_test_ens.txt"
#   , header = TRUE)
# plDF <- read.csv("mmc5.csv", header = TRUE, fill = TRUE)

require(RRHO)
require(reshape2)
require(ggplot2)

options(stringsAsFactors = FALSE)

# Cell types in fetal brain (Pollen et al S4)
plDF <- read.csv("../source/Pollen_S4.csv", header = TRUE, fill = TRUE)

# Seurat marker lists
mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_test.txt"
  , header = TRUE)

## Variables
graphCodeTitle <- "RRHO_DS002_003.R"
outGraphPfx <- "../analysis/graphs/RRHO_DS002_003_"

## Output directory
dir.create(dirname(outGraphPfx), showWarnings = FALSE, recursive = TRUE)
################################################################################

### Functions

# Function to find and change values
# Inputs:
#   [values] Vector to change
#   [from] Value to change from
#   [to] Value to change to, in the same order as [from]
Map_Values <- function (values, from, to) {
  for (i in 1:length(from)) {
    values[values == from[i]] <- to[i]
  }
  return(values)
}
################################################################################

# Add column of cluster annotations to marker gene data frame
# Cluster IDs
currentIDs <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
# Cluster annotations
newIDs <- c(
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
# Add column of cluster annotations
mkDF$CLUSTER_ANNOTATION <- Map_Values(mkDF$cluster, currentIDs, newIDs)

# Split marker gene DF by cluster
mkLDF <- split(mkDF[c("gene", "p_val")], mkDF$CLUSTER_ANNOTATION)
# Order marker genes by pvalue
mkLDF <- lapply(mkLDF, function(df) {df[order(df$p_val), ]})

# Select Pollen gene correlation columns
plDF <- plDF[ ,c(1,6:9)]

## Make list of RRHO output hypermats
# Initialize output list
hypermatL <- list()
# Counter of output list
k <- 0
# Loop through gene lists to compare and run RRHO
# Pollen marker gene lists
for (j in c(2:length(plDF))){
  
  # Marker gene lists
  for (i in c(1:length(mkLDF))){

    # Output list counter
    k <- k + 1
    
    # Name of list
    name_i <- names(mkLDF)[i]
    print(name_i)
    # Pollen
    name_j <- names(plDF)[j]
    print(name_j)
    
    # Select Pollen gene list and order
    df1 <- data.frame(GENE = plDF[ ,1], RANK = plDF[ ,j])
    df1 <- df1[order(-df1$RANK), ]
    df1$RANK <- c(length(plDF[ ,2]):1)
    df1 <- df1[! duplicated(df1$GENE), ]
    print(head(df1))
    
    # Select marker gene list
    df2 <- data.frame(GENE = mkLDF[[i]][ ,1]
      , RANK = c(length(mkLDF[[i]][ ,1]):1))
    df2 <- df2[! duplicated(df2$GENE), ]
    print(head(df2))
    
    # Run RRHO
    tmp <- RRHO(df2, df1, stepsize = 100
      , labels = c(paste0("_DS_", name_i), paste0("Pollen_", name_j))
      , alternative = "enrichment", plots = FALSE)
    print(str(tmp))
    # Add RRHO hypermat to output list
    hypermatL[[k]] <- tmp$hypermat
    # Set output list item name
    names(hypermatL)[k] <- paste0("Pollen_", name_j, " DS_", name_i)
  }
}

### Plot the RRHO

jetColors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F"
  , "yellow", "#FF7F00", "red", "#7F0000"))
colorMap = jetColors(100)

# Format for ggplot
ggLDF <- lapply(hypermatL, melt)
# Add column of cluster IDs for each gene list comparison
for(i in 1:length(ggLDF)) {ggLDF[[i]]$LABEL <- names(ggLDF)[i]}
# List of data frames to single data frame
ggDF <- do.call(rbind, ggLDF)
# Cap -log10 p-value at 100 for consistent scale on graph
ggDF$value[ggDF$value > 100] <- 100
# Format labels
ggDF$LABEL <- gsub("Correlation.", "", ggDF$LABEL)
ggDF$LABEL_POLLEN <- gsub(" DS_.*", "", ggDF$LABEL)
ggDF$LABEL_DS <- gsub("Pollen_.* DS", "DS", ggDF$LABEL)

# ggplot
ggO <- ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  facet_grid(LABEL_POLLEN~LABEL_DS) +
  geom_tile() +
  scale_fill_gradientn(colours = colorMap, limits = c(0, 100)
    , name = "-log P-value") +
  theme(strip.text.x = element_text(size = 8, angle = 90, face = "bold")) +
  theme(strip.text.y = element_text(size = 8, face = "bold")) +
  xlab("Drop-seq") +
  ylab("Pollen") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nDS-002-003 Seurat marker ranking and Pollen S4"
    , "\nRank rank hypergeometric overlap"
    , "\nStepsize = 100"))
ggsave(paste0(outGraphPfx, "Pollen.pdf"), ggO, width = 10, height = 8)
ggsave(paste0(outGraphPfx, "Pollen.png"), ggO, width = 10, height = 8)
################################################################################