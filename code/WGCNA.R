# Damon Polioudakis
# 2017-05-16
# WGCNA
################################################################################

rm(list = ls())
sessionInfo()
set.seed(27)

require(Seurat)
require(Matrix)
require(WGCNA)
require(cluster)
require(flashClust)
require(reshape2)
require(cowplot)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# disableWGCNAThreads() 

# WGCNA output
load("../analysis/WGCNA_Workspace_VarGenes.RData")
vgAdjTomLLL <- adjTomLLL
vgDendLLLL <- dendLLLL
load("../analysis/WGCNA_Workspace_AllGenes.RData")

# Expression and metadata data stored as Seurat object
load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")

# Variables
outGraph <- "../analysis/graphs/WGCNA_VarGenes_"
outAnalysis <- "../analysis/WGCNA_VarGenes_"
graphCodeTitle <- "WGCNA.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

Soft_Threshold <- function(exDF, powers) {
  exDF <- t(exDF)
  sft = pickSoftThreshold(exDF, powerVector =  powers, verbose =  5
    , corFnc = "bicor")
}

SoftThreshold_Plot <- function(sft, plot.title) {
  # Plot the results:
  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3])*sft$fitIndices[ ,2]
    , xlab = "Soft Threshold (power)"
    , ylab = "Scale Free Topology module Fit,signed R^2",type = "n"
    , main = paste(plot.title))
  text(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3])*sft$fitIndices[ ,2],
    labels = powers,cex = cex1,col = "red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90,col = "red")
  abline(h = 0.80,col = "blue")
  abline(h = 0.70,col = "orange")
  abline(h = 0.60,col = "green")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], sft$fitIndices[ ,5]
    , xlab = "Soft Threshold (power)"
    , ylab = "Mean Connectivity"
    , type = "n"
    , main = paste("Mean connectivity"))
  text(sft$fitIndices[ ,1], sft$fitIndices[ ,5], labels = powers, cex = cex1
    , col = "red")
}

# Module identification using hybrid tree cut:
# Function to construct modules
# Args: minimum module size, deep split
Make_Modules <- function (minModSize, deepSplit, geneTree, dissTOM) {
  print("Treecut arguments:")
  # print(c("minModSize"=minModSize,"cut.height"=cutHeightMergeME, "deepSplit"=ds))
  print(c(minModSize, deepSplit))
  tree = cutreeHybrid(dendro = geneTree, pamRespectsDendro = TRUE
    , minClusterSize = minModSize
    # , cutHeight = 0.967
    , deepSplit = deepSplit, distM = as.matrix(dissTOM))
  print("Table of genes per module:")
  print(table(tree$labels))
  tree$labels
}

# Merge modules based on ME function
# Args: Modules colors, Cut height to merge ME
Merge_Modules_ME <- function (exDF, genesModuleColor, cutHeightMergeME) {
  # Call an automatic merging function
  # merged: The merged module colors
  # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
  merged <- mergeCloseModules(exprData = exDF, colors = genesModuleColor,
    cutHeight = cutHeightMergeME, trapErrors = TRUE)
  print("Table of genes per module:")
  print(table(merged$colors))
  labels2colors(merged$colors)
}
################################################################################

### Filter expression matrix

# Mean expression of each gene to use for filtering
mnEx <- sort(rowMeans(seuratO@scale.data), decreasing = TRUE)
# Check top 10% of genes sorted by mean expression
names(mnEx[1:(length(mnEx)*0.1)])

# List of expression matrices after different filters
exLDF <- list(
  
  # Regression normalization
  # Variable genes + 1000 genes detected threshold
  VarGenes_1000genesDetected = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes
    , seuratO@data.info$nGene > 1000]
  # Variable genes + top 10% expressed genes
  # 479
  , VarGenes_Top10 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes &
      dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # Variable genes + top 20% expressed genes
  # 1140
  , VarGenes_Top20 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes &
      dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  # Variable genes + top 30% expressed genes
  #
  , VarGenes_Top30 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes &
      dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.3)]), ]
  # Variable genes + top 40% expressed genes
  #
  , VarGenes_Top40 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes &
      dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.4)]), ]
  # Variable genes + top 50% expressed genes
  # 
  , VarGenes_Top50 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% seuratO@var.genes &
      dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.5)]), ]
  # Top 10% expressed genes
  # 2320
  , Top10 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # Top 20% expressed genes
  # 4641
  , Top20 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  # Top 30% expressed genes
  #
  , Top30 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.3)]), ]
  # Top 40% expressed genes
  #
  , Top40 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.4)]), ]
  # Top 50% expressed genes
  #
  , Top50 = as.data.frame(as.matrix(seuratO@scale.data))[
    dimnames(seuratO@scale.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.5)]), ]
  
  # # Scale factor normalization
  # # Variable genes + 1000 genes detected threshold
  # , VarGenes_1000genesDetected_ScFa = as.data.frame(as.matrix(seuratO@data))[
  #   dimnames(seuratO@data)[[1]] %in% seuratO@var.genes
  #   , seuratO@data.info$nGene > 1000]
  # # Variable genes + top 10% expressed genes
  # # 479
  # , VarGenes_Top10_ScFa = as.data.frame(as.matrix(seuratO@data))[
  #   dimnames(seuratO@data)[[1]] %in% seuratO@var.genes &
  #     dimnames(seuratO@data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # # Variable genes + top 20% expressed genes
  # # 1140
  # , VarGenes_Top20_ScFa = as.data.frame(as.matrix(seuratO@data))[
  #   dimnames(seuratO@data)[[1]] %in% seuratO@var.genes &
  #     dimnames(seuratO@data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  # # Top 10% expressed genes
  # # 2320
  # , Top10_ScFa = as.data.frame(as.matrix(seuratO@data))[
  #   dimnames(seuratO@data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  # # Top 20% expressed genes
  # # 4641
  # , Top20_ScFa = as.data.frame(as.matrix(seuratO@data))[
  #   dimnames(seuratO@data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  
  #   # Raw counts
  #   # Variable genes + 1000 genes detected threshold
  #   , VarGenes_1000genesDetected_Raw = as.data.frame(as.matrix(seuratO@raw.data))[
  #     dimnames(seuratO@raw.data)[[1]] %in% seuratO@var.genes
  #     , seuratO@data.info$nGene > 1000]
  #   # Variable genes + top 10% expressed genes
  #   # 479
  #   , VarGenes_Top10_Raw = as.data.frame(as.matrix(seuratO@raw.data))[
  #     dimnames(seuratO@raw.data)[[1]] %in% seuratO@var.genes &
  #       dimnames(seuratO@raw.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  #   # Variable genes + top 20% expressed genes
  #   # 1140
  #   , VarGenes_Top20_Raw = as.data.frame(as.matrix(seuratO@raw.data))[
  #     dimnames(seuratO@raw.data)[[1]] %in% seuratO@var.genes &
  #       dimnames(seuratO@raw.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
  #   # Top 10% expressed genes
  #   # 2320
  #   , Top10_Raw = as.data.frame(as.matrix(seuratO@raw.data))[
  #     dimnames(seuratO@raw.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.1)]), ]
  #   # Top 20% expressed genes
  #   # 4641
  #   , Top20_Raw = as.data.frame(as.matrix(seuratO@raw.data))[
  #     dimnames(seuratO@raw.data)[[1]] %in% names(mnEx[1:(length(mnEx)*0.2)]), ]
)

# Subset for testing
# idx <- sample(1:ncol(exNmVgNg1000DF), 100)
# exNmVgNg1000DF <- exNmVgNg1000DF[ ,idx]
################################################################################

# ### Choosing the soft-thresholding power: analysis of network topology
# 
# # Choose a set of soft-thresholding powers
# powers <- c(1:30, seq(35, 100, 5))
# 
# # Call the network topology analysis function
# sftL <- lapply(exLDF, function(df) {
#   Soft_Threshold(df, powers)
# })
# 
# # Plot
# pdf(paste0(outGraph, "Power.pdf"), height = 6, width = 12)
# lapply(names(sftL), function(name){
#   sft <- sftL[[name]]
#   SoftThreshold_Plot(sft
#     , paste0("Scale independence ", name))
# })
# dev.off()
################################################################################

### Adjacency and TOM

# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.

# Soft powers to loop through
softPowers <- c(SP10 = 10, SP20 = 20, SP30 = 30, SP70 = 70)
# softPowers <- c(SP10 = 10, SP20 = 20)

# # Calculate adjacency, TOM, and dissTOM for different filtered expression matrices at different soft powers
# # Loop through expression matrices, then soft powers
# # Output list of list of adjacencies, TOMs, and dissTOMs
# # List 1 is split by filtered expression matrices
# # List 2 is split by soft powers
# # List 3 is adj, TOM, and dissTOM matrices
# adjTomLLL <- lapply(exLDF, function(df) {
#   
#   exDF <- df
#   
#   ldf <- lapply(softPowers, function(sp) {
#     
#     print(paste0("Adjacency and TOM for power: ", sp))
#     
#     print("Starting adjacency calculation...")
#     adjacency <- adjacency(t(df), power = sp, corFnc = "bicor"
#       , type = "signed")
#     print("Finished adjacency calculation...")
#     
#     print("Starting TOM calculation...")
#     TOM <- TOMsimilarity(adjacency, TOMType = "signed")
#     dissTOM <- 1-TOM
#     print("Finished TOM calculation...")
#     
#     return(list(adjacency = adjacency, TOM = TOM, dissTOM = dissTOM))
#   })
#   ldf$exDF <- exDF
#   return(ldf)
# })
# 
# save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))
# # save(adjTomLLL, file = paste0(outAnalysis, "SP", sp, "_Adjacency_TOM.rdat"))
################################################################################

### Make Modules

# Make trees and test different cutting and merging parameters
# Outputs nested list
# List 1 is split by filtered expression matrices
# List 2 is split by soft powers
# List 3 is list of dendrogram, module colors, and module making parameters
# List 4 Module colors list
# List 4 Module making parameters list
# Loop through list of adj and TOMs for different gene and cell sets
dendLLLL <- lapply(adjTomLLL[7:11], function(adjTomLL) {
  
  # Loop through list of adj and TOMs for different soft powers
  ll <- lapply(adjTomLL[! names(adjTomLL) == "exDF"], function(adjTomL) {
    
    # Make tree
    geneTree <- hclust(as.dist(adjTomL$dissTOM), method = "average")
    print(str(geneTree))
    
    # Test different parameters for constructing and merging modules
    # Define arguments to test for cutreeHybrid
    minModSizes <- c(30, 50, 100)
    deepSplits <- c(2, 4)
    cutHeightMergeMEs <- c(0, 0.1, 0.2, 0.25)
    
    # Cut tree and merge modules looping through lists of parameters
    modulesColors <- NULL
    moduleParameterLabels <- NULL
    # Module sizes loop
    for (minModSize in minModSizes) {
      # Deep split loop
      for (deepSplit in deepSplits) {
        # cutreeHybrid to make modules
        module <- Make_Modules(minModSize, deepSplit, geneTree, adjTomL$dissTOM)
        # Merge modules based off ME
        for (cutHeightMergeME in cutHeightMergeMEs) {
          # Test ME merge cut heights
          exDF <- adjTomLL$exDF
          modulesColors <- cbind(modulesColors,
            Merge_Modules_ME(t(exDF), module, cutHeightMergeME))
          # Make label from parameters used to make each module
          moduleParameterLabels <- c(moduleParameterLabels, paste(
            "MMS=", minModSize
            , " \nDS=", deepSplit
            , " \nMEcor=", cutHeightMergeME
          ))
        }
      }
    }
    print("Done making modules...")
    return(list(geneTree = geneTree
      , modulesColors = modulesColors
      , moduleParameterLabels = moduleParameterLabels))
  })
  return(ll)
})

save(list = ls(), file = paste0(outAnalysis, "Workspace_AllGenes.RData"))

# Plot dendrograms and modules
pdf(file = paste0(outGraph, "Dendrogram_Modules_AllGenes.pdf"), height = 25, width = 20)

# Loop list split by filtered expression matrices
lapply(names(dendLLLL), function(nameEx) {
  
  print(nameEx)
  dendLLL <- dendLLLL[[nameEx]]
  
  # Loop list split by soft power
  lapply(names(dendLLL), function(nameSp) {
    
    print(nameSp)
    dendLL <- dendLLL[[nameSp]]  
    
    # Plot
    plotDendroAndColors(
      dendro = dendLL$geneTree
      , colors = dendLL$modulesColors
      , groupLabels = dendLL$moduleParameterLabels
      , addGuide = TRUE
      , dendroLabels = FALSE
      , main = paste0("Dendrogram With Different Module Cutting Parameters"
        , "\n", nameEx, " ", nameSp)
    )
  })
})
dev.off()
################################################################################

# # Cluster
# geneTree <- hclust(as.dist(adjTomLLL[["VarGenes_1000genesDetected"]]$dissTOM)
#   , method = "average")
# print(str(geneTree))
# modulesColors <- NULL
# moduleParameterLabels <- NULL
# # cutreeHybrid
# module <- Make_Modules(minModSize = 100, deepSplit = 4, geneTree
#   , adjTomLLL[["VarGenes_1000genesDetected"]]$dissTOM)
# # ME merge cut heights
# modulesColors <- cbind(modulesColors
#   , Merge_Modules_ME(t(adjTomLLL[["VarGenes_1000genesDetected"]]$exDF)
#     , module, cutHeightMergeME = 0.1))

# ME
meL <- list(
  VariableGenes_Top20_SP10 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top20"]]$exDF)
    , vgDendLLLL[["VarGenes_Top20"]][["SP10"]]$modulesColors[ ,14])$eigengenes
  , VariableGenes_Top20_SP20 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top20"]]$exDF)
    , vgDendLLLL[["VarGenes_Top20"]][["SP20"]]$modulesColors[ ,14])$eigengenes
  , VariableGenes_Top20_SP30 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top20"]]$exDF)
    , vgDendLLLL[["VarGenes_Top20"]][["SP30"]]$modulesColors[ ,14])$eigengenes
  
  , VariableGenes_Top30_SP10 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top30"]]$exDF)
    , vgDendLLLL[["VarGenes_Top30"]][["SP10"]]$modulesColors[ ,14])$eigengenes
  , VariableGenes_Top30_SP20 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top30"]]$exDF)
    , vgDendLLLL[["VarGenes_Top30"]][["SP20"]]$modulesColors[ ,14])$eigengenes
  , VariableGenes_Top30_SP30 = moduleEigengenes(t(vgAdjTomLLL[["VarGenes_Top30"]]$exDF)
    , vgDendLLLL[["VarGenes_Top30"]][["SP20"]]$modulesColors[ ,14])$eigengenes
  
  , AllGenes_Top30_SP10 = moduleEigengenes(t(adjTomLLL[["Top30"]]$exDF)
    , dendLLLL[["Top30"]][["SP10"]]$modulesColors[ ,14])$eigengenes
  , AllGenes_Top30_SP20 = moduleEigengenes(t(adjTomLLL[["Top30"]]$exDF)
    , dendLLLL[["Top30"]][["SP20"]]$modulesColors[ ,14])$eigengenes
  , AllGenes_Top30_SP30 = moduleEigengenes(t(adjTomLLL[["Top30"]]$exDF)
    , dendLLLL[["Top30"]][["SP20"]]$modulesColors[ ,14])$eigengenes
) 
# Correlation of ME to each other
signif(cor(meL[[2]], use = "p"), 2)

# # Filter cluster IDs to cells with >1000 genes detected (used for clustering)
# cls <-seuratO@ident[seuratO@data.info$nGene > 1000]

ldf <- lapply(names(meL), function(name) {
  datME <- meL[[name]]
  # Correlation of MEs to clusters
  # Cluster IDs
  cls <- seuratO@ident
  ldf <- lapply(unique(cls), function(cl) {
    v1 <- rep(0, length(cls))
    v1[cls == cl] <- 1
    # table(v1)
    df <- data.frame(cor(datME, v1))
    colnames(df) <- cl
    return(df)
  })
  print(ldf[[1]])
  # Correlation of ME to combined excitatory 1 and 2
  v1 <- rep(0, length(cls))
  v1[cls %in% c("Excitatory Upper Layer Neuron 1"
    , "Excitatory Upper Layer Neuron 2")] <- 1
  df <- data.frame(cor(datME, v1))
  colnames(df) <- "Excitatory 1 and 2 combined"
  ldf[["Combined"]] <- df
  print(df)
  # Format for ggplot
  ggDF <- t(round(do.call(cbind, ldf), 2))
  ggDF <- melt(ggDF)
  ggDF$name <- name
  return(ggDF)
})
ggDF <- do.call("rbind", ldf)
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  facet_wrap(~name, ncol = 3, scales = "free_x") +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b"
    , space = "Lab", name = "Correlation", limits = c(-1, 1)) +
  geom_text(data = ggDF, aes(x = Var2, y = Var1, label = value)
    , color = "black", size = 2.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # theme(axis.text.y = element_text(size = 8)) +
  ylab("Seurat Cluster") +
  xlab("ME") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nME and Cluster ID pearson correlation"
    , "\n"))
ggsave(paste0(outGraph, "MeClusterCorr_Heatmap.pdf"), height = 10, width = 17)


# Check coding cluster membership as a 0 1 vector has same correlation as factor in linear model
# Pearson correlation should = R
v1 <- rep(0, length(cl))
v1[cl == "oRG"] <- 1
cor(datME, v1)
mod <- lm(v1~datME$MEblack)
v2 <- rep("No", length(cl))
v2[cl == "oRG"] <- "oRG"
v2 <- as.factor(v2)
mod <- lm(datME$MEblack~v2)
summary(mod)
# Heatmap of correlation

# Heatmap of known markers module membership


# Code from WGCNA tutorial
# # Set up variables to contain the module-trait correlations
# moduleTraitCor = list();
# moduleTraitPvalue = list();
# # Calculate the correlations
# for (set in 1:nSets)
# {
#   moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
#   moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
# }