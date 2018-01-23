# Damon Polioudakis
# 2018-01-16
# WGCNA test soft power and dendrogram cutting parameters
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
require(biomaRt)
source("Function_Library.R")

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# disableWGCNAThreads()

# WGCNA output
# load("../analysis/WGCNA_Workspace_VarGenes.RData")
# vgAdjTomLLL <- adjTomLLL
# vgDendLLLL <- dendLLLL
# load("../analysis/WGCNA_Workspace_AllGenes.RData")

# Expression and metadata data stored as Seurat object
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
noCentExM <- ssNoCentExM
centSO <- ssCentSO

# Variables
outGraph <- "../analysis/graphs/WGCNA/DS2-11/Test_Parameters/WGCNA_Test_Parameters_"
outAnalysis <- "../analysis/analyzed_data/WGCNA/DS2-11/Test_Parameters/WGCNA_Test_Parameters_"
graphCodeTitle <- "WGCNA_Test_Parameters.R DS2-11"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outAnalysis), recursive = TRUE)

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

# Merge modules based on ME function
# Args: Modules colors, Cut height to merge ME
Merge_Modules_ME <- function(exM, genesModuleColor, cutHeightMergeME) {
  # Call an automatic merging function
  # merged: The merged module colors
  # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
  merged <- mergeCloseModules(exprData = exM, colors = genesModuleColor,
    cutHeight = cutHeightMergeME, trapErrors = TRUE)
  print("Table of genes per module:")
  print(table(merged$colors))
  labels2colors(merged$colors)
}
################################################################################

# ### Choosing the soft-thresholding power: analysis of network topology
#
# # Choose a set of soft-thresholding powers
# powers <- c(seq(2, 30, 2), seq(35, 100, 5))
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

WGCNA_Adjacency_TOM <- function(exM, softPower = 10) {

  print(paste0("Adjacency and TOM for power: ", softPower))

  print("Starting adjacency calculation...")
  adjacency <- adjacency(t(exM), power = softPower, corFnc = "bicor"
    , type = "signed")
  print("Finished adjacency calculation...")

  print("Starting TOM calculation...")
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM <- 1-TOM
  print("Finished TOM calculation...")

  return(list(adjacency = adjacency, TOM = TOM, dissTOM = dissTOM))
}
################################################################################

### Make Modules

WGCNA_Make_Modules <- function(adj_TOM_disTOM_L, exM) {
  print("WGCNA_Make_Modules")
  # Make tree
  geneTree <- hclust(as.dist(adj_TOM_disTOM_L$dissTOM), method = "average")
  print(str(geneTree))

  # Test different parameters for constructing and merging modules
  # Define arguments to test for cutreeHybrid
  minModSizes <- c(30, 50, 100)
  # deepSplits <- c(2, 4)
  cutHeightMergeMEs <- c(0, 0.1, 0.2, 0.25)
  # minModSizes <- c(30, 50)
  # deepSplits <- c(2)
  # cutHeightMergeMEs <- c(0)

  # Cut tree and merge modules looping through lists of parameters
  modulesColors <- NULL
  moduleParameterLabels <- NULL
  # Module sizes loop
  for (minModSize in minModSizes) {
    # cutreeHybrid to make modules
    geneTreeCut <- cutreeHybrid(dendro = geneTree, pamRespectsDendro = TRUE
      , minClusterSize = minModSize
      , distM = as.matrix(adj_TOM_disTOM_L$dissTOM))
    # Merge modules based off ME
    for (cutHeightMergeME in cutHeightMergeMEs) {
      # Test ME merge cut heights
      modulesColors <- cbind(modulesColors,
        Merge_Modules_ME(t(exM), geneTreeCut$labels, cutHeightMergeME))
      # Make label from parameters used to make each module
      moduleParameterLabels <- c(moduleParameterLabels, paste(
        "MMS=", minModSize
        # , " \nDS=", deepSplit
        , " \nMEcor=", cutHeightMergeME
      ))
    }
  }

  print("Done making modules...")

  return(list(geneTree = geneTree
    , modulesColors = modulesColors
    , moduleParameterLabels = moduleParameterLabels))
}

WGCNA_Plot_Dendrogram_Modules <- function(dendro_modules, name) {
  print("WGCNA_Plot_Dendrogram_Modules")
# Plot dendrograms and modules
  pdf(file = paste0(outGraph, "Dendrogram_Modules_", name, ".pdf")
    , height = 15, width = 20)
  # Plot
  plotDendroAndColors(
    dendro = dendro_modules$geneTree
    , colors = dendro_modules$modulesColors
    , groupLabels = dendro_modules$moduleParameterLabels
    , addGuide = TRUE
    , dendroLabels = FALSE
    , main = paste0("Dendrogram With Different Module Cutting Parameters"
      , "\n", name)
  )
  dev.off()
}
################################################################################

WGCNA_Run_Parameters_Test <- function(
  exM
  , dendroName
  , softPowers = c(SP10 = 10, SP15 = 15)) {

  print("WGCNA_Run_Parameters_Test")

  dendro_modules_L <- lapply(softPowers, function(softPower){
    adj_TOM_disTOM_L <- WGCNA_Adjacency_TOM(exM, softPower = softPower)
    dendro_modules <- WGCNA_Make_Modules(adj_TOM_disTOM_L, exM)
    WGCNA_Plot_Dendrogram_Modules(
      dendro_modules
      , name = paste0(dendroName, "_SoftPower", softPower)
    )
    return(dendro_modules)
  })
  print("Done running WGCNA with different parameters")
  return(dendro_modules_L)
}
################################################################################

### Correlate module eigengenes and Seurat cluster membership

WGCNA_Correlate_ME_Seurat_Cluster_Membership <- function(
  exM
  , modulesColors
  , graphTitle) {

  print("WGCNA_Correlate_ME_Seurat_Cluster_Membership")

  meM <- moduleEigengenes(t(exM), modulesColors)$eigengenes

  ldf <- lapply(names(meM), function(name) {
    datME <- meM[[name]]
    # Correlation of MEs to clusters
    # Cluster IDs
    cls <- centSO@ident
    ldf <- lapply(unique(cls), function(cl) {
      v1 <- rep(0, length(cls))
      v1[cls == cl] <- 1
      # table(v1)
      df <- data.frame(cor(datME, v1))
      colnames(df) <- cl
      return(df)
    })
    # Format for ggplot
    ggDF <- t(round(do.call(cbind, ldf), 2))
    ggDF <- melt(ggDF)
    ggDF$name <- name
    return(ggDF)
  })
  ggDF <- do.call("rbind", ldf)
  ggDF$Var1 <- as.factor(ggDF$Var1)
  gg <- ggplot(ggDF, aes(x = Var1, y = name, fill = value)) +
    # facet_wrap(~name, ncol = 3, scales = "free_x") +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b"
      , space = "Lab", name = "Correlation", limits = c(-1, 1)) +
    geom_text(data = ggDF, aes(x = Var1, y = name, label = value)
      , color = "black", size = 2.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(axis.text.y = element_text(size = 8)) +
    xlab("Seurat Cluster") +
    ylab("ME") +
    ggtitle(graphTitle)
  return(gg)
}

WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper <- function(exM = exM
  , dendro_modules_LL = dendro_modules_LL) {
  # Correlate module ME to Seurat cluster membership
  ggL <- lapply(names(dendro_modules_LL), function(softPower) {
    print(str(dendro_modules_LL))
    modulesColors <- dendro_modules_LL[[softPower]]$modulesColors[ ,2]
    WGCNA_Correlate_ME_Seurat_Cluster_Membership(
      exM = exM
      , modulesColors = modulesColors
      , graphTitle = softPower
    )
  })
  pg <- Plot_Grid(ggPlotsL = ggL, ncol = 1, rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\nME and Cluster ID pearson correlation"
      , "\n")
  )
  return(pg)
}
################################################################################

### Module eigengene expression by Seurat cluster

WGCNA_ME_Mean_Seurat_Cluster <- function(
  exM
  , modulesColors
  , graphTitle) {

  print("WGCNA_ME_Mean_Seurat_Cluster")

  meM <- moduleEigengenes(t(exM), modulesColors)$eigengenes
  meM$CellID <- colnames(exM)
  meM$Cluster <- centSO@ident[match(meM$CellID, names(centSO@ident))]

  df1 <- melt(meM)
  df1 <- aggregate(value~Cluster+variable, df1, mean)

  gg <- ggplot(df1, aes(x = Cluster, y = variable, fill = value)) +
    # facet_wrap(~name, ncol = 3, scales = "free_x") +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b"
      , space = "Lab", name = "ME value") +
    geom_text(data = df1, aes(
        x = Cluster, y = variable, label = round(value, 2)
      ), color = "black", size = 2.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(axis.text.y = element_text(size = 8)) +
    xlab("Seurat Cluster") +
    ylab("Mean ME value") +
    ggtitle(graphTitle)
  return(gg)
}

WGCNA_ME_Mean_Seurat_Cluster_Wrapper <- function(exM = exM
  , dendro_modules_LL = dendro_modules_LL) {
  # Correlate module ME to Seurat cluster membership
  ggL <- lapply(names(dendro_modules_LL), function(softPower) {
    print(softPower)
    modulesColors <- dendro_modules_LL[[softPower]]$modulesColors[ ,2]
    WGCNA_ME_Mean_Seurat_Cluster(
      exM = exM
      , modulesColors = modulesColors
      , graphTitle = softPower
    )
  })
  pg <- Plot_Grid(ggPlotsL = ggL, ncol = 1, rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\nMean ME by Seurat cluster"
      , "\n")
  )
  return(pg)
}
################################################################################

### Overlap with Kang modules

WGCNA_Overlap_KangMods <- function(exM, dendro_modules_LL, softPower){

  print("WGCNA_Overlap_KangMods")

  load("../neurogenesis/orig.data/InVivoData/net2_power16_cutHeigt0.15.RData")
  load("../neurogenesis/orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
  KangAnnotRAW = read.csv("../neurogenesis/orig.data/InVivoData/annot.csv", row.names = 1)
  KangAnnotnet2 = KangAnnotRAW[which(
    rownames(KangAnnotRAW) %in% colnames(datExpr)), ]

  kangMods <- KangAnnotnet2
  kangMods$Modules <- net2$colors

  ## First, we get all the genes from Gencode v19
  getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
  geneAnno1 <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)
  geneAnno1 <- geneAnno1[geneAnno1[,"gene_biotype"]=="protein_coding",]; geneAnno1 <- geneAnno1[!duplicated(geneAnno1[,"hgnc_symbol"]) & geneAnno1[,"hgnc_symbol"]!="",] ## 19163 with hgnc symbols after removing duplicates. We use gene symbols here because the RDNV data comes in gene symbol format. Otherwise, I prefer to use ENSG IDs.

  ## Construct a matrix to store the gene lists - this starts as a matrix with gene information and list membership information, and new lists are appended. If the new list contains the gene, it gets a 1, if not 0, if the gene is not even in the study (not in the background) we mark an NA
  metaMat <- matrix(NA,ncol=1,nrow=nrow(geneAnno1)) ## We can add columns as necessary
  rownames(metaMat) <- geneAnno1[,"hgnc_symbol"]
  ENSGIDs <- geneAnno1[,"ensembl_gene_id"]
  geneLength <- geneAnno1[,"end_position"]-geneAnno1[,"start_position"]
  names(geneLength) <- geneAnno1[,"hgnc_symbol"] ## For use later in making the covariate matrix
  metaX = geneAnno1[match(rownames(metaMat), geneAnno1$hgnc_symbol),"chromosome_name"]

  # Add Kang modules to metaMat
  for(i in 1:length(unique(kangMods$Modules))){
    gename = kangMods$SYMBOL[kangMods$Modules == unique(kangMods$Modules)[i]]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = paste0("Kang_", unique(kangMods$Modules)[i])
  }

  # Add Drop-seq modules to metaMat
  df1 <- data.frame(Gene = row.names(exM)
    , Module = dendro_modules_LL[[softPower]]$modulesColors[ ,2]
  )
  for(i in 1:length(unique(df1$Module))){
    gename = df1$Gene[df1$Module == unique(df1$Module)[i]]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = paste0("DS_", unique(df1$Module)[i])
  }

  metaMat[metaMat>0] <- 1
  table(as.numeric(data.matrix(metaMat))) ## Confirm all values are 1s and 0s

  enrichOR <- enrichP <- matrix(NA,nrow=ncol(metaMat),ncol=ncol(metaMat))
  colnames(enrichOR) <- rownames(enrichP) <- rownames(enrichOR) <- rownames(enrichP) <- colnames(metaMat)


  for (i in 40:ncol(metaMat)) { ## Query Lists; Loop through columns of the data
    for (j in 2:39) { ## Biology/Interpretation Lists
      # print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
      dat1 <- as.numeric(metaMat[,i])
      dat2 <- as.numeric(metaMat[,j])
      glm.out <- glm(dat1~dat2,family=binomial)
      keepdat <- !is.na(dat1)&!is.na(dat2)
      ## Compute glm that asks how well dat2 predicts dat1 when controlling for thiscovar

      enrichP[i,j] <- summary(glm.out)$coefficients[2,4]
      enrichOR[i,j] <- summary(glm.out)$coefficients[2,1]
    }
  }

  Pmat <- enrichP
  diag(Pmat) <- 0
  Pmat[upper.tri(Pmat)] <- t(Pmat)[upper.tri(Pmat)] ## Fill in the other half, keep a symmetric matrix

  ORmat <- enrichOR
  diag(ORmat) <- Inf
  ORmat[upper.tri(Pmat)] <- t(ORmat)[upper.tri(ORmat)] ## Fill in the other half, keep a symmetric matrix

  keepcols2 = c(2:39);        #Interpretation Lists
  keepcols1 = c(40:ncol(metaMat));  #Target Lists

  subORmat <- ORmat[keepcols1,keepcols2]
  subPmat <- Pmat[keepcols1,keepcols2]

  ## FDR correction
  subPmat <- p.adjust(subPmat, method="BH")   ##Correcting for all columns and rows

  dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
  rmval <- subPmat > 0.05 ## Display only values with p > 0.05
  dispMat[rmval] <- 0
  rownames(dispMat) <- rownames(subORmat)

  FDRmat = subPmat; FDRfilt <- FDRmat<0.05

  # ggplot
  df1 <- melt(dispMat)
  df1$FDR <- melt(FDRmat)[ ,1]
  ggplot(df1, aes(x = Var2, y = Var1, fill = -1*log(FDR, 10))) +
    # facet_wrap(~name, ncol = 3, scales = "free_x") +
    geom_tile() +
    scale_fill_distiller(type = "seq", palette = "Blues", direction = 1, name = "-log10 FDR") +
    geom_text(data = df1, aes(
        x = Var2, y = Var1, label = round(value, 1)
      ), color = "black", size = 2.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # theme(axis.text.y = element_text(size = 8)) +
    xlab("Kang") +
    ylab("DS2-11") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nEnrichment of Kang modules"
      , "log2 odds ratio"))
}
################################################################################

# Execute script

Run_Top_1000_Genes <- function() {
  print("Run_Top_1000_Genes")
  # Mean expression of each gene to use for filtering
  mnEx <- sort(rowMeans(noCentExM), decreasing = TRUE)
  # Filter expression matrix to genes of interest
  exM <- as.data.frame(as.matrix(centSO@scale.data))[
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:1000]), ]

  # WGCNA test parameters
  dendro_modules_LL <- WGCNA_Run_Parameters_Test(
    exM = exM
    , dendroName = "Top_1000_genes"
    , softPowers = c(SP10 = 10, SP15 = 15)
  )

  # Correlate module ME to Seurat cluster membership
  WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeClusterCorr_Heatmap_Top_1000_genes.pdf")
    , height = 9, width = 9)

  WGCNA_ME_Mean_Seurat_Cluster_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeMeanByCluster_Heatmap_Top_1000_genes.pdf")
    , height = 9, width = 9)

  WGCNA_Overlap_KangMods(
    exM = exM
    , dendro_modules_LL = dendro_modules_LL
    , softPower = "SP15")
  ggsave(paste0(outGraph, "Enrichment_KangModules_Top_1000_genes.pdf")
    , height = 9, width = 12)
}

Run_Top_5000_Genes <- function() {
  print("Run_Top_5000_Genes")
  # Mean expression of each gene to use for filtering
  mnEx <- sort(rowMeans(noCentExM), decreasing = TRUE)
  # Filter expression matrix to genes of interest
  exM <- as.data.frame(as.matrix(centSO@scale.data))[
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:5000]), ]

  # WGCNA test parameters
  dendro_modules_LL <- WGCNA_Run_Parameters_Test(
    exM = exM
    , dendroName = "Top_5000_genes"
    , softPowers = c(SP10 = 10, SP15 = 15)
  )

  # Correlate module ME to Seurat cluster membership
  WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeClusterCorr_Heatmap_Top_5000_genes.pdf")
    , height = 9, width = 9)

  WGCNA_Overlap_KangMods(
    exM = exM
    , dendro_modules_LL = dendro_modules_LL
    , softPower = "SP15")
  ggsave(paste0(outGraph, "Enrichment_KangModules_Top_5000_genes.pdf")
    , height = 9, width = 12)
}

Run_Top_10000_Genes <- function() {
  print("Run_Top_10000_Genes")
  # Mean expression of each gene to use for filtering
  mnEx <- sort(rowMeans(noCentExM), decreasing = TRUE)
  # Filter expression matrix to genes of interest
  exM <- as.data.frame(as.matrix(centSO@scale.data))[
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:10000]), ]

  # WGCNA test parameters
  dendro_modules_LL <- WGCNA_Run_Parameters_Test(
    exM = exM
    , dendroName = "Top_10000_genes"
    , softPowers = c(SP10 = 10, SP15 = 15)
  )

  # Correlate module ME to Seurat cluster membership
  WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeClusterCorr_Heatmap_Top_10000_genes.pdf")
    , height = 9, width = 9)

  WGCNA_Overlap_KangMods(
    exM = exM
    , dendro_modules_LL = dendro_modules_LL
    , softPower = "SP15")
  ggsave(paste0(outGraph, "Enrichment_KangModules_Top_10000_genes.pdf")
    , height = 9, width = 12)
}

Run_Seurat_Variable_Genes <- function() {
  print("Run_Seurat_Variable_Genes")

  # Filter expression matrix to genes of interest
  exM <- as.data.frame(as.matrix(centSO@scale.data))[
      dimnames(centSO@scale.data)[[1]] %in% centSO@var.genes, ]

  # WGCNA test parameters
  dendro_modules_LL <- WGCNA_Run_Parameters_Test(
    exM = exM
    , dendroName = "Seurat_Variable_genes"
    , softPowers = c(SP10 = 10, SP15 = 15)
  )

  # Correlate module ME to Seurat cluster membership
  WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeClusterCorr_Heatmap_Seurat_Variable_genes.pdf")
    , height = 9, width = 9)

  WGCNA_Overlap_KangMods(
    exM = exM
    , dendro_modules_LL = dendro_modules_LL
    , softPower = "SP15")
  ggsave(paste0(outGraph, "Enrichment_KangModules_Seurat_Variable_genes.pdf")
    , height = 9, width = 12)
}

Run_Top_20000_Genes <- function() {
  print("Run_Top_20000_Genes")
  # Mean expression of each gene to use for filtering
  mnEx <- sort(rowMeans(noCentExM), decreasing = TRUE)
  # Filter expression matrix to genes of interest
  exM <- as.data.frame(as.matrix(centSO@scale.data))[
      dimnames(centSO@scale.data)[[1]] %in% names(mnEx[1:20000]), ]

  # WGCNA test parameters
  dendro_modules_LL <- WGCNA_Run_Parameters_Test(
    exM = exM
    , dendroName = "Top_20000_genes"
    , softPowers = c(SP10 = 10, SP15 = 15)
  )

  # Correlate module ME to Seurat cluster membership
  WGCNA_Correlate_ME_Seurat_Cluster_Membership_Wrapper(exM = exM,
    dendro_modules_LL = dendro_modules_LL)
  ggsave(paste0(outGraph, "MeClusterCorr_Heatmap_Top_20000_genes.pdf")
    , height = 9, width = 9)

  WGCNA_Overlap_KangMods(
    exM = exM
    , dendro_modules_LL = dendro_modules_LL
    , softPower = "SP15")
  ggsave(paste0(outGraph, "Enrichment_KangModules_Top_20000_genes.pdf")
    , height = 9, width = 12)
}


Main_Function <- function () {
  Run_Top_1000_Genes()
  Run_Top_5000_Genes()
  Run_Top_10000_Genes()
  Run_Seurat_Variable_Genes()
  Run_Top_20000_Genes()
}

Main_Function()
################################################################################


# # Check coding cluster membership as a 0 1 vector has same correlation as factor in linear model
# # Pearson correlation should = R
# v1 <- rep(0, length(cl))
# v1[cl == "oRG"] <- 1
# cor(datME, v1)
# mod <- lm(v1~datME$MEblack)
# v2 <- rep("No", length(cl))
# v2[cl == "oRG"] <- "oRG"
# v2 <- as.factor(v2)
# mod <- lm(datME$MEblack~v2)
# summary(mod)
