# Damon Polioudakis
# 2017-11-01
# Separate cycling vRG and oRG

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(reshape2)
require(ggplot2)
require(cowplot)
# require(pheatmap)
require(RColorBrewer)
require(viridis)
# require(fdrtool)
require(monocle)
require(WGCNA)
require(gridExtra)
require(ggpubr)
source("Function_Library.R")
source("Seurat_Cluster_Cycling_vRGoRG_Functions.R")

options(stringsAsFactors = FALSE)

# Keep CC genes from variable gene list used for clustering
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TESTcluster0278910_seuratO.Robj")
# centSO <- ssCentSO; rm(ssCentSO)
# noCentExM <- ssNoCentExM; rm(ssNoCentExM)

# Cell cycle markers used for phase determination (Tirosh et al. 2016)
ccGenes <- readLines(con = "../source/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
sGenes <- ccGenes[1:43]
g2mGenes <- ccGenes[44:98]

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv", header = TRUE
  , fill = TRUE)

## Variables
graphCodeTitle <- "Seurat_Cluster_Cycling_vRGoRG_DS2-11.R"
outGraph <- "../analysis/graphs/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outTable <- "../analysis/tables/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"
outData <- "../analysis/analyzed_data/Seurat_Cluster_DS2-11_Cycling_vRGoRG/Seurat_Cluster_Cycling_vRGoRG_DS2-11_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outData), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)
################################################################################

### DE RG(cluster 7,9) vs IP(2); RG(7,9) vs Neuron(0); IP(2) vs Neuron(0)

DE_Clusters_Vs_Clusters <- function(
  clusters1, clusters2, positive_DE_Label, negative_DE_label) {

  cell_IDs_1 <- names(centSO@ident)[centSO@ident %in% c(clusters1)]
  cell_IDs_2 <- names(centSO@ident)[centSO@ident %in% c(clusters2)]

  deDF <- DE_By_CellIDs(
    cell_IDs_1 = cell_IDs_1
    , cell_IDs_2 = cell_IDs_2
    , positive_DE_label = positive_DE_Label
    , negative_DE_label = negative_DE_label)

  return(deDF)
}

DE_By_CellIDs <- function(
  cell_IDs_1, cell_IDs_2, positive_DE_label = NA, negative_DE_label = NA) {
  ids <- c(cell_IDs_1, cell_IDs_2)
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% ids]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% ids
    , c("nUMI", "librarylab", "individual", "CELL")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$groups <- "cell_IDs_1"
  termsDF$groups[termsDF$CELL %in% cell_IDs_2] <- "cell_IDs_2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ groups+nUMI+librarylab+individual")
  # Format LM DE
  deDF <- data.frame(
    Log2_FC_Group1_vs_Group2 = deLM$coefmat[ ,"groupscell_IDs_2"]
    , Pvalue = deLM$pvalmat[ ,"groupscell_IDs_2"]
  )
  deDF$Gene = row.names(deDF)
  # Make DE cell_IDs_1 positive fold change
  deDF$Log2_FC_Group1_vs_Group2 <- deDF$Log2_FC_Group1_vs_Group2 * -1
  # Label DE group
  deDF$DE_Group <- NA
  deDF$DE_Group[deDF$Log2_FC_Group1_vs_Group2 > 0] <- positive_DE_label
  deDF$DE_Group[deDF$Log2_FC_Group1_vs_Group2 < 0] <- negative_DE_label
  # Order by fold change
  deDF <- deDF[order(deDF$Log2_FC_Group1_vs_Group2), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method = "BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))
  return(deDF)
}

# DE cell types
# RG vs IP clusters
de_RG_v_IP_DF <- DE_Clusters_Vs_Clusters(
  clusters1 = c(7,9), clusters2 = 2
  , positive_DE_Label = "RG", negative_DE_label = "IP"
)
# RG vs Migrating Neuron clusters
de_RG_v_Ne_DF <- DE_Clusters_Vs_Clusters(
  clusters1 = c(7,9), clusters2 = 0
  , positive_DE_Label = "RG", negative_DE_label = "Neuron"
)
# IP vs Migrating Neuron clusters
de_IP_v_Ne_DF <- DE_Clusters_Vs_Clusters(
  clusters1 = 2, clusters2 = 0
  , positive_DE_Label = "IP", negative_DE_label = "Neuron"
)

# Format and save
de_RG_v_IP_DF$Comparison <- "RG_vs_IP"
de_RG_v_Ne_DF$Comparison <- "RG_vs_Neuron"
de_IP_v_Ne_DF$Comparison <- "IP_vs_Neuron"
deDF <- rbind(de_RG_v_IP_DF, de_RG_v_Ne_DF, de_IP_v_Ne_DF)
# Save as csv
write.csv(deDF, file = paste0(outTable, "DE_CellTypes.csv")
  , quote = FALSE)


## DE transition states

Format_Cell_IDs_1_2 <- function(
  cellIDs1, cellIDs2, comparison_label, positive_DE_label, negative_DE_label){
  de_to_run_DF <- data.frame(
    Cell_IDs = c(cellIDs1, cellIDs2)
    , Group = c(rep("cell_IDs_1", length(cellIDs1))
      , rep("cell_IDs_2", length(cellIDs2)))
    , Comparison = comparison_label
    , Positive_DE_Label = positive_DE_label
    , Negative_DE_Label = negative_DE_label
  )
  return(de_to_run_DF)
}

# RG+ vs RG+IP+
# S phase
Subset_Cell_IDs_RG_RGIP_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(7,9,8)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(8)])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "RG_vs_RGIP_Sphase"
    , positive_DE_label = "RG"
    , negative_DE_label = "RGIP"
  )
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_RG_RGIP_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(7,9,10)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(10)])
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "RG_vs_RGIP_G2Mphase"
    , positive_DE_label = "RG"
    , negative_DE_label = "RGIP"
  )
  return(de_to_run_DF)
}

# RG+ vs RG+Neuron+
# S phase
Subset_Cell_IDs_RG_RGNeuron_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(7,9,8)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.5 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(8)])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "RG_vs_RGNeuron_Sphase"
    , positive_DE_label = "RG"
    , negative_DE_label = "RGNeuron"
  )
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_RG_RGNeuron_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP < 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(7,9,10)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(10)])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "RG_vs_RGNeuron_G2Mphase"
    , positive_DE_label = "RG"
    , negative_DE_label = "RGNeuron"
  )
  return(de_to_run_DF)
}

# IP+ vs IP+Neuron+
# S phase
Subset_Cell_IDs_IP_IPNeuron_Sphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(2,8)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(8)])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "IP_vs_IPNeuron_Sphase"
    , positive_DE_label = "IP"
    , negative_DE_label = "IPNeuron"
  )
  return(de_to_run_DF)
}
# G2/M phase
Subset_Cell_IDs_IP_IPNeuron_G2Mphase <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  # Cell IDs 1
  cellIDs1 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron < 0.5
    ]
  cellIDs1 <- intersect(
    cellIDs1, names(centSO@ident)[centSO@ident %in% c(2,10)])
  # Cell IDs 2
  cellIDs2 <- row.names(mk_exp_DF)[
    mk_exp_DF$RG < 0.5 &
    mk_exp_DF$IP > 0.5 &
    mk_exp_DF$Neuron > 0.5
    ]
  cellIDs2 <- intersect(
    cellIDs2, names(centSO@ident)[centSO@ident %in% c(10)])
  #
  de_to_run_DF <- Format_Cell_IDs_1_2(
    cellIDs1 = cellIDs1
    , cellIDs2 = cellIDs2
    , comparison_label = "IP_vs_IPNeuron_G2Mphase"
    , positive_DE_label = "IP"
    , negative_DE_label = "IPNeuron"
  )
  return(de_to_run_DF)
}

# Run DE for transition states
de_to_run_DF <- rbind(
  Subset_Cell_IDs_RG_RGIP_Sphase()
  , Subset_Cell_IDs_RG_RGIP_G2Mphase()
  , Subset_Cell_IDs_RG_RGNeuron_Sphase()
  , Subset_Cell_IDs_RG_RGNeuron_G2Mphase()
  , Subset_Cell_IDs_IP_IPNeuron_Sphase()
  , Subset_Cell_IDs_IP_IPNeuron_G2Mphase()
)
transition_state_DE_DF <- do.call("rbind"
  , lapply(unique(de_to_run_DF$Comparison)
    , function(comparison){

    subset_de_to_run_DF <- de_to_run_DF[de_to_run_DF$Comparison == comparison, ]

    cellIDs1 <- subset_de_to_run_DF$Cell_IDs[
      subset_de_to_run_DF$Group == "cell_IDs_1"
      ]
    cellIDs2 <- subset_de_to_run_DF$Cell_IDs[
      subset_de_to_run_DF$Group == "cell_IDs_2"
      ]

    deDF <- DE_By_CellIDs(
      cell_IDs_1 = cellIDs1
      , cell_IDs_2 = cellIDs2
      , positive_DE_label = subset_de_to_run_DF$Positive_DE_Label[1]
      , negative_DE_label = subset_de_to_run_DF$Negative_DE_Label[1]
    )
    deDF$Comparison <- comparison

    return(deDF)
  })
)
# Save as csv
write.csv(transition_state_DE_DF, file = paste0(outTable, "DE_CellTransitionStates.csv")
  , quote = FALSE)
################################################################################

### ME of DE genes in RG+ IP+ Neuron+

## Functions

Make_Gene_List_for_moduleEigengenes <- function(
  foldChange, exM, class1_label, class2_label){

  names(foldChange) <- foldChange
  genesDFL <- lapply(foldChange, function(foldChange){

    genesDF <- data.frame(Genes = row.names(exM)
      , Class = rep("Neither", length(row.names(exM))))
    genesDF$Class[genesDF$Genes %in% deDF$Gene[
      deDF$Log2_FC_Group1_vs_Group2 > foldChange & deDF$FDR < 0.05]] <- class2_label
    genesDF$Class[genesDF$Genes %in% deDF$Gene[
      deDF$Log2_FC_Group1_vs_Group2 < -(foldChange) & deDF$FDR < 0.05]] <- class1_label

    return(genesDF)
  })
  return(genesDFL)
}

Calculate_Module_Eigengene <- function(genesDF, exM){
  print("Calculate_Module_Eigengene")
  subset_exM <- exM[genesDF$Class != "Neither", ]
  genesDF <- genesDF[genesDF$Class != "Neither", ]
  meM <- moduleEigengenes(t(subset_exM), genesDF$Class)$eigengenes
  return(meM)
}

Format_moduleEigengenes_Output <- function(meM, exM, genesDF, mk_exp_DF){
  print("Format_moduleEigengenes_Output")
  meM$CellID <- colnames(exM[genesDF$Class != "Neither", ])
  me_markerFlag_DF <- Marker_Expression_Flag(meM, mk_exp_DF)
  me_markerFlag_DF <- me_markerFlag_DF[ ,c(1:6)]
  me_markerFlag_DF <- melt(me_markerFlag_DF)
  me_markerFlag_DF$Number_DE_Genes <- sum(genesDF$Class != "Neither")
  return(me_markerFlag_DF)
}

T_test_eigengene <- function(me_markerFlag_DF){
  print("T_test_eigengene")
  pvals_DF <- sapply(split(me_markerFlag_DF, me_markerFlag_DF$variable), function(df1) {
    combinations <- combn(unique(df1$Cell_Subset), 2)
    pvals <- apply(combinations, 2, function(y){
      results <- t.test(df1$value[df1$Cell_Subset == y[1]]
        , df1$value[df1$Cell_Subset == y[2]])
      return(results$p.value)
    })
    names(pvals) <- paste(combinations[1, ], "vs", combinations[2, ], "\npvalue")
    pvals <- signif(pvals, 2)
    return(pvals)
  })
  pvals_DF <- as.data.frame(pvals_DF)
  return(pvals_DF)
}

Calculate_ME_For_CellTypeDE_Genes <- function(
  cellIDs, deDF, mk_exp_DF, foldChange, class1_label, class2_label) {

  print("Calculate_ME_For_CellTypeDE_Genes")

  exM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

  genesDFL <- Make_Gene_List_for_moduleEigengenes(
    foldChange = foldChange, exM = exM
    , class1_label = class1_label, class2_label = class2_label)

  print(str(genesDFL))

  me_markerFlag_DF <- lapply(names(genesDFL), function(name){

    print(name)

    genesDF <- genesDFL[[name]]
    meM <- Calculate_Module_Eigengene(
      genesDF = genesDF
      , exM = exM
    )
    me_markerFlag_DF <- Format_moduleEigengenes_Output(
      meM = meM
      , exM = exM
      , genesDF = genesDF
      , mk_exp_DF
    )
    me_markerFlag_DF$Fold_Change_Cutoff <- name
    # Format ME; eg MERG to RG eigengene
    me_markerFlag_DF$variable <- paste(gsub("ME", "", me_markerFlag_DF$variable)
      , "eigengene")
    return(me_markerFlag_DF)
  })
  me_markerFlag_DF <- do.call("rbind", me_markerFlag_DF)

  return(me_markerFlag_DF)
}

## Functions to run

ME_CellType_RG_to_IP_8 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  print(cellIDs)
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
    , mk_exp_DF = mk_exp_DF, class1_label = "IP", class2_label = "RG"
    , foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellType_RG_to_IP_10 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
    , mk_exp_DF = mk_exp_DF, class1_label = "IP", class2_label = "RG"
    , foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellType_RG_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,8)])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
    , mk_exp_DF = mk_exp_DF, class1_label = "Neuron", class2_label = "RG"
    , foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellType_RG_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,10)])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
    , mk_exp_DF = mk_exp_DF
    , class1_label = "Neuron", class2_label = "RG", foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}


ME_CellType_IP_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,8)])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "IP_vs_Neuron", ]
    , mk_exp_DF = mk_exp_DF
    , class1_label = "Neuron", class2_label = "IP", foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellType_IP_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == c(0,10)])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeDE_Genes(
    cellIDs = cellIDs, deDF = deDF[deDF$Comparison == "IP_vs_Neuron", ]
    , mk_exp_DF = mk_exp_DF
    , class1_label = "Neuron", class2_label = "IP", foldChange = c(0, 0.25, 0.5)
  )
  return(me_markerFlag_DF)
}


## Run

ME_CellType_L <- list(
  RG_to_IP_8 = ME_CellType_RG_to_IP_8()
  , RG_to_IP_10 = ME_CellType_RG_to_IP_10()
  , RG_to_Neuron_08 = ME_CellType_RG_to_Neuron_08()
  , RG_to_Neuron_010 = ME_CellType_RG_to_Neuron_010()
  , IP_to_Neuron_08 = ME_CellType_IP_to_Neuron_08()
  , IP_to_Neuron_010 = ME_CellType_IP_to_Neuron_010()
)

save(ME_CellType_L, file = paste0(outData, "ME_CellType.rdata"))
################################################################################

### ME of cell type enriched genes in RG+ IP+ Neuron+

Make_Gene_List_for_moduleEigengenes_CellTypeEnriched <- function(
  foldChange
  , exM
  , seurat_clusters
  , class_label
  ){

  print("Make_Gene_List_for_moduleEigengenes_CellTypeEnriched")

  names(foldChange) <- foldChange
  genesDFL <- lapply(foldChange, function(foldChange){

    # Subset genes to top enriched for Seurat clusters of interest
    genes <- cluster_DE_DF$Gene[
      cluster_DE_DF$Cluster %in% seurat_clusters
      & cluster_DE_DF$Log2_Fold_Change > foldChange
      ]

    genesDF <- data.frame(
      Genes = row.names(exM)
      , Class = rep("Neither", length(row.names(exM)))
    )
    genesDF$Class[genesDF$Genes %in% genes] <- class_label

    return(genesDF)
  })
  return(genesDFL)
}

Calculate_ME_For_CellTypeEnriched_Genes <- function(
  cellIDs
  , mk_exp_DF
  , foldChange
  , seurat_clusters
  , class_label
  ) {

  print("Calculate_ME_For_CellTypeEnriched_Genes")

  exM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

  genesDFL <- Make_Gene_List_for_moduleEigengenes_CellTypeEnriched(
    foldChange = foldChange
    , exM = exM
    , seurat_clusters = seurat_clusters
    , class_label = class_label
  )

  print(str(genesDFL))

  me_markerFlag_DF <- lapply(names(genesDFL), function(name){

    print(name)

    genesDF <- genesDFL[[name]]
    print(table(genesDF$Class))
    meM <- Calculate_Module_Eigengene(
      genesDF = genesDF
      , exM = exM
    )
    me_markerFlag_DF <- Format_moduleEigengenes_Output(
      meM = meM
      , exM = exM
      , genesDF = genesDF
      , mk_exp_DF
    )
    me_markerFlag_DF$Fold_Change_Cutoff <- name
    # Format ME; eg MERG to RG eigengene
    me_markerFlag_DF$variable <- paste(gsub("ME", "", me_markerFlag_DF$variable)
      , "eigengene")
    return(me_markerFlag_DF)
  })
  me_markerFlag_DF <- do.call("rbind", me_markerFlag_DF)

  return(me_markerFlag_DF)
}

## Functions to run

ME_CellTypeEnrichedRG_RG_to_IP_8 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = c(7,9)
    , class_label = "RG"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedIP_RG_to_IP_8 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 2
    , class_label = "IP"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedRG_RG_to_IP_10 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = c(7,9)
    , class_label = "RG"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedIP_RG_to_IP_10 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$Neuron < 0.25 | mk_exp_DF$IP > 0.5 & mk_exp_DF$Neuron < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 2
    , class_label = "IP"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedRG_RG_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = c(7,9)
    , class_label = "RG"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedNeuron_RG_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 0
    , class_label = "Neuron"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedRG_RG_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = c(7,9)
    , class_label = "RG"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedNeuron_RG_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$RG > 0.5 & mk_exp_DF$IP < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$IP < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 0
    , class_label = "Neuron"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedIP_IP_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 2
    , class_label = "IP"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedNeuron_IP_to_Neuron_08 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 0
    , class_label = "Neuron"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedIP_IP_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 2
    , class_label = "IP"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnrichedNeuron_IP_to_Neuron_010 <- function(){
  mk_exp_DF <- Average_MarkersExp_Per_Cell(
    exM = noCentExM, seuratO = centSO)
  cellIDs <- row.names(mk_exp_DF)[
    mk_exp_DF$IP > 0.5 & mk_exp_DF$RG < 0.25 | mk_exp_DF$Neuron > 0.5 & mk_exp_DF$RG < 0.25]
  cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
  me_markerFlag_DF <- Calculate_ME_For_CellTypeEnriched_Genes(
    cellIDs = cellIDs
    , mk_exp_DF = mk_exp_DF
    , seurat_clusters = 0
    , class_label = "Neuron"
    , foldChange = c(0.2, 0.4, 0.5)
  )
  return(me_markerFlag_DF)
}

ME_CellTypeEnriched_L <- list(
  RG_to_IP_8 = rbind(
      ME_CellTypeEnrichedRG_RG_to_IP_8()
    , ME_CellTypeEnrichedIP_RG_to_IP_8()
  )
  , RG_to_IP_10 = rbind(
    ME_CellTypeEnrichedRG_RG_to_IP_10()
    , ME_CellTypeEnrichedIP_RG_to_IP_10()
  )
  , RG_to_Neuron_08 = rbind(
    ME_CellTypeEnrichedRG_RG_to_Neuron_08()
    , ME_CellTypeEnrichedNeuron_RG_to_Neuron_08()
  )
  , RG_to_Neuron_010 = rbind(
    ME_CellTypeEnrichedRG_RG_to_Neuron_010()
    , ME_CellTypeEnrichedNeuron_RG_to_Neuron_010()
  )
  , IP_to_Neuron_08 = rbind(
    ME_CellTypeEnrichedIP_IP_to_Neuron_08()
    , ME_CellTypeEnrichedNeuron_IP_to_Neuron_08()
  )
  , IP_to_Neuron_010 = rbind(
    ME_CellTypeEnrichedIP_IP_to_Neuron_010()
    , ME_CellTypeEnrichedNeuron_IP_to_Neuron_010()
  )
)

save(ME_CellTypeEnriched_L, file = paste0(outData, "ME_CellTypeEnriched.rdata"))
################################################################################

### Correlation of cluster mean expression profile and RG IP cells
print("### Correlation of cluster EG and RG IP cells")

df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)

cellSubsetsL <- list(
  IP_G2M = df1$IP > 0.5 & df1$PHASE %in% "G2M"
  , RG_G2M = df1$RG > 0.5 & df1$PHASE %in% "G2M"

  , IP_S = df1$IP > 0.5 & df1$PHASE %in% "S"
  , RG_S = df1$RG > 0.5 & df1$PHASE %in% "S"

  , IP_G1 = df1$IP > 0.5 & df1$PHASE %in% "G1"
  , RG_G1 = df1$RG > 0.5 & df1$PHASE %in% "G1"

  , IP = df1$IP > 0.5
  , RG = df1$RG > 0.5

  , RG_IPn_G2M = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "G2M"
  , RG_IP_G2M = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G2M"
  , RGn_IP_G2M = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G2M"

  , RG_IPn_S = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "S"
  , RG_IP_S = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "S"
  , RGn_IP_S = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "S"

  , RG_IPn_G1 = df1$RG > 0.5 & df1$IP < 0.5 & df1$PHASE %in% "G1"
  , RG_IP_G1 = df1$RG > 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G1"
  , RGn_IP_G1 = df1$RG < 0.5 & df1$IP > 0.5 & df1$PHASE %in% "G1"

  , RG_IP = df1$RG > 0.5 & df1$IP > 0.5
  , RG_IPn = df1$RG > 0.5 & df1$IP < 0.5
  , RGn_IP = df1$RG < 0.5 & df1$IP > 0.5
)
egL <- lapply(names(cellSubsetsL), function(name) {
  cellSubset <- cellSubsetsL[[name]]
  exM <- noCentExM
  mns <- rowMeans(exM)
  genes <- names(sort(mns, decreasing = TRUE))[1:2500]
  exM <- exM[row.names(exM) %in% genes, cellSubset]
  rowMeans(exM)
  # print(str(exM))
  # colors <- rep("TRUE", nrow(exM))
  # eg <- moduleEigengenes(t(exM), colors)$eigengenes
  # eg <- eg["METRUE"]
  # colnames(eg) <- name
  # return(eg)
})
names(egL) <- names(cellSubsetsL)
corM <- round(cor(data.frame(egL)), 2)

write.csv(corM
  , file = paste0(outTable, "RG_IP_subsets_Mean_correlation_matrix.csv")
  , quote = FALSE)

# randomIDs <- sample(colnames(noCentExM), size = length(cellIDs), replace = FALSE)
# exM <- noCentExM[ ,colnames(noCentExM) %in% randomIDs]
# mnEx <- rowMeans(exM)
# cor(mnEx, eg$MEIPC)
################################################################################

### vRG and oRG gene expression by cluster
print("### vRG and oRG gene expression by cluster")

ggL <- lapply(c("vRG", "oRG", "RG"), function(grouping) {
  genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% grouping]
  genes <- genes[! duplicated(genes)]
  gg <- Gene_Expression_Facet_By_Cluster_ViolinPlot(
    genes = genes
    , exprM = noCentExM[ ,colnames(noCentExM) %in% names(centSO@ident)[centSO@ident %in% c(7,8,9,10)]]
    , clusterIDs = centSO@ident[centSO@ident %in% c(7,8,9,10)]
    , geneOrder = genes
    , ncol = 4
    , ggtitle = grouping)
  gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(gg)
})
Plot_Grid(ggPlotsL = ggL
  , ncol = 1
  , title = paste0(graphCodeTitle
    , "\n\nvRG, oRG, and RG markers expression by cluster")
  , rel_height = 0.1)
ggsave(paste0(outGraph, "vRG_oRG_RG_Expression_ViolinPlot.png")
  , width = 13, height = 11)
################################################################################

### tSNE of subsets of cells

tSNE_CCgenes_Markers <- function(cellIDs, title, df1) {

  # Subset to cells of interest
  ssCentSO <- SubsetData(centSO, cells.use = cellIDs)
  ssCentSO@raw.data <- ssCentSO@raw.data[
    ,colnames(ssCentSO@raw.data) %in% cellIDs]
  ssNoCentExM <- noCentExM[ ,colnames(noCentExM) %in% cellIDs]

  # ssCentSO@data <- ssNoCentExM
  # ssCentSO <- ScaleData(ssCentSO)

  # tSNE
  ssCentSO <- RunTSNE(ssCentSO, dims.use = 1:8, do.fast = TRUE, perplexity = 20)

  # Format Macosko cell cycle genes table
  genesGroupDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(genesGroupDF) <- c("Grouping", "Gene.Symbol")
  genesGroupDF$Gene.Symbol <- gsub(" *", "", genesGroupDF$Gene.Symbol)
  genesGroupDF <- genesGroupDF[! genesGroupDF$Gene.Symbol == "", ]
  genesGroupDF <- genesGroupDF[! is.na(genesGroupDF$Grouping), ]
  genesGroupDF$Grouping <- gsub(" *", "", genesGroupDF$Grouping)
  genesGroupDF$Grouping <- factor(genesGroupDF$Grouping, levels = unique(genesGroupDF$Grouping))

  # Gene groups to plot - add known markers to cell cycle genes table
  genesGroupDF <- rbind(genesGroupDF
    , kmDF[kmDF$Grouping %in% c("RG", "IP", "Neuron", "vRG", "oRG") ,c(3,2)])

  # Feature plot
  ggL <- FeaturePlot(
    genes = genesGroupDF$Gene.Symbol
    , tsneDF = as.data.frame(ssCentSO@dr$tsne@cell.embeddings)
    , seuratO = ssCentSO
    , exM = ssCentSO@scale.data
    , centScale = TRUE
    , limLow = -1.5, limHigh = 1.5
    , geneGrouping = genesGroupDF$Grouping)
  ggL <- lapply(ggL, function(gg){
    gg + geom_point(size = 0.5)
  })

  df2 <- as.data.frame(ssCentSO@dr$tsne@cell.embeddings)
  df2$CellID <- row.names(df2)
  df2 <- Marker_Expression_Flag(df2, df1)
  print(str(df2))
  gg <- ggplot(df2, aes(x = tSNE_1, y = tSNE_2, color = Cell_Subset)) +
    geom_point(size = 0.5)
    # xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
    # ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
  ggL <- append(ggL, list(gg))
  # gg <- ggplot(df2, aes(x = tSNE_1, y = tSNE_2, color = Cell_Subset_075)) +
  #   geom_point()
  #   # xlab(paste0(PCx, " (", round(varExpL[[PCx]]*100, 2), "%)")) +
  #   # ylab(paste0(PCy, " (", round(varExpL[[PCy]]*100, 2), "%)"))
  # ggL <- append(ggL, list(gg))

  pg <- Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.2, align = 'v'
    , axis = 'r'
    , title = paste0(graphCodeTitle
      , "\n\n", title
      , "\nExpression of CC phase genes and cell type marker genes"
      , "\nNormalized centered scaled expression"
      , "\nGene list from Macosko et al. 2016")
    )
    return(pg)
}

# Identify RG+ and / or IP+ but Neuron- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 8])
# tSNE and plot
tSNE_CCgenes_Markers(cellIDs = cellIDs
  , title = "Subset cluster 8 RG+ IP+ Neuron- (>0.5, >0.5, <0.25)"
  , df1 = df1)
ggsave(paste0(outGraph, "RG_IP_Neuronn_Cluster8_FeaturePlot.png")
  , width = 12, height = 12, limitsize = FALSE)

# Identify RG+ and / or IP+ but Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident)[centSO@ident == 10])
# tSNE and plot
tSNE_CCgenes_Markers(cellIDs = cellIDs
  , title = "Subset cluster 8 RG+ IP+ Neuron- (>0.5, >0.5, <0.25)"
  , df1 = df1)
ggsave(paste0(outGraph, "RG_IP_Neuronn_Cluster10_FeaturePlot.png")
  , width = 12, height = 9, limitsize = FALSE)
################################################################################

# ### Heatmap and heirarchical clustering of RG and IP cells
# print("### Heatmap and heirarchical clustering of RG and IP cells")
#
# ## Heatmap and heirarchical clustering of RG and IP G2M cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.5 | df1$IP > 0.5 & df1$Neuron < 0.5]
# cellIDs <- intersect(cellIDs, row.names(df1)[df1$PHASE == "G2M"])
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-2, 2, by = 0.1)
#
# png(paste0(outGraph, "RG_IP_G2M_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
# )
# dev.off()
#
#
# ## Heatmap and heirarchical clustering of RG cluster 10 cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
# df1 <- df1[row.names(df1) %in% cellIDs, ]
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# # annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# # annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-3, 3, by = 0.1)
#
# png(paste0(outGraph, "RG_Cluster10_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
#   , treeheight_col = 150
#   , cutree_cols = 10
# )
# dev.off()
#
#
# ## Heatmap and heirarchical clustering of RG cluster 8 cells
#
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5]
# cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
# df1 <- df1[row.names(df1) %in% cellIDs, ]
# exM <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]
# mns <- rowMeans(exM)
# genes <- names(sort(mns, decreasing = TRUE))[1:2500]
# exM <- exM[row.names(exM) %in% genes, ]
#
# annotation_col <- data.frame(Type = rep(NA, nrow(df1)))
# row.names(annotation_col) <- row.names(df1)
# annotation_col$Type[df1$RG > 0.5] <- "RG+"
# annotation_col$Type[df1$IP > 0.5] <- "IP+"
# annotation_col$Type[df1$IP > 0.5 & df1$Neuron > 0.5] <- "IP+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$Neuron > 0.5] <- "RG+ Neuron+"
# annotation_col$Type[df1$RG > 0.5 & df1$IP > 0.5] <- "RG+ IP+"
# # annotation_col <- annotation_col[df1$PHASE == "G2M", , drop = FALSE]
# # annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]
#
# annotation_col$vRGoRG <- NA
# ids <- row.names(df1)[df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "oRG+"
# ids <- row.names(df1)[df1$vRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+"
# ids <- row.names(df1)[df1$vRG > 0.5 & df1$oRG > 0.5]
# annotation_col$vRGoRG[row.names(annotation_col) %in% ids] <- "vRG+ oRG+"
#
# idx <- match(row.names(annotation_col), row.names(df1))
# annotation_col$G2M <- df1$G2Mscore[idx]
# annotation_col$Sscore <- df1$Sscore[idx]
#
# exM[exM > 3] <- 3
# exM[exM < -3] <- 3
# breaks <- seq(-3, 3, by = 0.1)
#
# png(paste0(outGraph, "RG_Cluster8_pheatmap.png"), width = 9, height = 9, units = "in", res = 300)
# pheatmap(exM,
#   cluster_row = TRUE
#   , cluster_cols = TRUE
#   , annotation_col = annotation_col
#   , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
#   , breaks = breaks
#   , show_rownames = FALSE
#   , show_colnames = FALSE
#   , treeheight_col = 150
#   , cutree_cols = 10
# )
# dev.off()
################################################################################

### PCA of RG+ and / or IP+ cells
print("### PCA of RG+ and / or IP+ cells")

exLM <- list()

# Identify cluster 8 cells
cellIDs <- names(centSO@ident[centSO@ident == 8])
exLM[["Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify cluster 10 cells
cellIDs <- names(centSO@ident[centSO@ident == 10])
exLM[["Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or IP+ but Neuron- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IP_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or IP+ but Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 | df1$IP > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IP_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuron_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ and / or Neuron+ but IP- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 | df1$Neuron > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuron_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+ IP- Neuron- cluster 8 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 8]))
exLM[["RG_IPn_Neuronn_Cluster8"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]

# Identify RG+  IP- Neuron- cluster 10 cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs, names(centSO@ident[centSO@ident == 10]))
exLM[["RG_IPn_Neuronn_Cluster10"]] <- centSO@scale.data[ ,colnames(centSO@scale.data) %in% cellIDs]


# Plot mean expression rank, standard deviation, and variance
pgL <- lapply(names(exLM), function(name) {
  exM <- exLM[[name]]
  MeanExprRank_Stdev_Variance_ScatterPlot(
    exM = exM
    , title = name
  )
})
# Combine plots from each cell subset
Plot_Grid(ggPlotsL = pgL, ncol = 1, rel_height = 0.1
  , title = paste0(
    graphCodeTitle
    , "\n\nSubset cells by RG+ and / or IP+ expression and G2M or S phase"))
ggsave(paste0(outGraph, "RG_IP_MeanExpr_StdDev.png"), width = 13, height = 20)

# Subset to top expressed genes
exLM <- lapply(exLM, function(exM){
  mns <- rowMeans(exM)
  genes <- names(sort(mns, decreasing = TRUE))[1:10000]
  # genes <- names(sort(mns, decreasing = TRUE))[1:1000]
  exM <- exM[row.names(exM) %in% genes, ]
  return(exM)
})
# Subset to high variance genes
exLM <- lapply(exLM, function(exM){
  varDF <- data.frame(apply(exM, 1, var))
  # ggplot(varDF, aes(x = varDF[,1])) +
  #   geom_histogram(binwidth = 0.1)
  # ggsave(paste0(outGraph, "VarHist.png"))
  exM <- exM[row.names(exM) %in% row.names(varDF)[varDF > 0.3], ]
  return(exM)
})

# PCA
pcaL <- lapply(exLM, function(exM){
  pca <- prcomp(t(exM), center = FALSE)
  print(head((pca$sdev)^2 / sum(pca$sdev^2)*100))
  return(pca)
})

# Plot PCA

outGraphPCA <- paste0(
  dirname(outGraph)
  , "/PCA_Top10000_Var03_NoCent/"
  , basename(outGraph)
)
dir.create(dirname(outGraphPCA), recursive = TRUE)

# Plot PCA loadings
lapply(names(pcaL), function(name){
  pca <- pcaL[[name]]
  Prcomp_Loadings_Plot(pca = pca, nGenes = c(1:20), nPCs = c(1:8)
    , title = paste0(graphCodeTitle, "\n\n", name, "\nGenes with highest PC loadings"))
  ggsave(paste0(outGraphPCA, name, "_PCAloadings.pdf"), width = 13
    , height = 20, limitsize = FALSE)
})

# Format for ggplot
ldf <- lapply(pcaL, function(pca){
  PCA_Format_For_GGplot(pca)
})

# mns <- rowMeans(noCentExM[ ,centSO@ident == 7])
# corCluster7 <- apply(noCentExM, 2, function(col){cor(col, mns)})
#
# mns <- rowMeans(noCentExM[ ,centSO@ident == 9])
# corCluster9 <- apply(noCentExM, 2, function(col){cor(col, mns)})
#
# ldf <- lapply(ldf, function(df){
#
#   idx <- match(row.names(df), names(corCluster7))
#   df$Cor_Cluster7 <- corCluster7[idx]
#
#   idx <- match(row.names(df), names(corCluster9))
#   df$Cor_Cluster9 <- corCluster7[idx]
#
#   return(df)
# })

# PCA plots
lapply(names(ldf), function(name){

  pcaL <- ldf[[name]]

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel05.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset_025")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.25 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel025.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Cell_Subset_075")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG+ and / or IP+"
      , "\n+ = > 0.75 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_MarkerLabel075.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG_oRG_Subset")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG+ and / or oRG+"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRGoRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "S_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by S phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Sscore.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "G2M_Score")
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by G2M phase score"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_G2Mscore.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "vRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by vRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_vRGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "oRG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_oRG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "oRG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by oRG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_oRGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "RG", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_RG.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "RG_PollenS3", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by RG Pollen S3 expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_RGPollenS3.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "IP", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_IP.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 1)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Neuron.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "IP", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by IP expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_IP025scale.png"), width = 14, height = 12)

  ggL <- PCA_Plot_PC1to8(pcaL = pcaL, colorBy = "Neuron", limLow = 0, limHigh = 0.25)
  Plot_Grid(ggPlotsL = ggL
    , ncol = 3
    , rel_height = 0.1
    , title = paste0(graphCodeTitle
      , "\n\n", name, " PCA"
      , "\nColored by Neuron expression"
      , "\n+ = > 0.5 log normalized expression")
  )
  ggsave(paste0(outGraphPCA, name, "_PCA_Neuron025scale.png"), width = 14, height = 12)

  # # PCA plots
  # lapply(names(ldf), function(name){
  #
  # df <- ldf[[name]]

  # ggL <- PCA_Plot_PC1to8(pcaDF = df, colorBy = "Cor_Cluster7", limLow = -1, limHigh = 1)
  # Plot_Grid(ggPlotsL = ggL
  #   , ncol = 3
  #   , rel_height = 0.1
  #   , title = paste0(graphCodeTitle
  #     , "\n\n", name, " PCA"
  #     , "\nColored by correlation to Cluster 7 (oRG) mean expression profile"
  #     , "\n+ = > 0.5 log normalized expression")
  # )
  # ggsave(paste0(outGraphPCA, name, "_PCA_CorCluster7.png"), width = 14, height = 12)
  #
  # ggL <- PCA_Plot_PC1to8(pcaDF = df, colorBy = "Cor_Cluster9", limLow = -1, limHigh = 1)
  # Plot_Grid(ggPlotsL = ggL
  #   , ncol = 3
  #   , rel_height = 0.1
  #   , title = paste0(graphCodeTitle
  #     , "\n\n", name, " PCA"
  #     , "\nColored by correlation to Cluster 9 (vRG) mean expression profile"
  #     , "\n+ = > 0.5 log normalized expression")
  # )
  # ggsave(paste0(outGraphPCA, name, "_PCA_CorCluster9.png"), width = 14, height = 12)

})
################################################################################

### Expression of vRG, oRG, RG markers versus cell cycle score
print("### Expression of vRG, oRG, RG markers by cell cycle phase")

# Violin plots of Seurat cell cycle scores
df1 <- centSO@meta.data[c("S.Score", "G2M.Score", "res.0.6")]
df1$res.0.6 <- factor(df1$res.0.6, levels = sort(as.numeric(unique(df1$res.0.6))))
ggplot(df1, aes(x = res.0.6, y = S.Score)) +
  geom_violin() +
  geom_jitter(size = 0.01, alpha = 0.25) +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nS phase score by cluster"))
ggsave(paste0(outGraph, "CellCycleSphaseScore_By_Cluster_Violin.png")
  , width = 7, height = 6)
ggplot(df1, aes(x = res.0.6, y = G2M.Score)) +
  geom_violin() +
  geom_jitter(size = 0.01, alpha = 0.25) +
  xlab("Cluster") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nG2M phase score by cluster"))
ggsave(paste0(outGraph, "CellCycleG2MScore_By_Cluster_Violin.png")
  , width = 7, height = 6)

df <- data.frame(Phase = rep("G0", ncol(noCentExM)))
df$Phase[centSO@meta.data["S.Score"] > 0.2] <- "S phase"
df$Phase[centSO@meta.data["G2M.Score"] > 0.5] <- "G2/M"

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "RG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$RG <- mnEx

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "vRG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$vRG <- mnEx

genes <- kmDF$Gene.Symbol[kmDF$Grouping %in% "oRG"]
exM <- noCentExM[row.names(noCentExM) %in% genes, ]
mnEx <- colMeans(exM)
df$oRG <- mnEx

# Subset to RG+ cells
df <- df[df$RG > 0.5, ]

# Format for ggplot
df <- melt(df)
df$Phase <- factor(df$Phase, levels = c("G0", "S phase", "G2/M"))

# Plot
ggplot(df, aes(x = variable, y = value, fill = Phase)) +
  geom_boxplot() +
  xlab("Marker gene group") +
  ylab("Mean normalized expression") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMean expression of vRG, oRG, or RG markers in RG+ cells by cell cycle"
    , "\nSubset to cells to > 0.5 mean normalized expression of RG markers")
  )
ggsave(paste0(outGraph, "RGexpr_By_CellCycle.png"), width = 7, height = 6)

# Plot
ggplot(df, aes(x = variable, y = value, fill = Phase)) +
  geom_violin() +
  geom_jitter(size = 0.25, alpha = 0.25) +
  xlab("Marker gene group") +
  ylab("Mean normalized expression") +
  ggtitle(paste0(graphCodeTitle
    , "\n\nMean expression of vRG, oRG, or RG markers in RG+ cells by cell cycle"
    , "\nSubset to cells to > 0.5 mean normalized expression of RG markers")
  )
ggsave(paste0(outGraph, "RGexpr_By_CellCycle_Violin.png"), width = 7, height = 6)
################################################################################

### DE RG vs IP vs IP+ RG+ analysis

## Heatmap of RG vs IP DE genes in IP+ and Neuron- or RG+ and Neuron- cells

# S phase cells

# Subset expression matrix to DE genes and IP+ and Neuron- or RG+ and Neuron- cells
# DE genes
genes <- row.names(deDF)[deDF$FDR < 0.05]
# IP+ and Neuron- or RG+ and Neuron- cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
ids <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$Neuron < 0.25
  | df1$CLUSTER == 8 & df1$IP > 0.5 & df1$Neuron < 0.25]
# Subset
exM <- as.matrix(centSO@scale.data)
exM <- exM[row.names(exM) %in% genes, colnames(exM) %in% ids]

# Heatmap column color bars
df2 <- df1[row.names(df1) %in% ids, ]
annotation_col <- data.frame(Type = rep(NA, nrow(df2)))
row.names(annotation_col) <- row.names(df2)
annotation_col$Type[df2$RG > 0.5] <- "RG+"
annotation_col$Type[df2$IP > 0.5] <- "IP+"
annotation_col$Type[df2$RG > 0.5 & df2$IP > 0.5] <- "RG+ IP+"
annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

# Heatmap row color bars
annotation_row <- data.frame(DE = rep(NA, nrow(exM)))
row.names(annotation_row) <- row.names(exM)
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC > 0]] <- "RG"
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC < -0]] <- "IP"

breaks <- seq(-2, 2, by = 0.1)

png(paste0(outGraph, "DE_RG_IP_S_pheatmap.png")
  , width = 9, height = 9, units = "in", res = 300)
pheatmap(exM,
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , annotation_row = annotation_row
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
)
dev.off()

# G2M phase cells

# Subset expression matrix to DE genes and IP+ and Neuron- or RG+ and Neuron- cells
# DE genes
genes <- row.names(deDF)[deDF$FDR < 0.05]
# IP+ and Neuron- or RG+ and Neuron- cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
ids <- row.names(df1)[df1$CLUSTER == 10 & df1$RG > 0.5 & df1$Neuron < 0.25
  | df1$CLUSTER == 10 & df1$IP > 0.5 & df1$Neuron < 0.25]
# Subset
exM <- as.matrix(centSO@scale.data)
exM <- exM[row.names(exM) %in% genes, colnames(exM) %in% ids]

# Heatmap column color bars
df2 <- df1[row.names(df1) %in% ids, ]
annotation_col <- data.frame(Type = rep(NA, nrow(df2)))
row.names(annotation_col) <- row.names(df2)
annotation_col$Type[df2$RG > 0.5] <- "RG+"
annotation_col$Type[df2$IP > 0.5] <- "IP+"
annotation_col$Type[df2$RG > 0.5 & df2$IP > 0.5] <- "RG+ IP+"
annotation_col <- annotation_col[! is.na(annotation_col$Type), , drop = FALSE]

# Heatmap row color bars
annotation_row <- data.frame(DE = rep(NA, nrow(exM)))
row.names(annotation_row) <- row.names(exM)
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC > 0]] <- "RG"
annotation_row$DE[row.names(annotation_row) %in% row.names(deDF)[
  deDF$FDR < 0.05 & deDF$Log_FC < -0]] <- "IP"

breaks <- seq(-2, 2, by = 0.1)

png(paste0(outGraph, "DE_RG_IP_G2M_pheatmap.png")
  , width = 9, height = 9, units = "in", res = 300)
pheatmap(exM,
  cluster_row = TRUE
  , cluster_cols = TRUE
  , annotation_col = annotation_col
  , annotation_row = annotation_row
  , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks))
  , breaks = breaks
  , show_rownames = FALSE
  , show_colnames = FALSE
)
dev.off()


## DE RG+ vs RG+ IP+ S phase

DE_CellGroup1_Vs_CellGroup2 <- function(ids1, ids2) {

  # Subset expression matrix
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% c(ids1, ids2)]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% c(ids1, ids2)
    , c("nUMI", "librarylab", "individual", "res.0.6")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$groups <- "ids1"
  termsDF$groups[row.names(termsDF) %in% ids2] <- "ids2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ groups+nUMI+librarylab+individual")

  print(head(deLM$coefmat))

  # Format LM DE
  deDF <- data.frame(Log_FC = deLM$coefmat[ ,2]
    , Pvalue = deLM$pvalmat[ ,2])
  deDF <- deDF[order(deDF$Log_FC), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1

  # FDR correct
  # NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
  corrected <- fdrtool(deDF$Pvalue, statistic = "pvalue", plot = FALSE)
  deDF$FDR <- corrected$lfdr
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  print(head(deDF))

  return(deDF)
}

# Mean expression of marker genes to use to subset groups of cells
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)

# S phase RG+ Neuron- vs RG+ IP+ Neuron-
ids1 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron < 0.25]
ids2 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron < 0.25]
# DE
deRGvsRgIpDF <- DE_CellGroup1_Vs_CellGroup2(ids1, ids2)
# Save as csv
write.csv(deRGvsRgIpDF, file = paste0(outTable, "DE_RGvsRGIP_S.csv")
  , quote = FALSE)
deRGvsRgIpDF <- read.csv(paste0(outTable, "DE_RGvsRGIP_S.csv")
  , header = TRUE, row.names = 1)


# S phase RG+ Neuron- vs RG+ IP+ Neuron-
ids1 <- row.names(df1)[df1$CLUSTER == 8 & df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron < 0.25]
ids2 <- row.names(df1)[df1$CLUSTER == 8 & df1$IP > 0.5 & df1$RG < 0.25 & df1$Neuron < 0.25]
# DE
deIPvsRgIpDF <- DE_CellGroup1_Vs_CellGroup2(ids1, ids2)
# Save as csv
write.csv(deIPvsRgIpDF, file = paste0(outTable, "DE_IPvsRGIP_S.csv")
  , quote = FALSE)
deIPvsRgIpDF <- read.csv(paste0(outTable, "DE_IPvsRGIP_S.csv")
  , header = TRUE, row.names = 1)


intersect(
  row.names(deDF)[deDF$Log_FC > 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC > 0 & deDF$FDR < 0.05]
)

intersect(
  row.names(deDF)[deDF$Log_FC < 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC < 0 & deDF$FDR < 0.05]
)

intersect(
  row.names(deDF)[deDF$Log_FC > 0 & deDF$FDR < 0.05]
  , row.names(deRGvsRgIpDF)[deRGvsRgIpDF$Log_FC < 0 & deDF$FDR < 0.05]
)



intersect(names(deLM1$coefmat[ ,"clustersclusters2"])[deLM1$coefmat[ ,"clustersclusters2"] < 0.5]
, names(deLM$coefmat[,"groupsids2"])[deLM$coefmat[,"groupsids2"] > 0.5])
################################################################################

### S phase gene expression

# Sort genes by mean expression
mnEx <- rowMeans(noCentExM[row.names(noCentExM) %in% sGenes, ])
mnEx <- sort(mnEx, decreasing = TRUE)
genes <- names(mnEx)

# Number of cells each marker expressed in
m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
df1 <- melt(m1)
df1$Cluster <- centSO@ident[match(df1$Var2, names(centSO@ident))]
df1 <- aggregate(value~Var1+Cluster, df1, sum)
df1$Var1 <- factor(df1$Var1, levels = genes)

# Plot number of cells each marker expressed in
ggplot(df1, aes(x = Cluster, y = value)) +
  facet_wrap(~Var1, scales ="free_x", ncol = 3) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cells expressing genes"
    , "\n(> 0.5 normalized expression)"))
ggsave(paste0(outGraph, "Sgenes_Number_Barplot.png")
  , width = 12, height = 34)

# Feature plot
ggL <- FeaturePlot(
  genes = genes
  , tsneDF = as.data.frame(centSO@dr$tsne@cell.embeddings)
  , seuratO = centSO
  , exM = noCentExM
  , limLow = -1, limHigh = 3)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.05, align = 'v', axis = 'r'
  , title = paste0(graphCodeTitle
    , "\n\nExpression of S phase genes sorted by mean expression"
    , "\nNormalized expression"
    , "\nGene list from Tirosh et al. 2016"))
ggsave(paste0(outGraph, "Sgenes_FeaturePlot.png"), width = 20, height = 70
  , limitsize = FALSE)
################################################################################

### Intersection of S phase genes, RG markers, IP markers for FISH probes

## tSNE colored by intersection and heatmap of numbers of intersections

# Genes to intersect
genes <- c("PCNA", "EOMES", "SOX2", "PAX6", "HOPX", "CRYAB")

# Number of cells each marker expressed in
m1 <- noCentExM[row.names(noCentExM) %in% genes, ] > 0.5
df1 <- melt(m1)
df1$Cluster <- centSO@ident[match(df1$Var2, names(centSO@ident))]
df1 <- aggregate(value~Var1+Cluster, df1, sum)

# Plot number of cells each marker expressed in
ggplot(df1, aes(x = Cluster, y = value)) +
  facet_wrap(~Var1, scales ="free_x") +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste0(graphCodeTitle
    , "\n\nNumber of cells expressing genes"
    , "\n(> 0.5 normalized expression)"))
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Number_Barplot.png")
  , width = 12, height = 10)

# tSNE
ggL <- Intersection_tSNE_Plots(genes)
gg1 <- TSNE_Plot(centSO) + theme(legend.position = "none")
ggL <- append(list(gg1), ggL)
Plot_Grid(ggPlotsL = ggL, ncol = 3, rel_height = 0.05, align = 'v', axis = 'r'
  , title = paste0(paste0(graphCodeTitle
    , "\n\ntSNE plot colored by intersection of expression of gene A and gene B"
    , "\n(> 0.5 normalized expression)"))
)
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Intersection_tSNE.png")
  , width = 20, height = 38)

# Heatmap
Number_Of_Cells_Intersection_Heatmap(
  genes = genes
  , title = paste0(graphCodeTitle
    , "\n\nNumber of cells expressing both gene A and gene B"
    , "\n(> 0.5 normalized expression)")
)
ggsave(paste0(outGraph, "RG_IP_Sphase_Marker_Intersection_Heatmap.png")
  , width = 7, height = 7)
################################################################################

### Number of cells positive for RG, IP, or Neuron markers

df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)

df <- data.frame(Cell_Subset = rep(NA, nrow(df1)))
row.names(df) <- row.names(df1)

df$Cell_Subset <- NA
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$Neuron > 0.5]] <- "Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5]] <- "RG"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5]] <- "IP"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$Neuron > 0.5]] <- "RG Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5]] <- "RG IP"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$IP > 0.5 & df1$Neuron > 0.5]] <- "IP Neuron"
df$Cell_Subset[rownames(df) %in% row.names(df1)[df1$RG > 0.5 & df1$IP > 0.5 & df1$Neuron > 0.5]] <- "RG IP Neuron"

table(df)
################################################################################
