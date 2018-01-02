# Damon Polioudakis
# 2017-12-09
# Analyize Monocle output using Telley gene lists
################################################################################

rm(list = ls())

require(monocle)
require(Seurat)
require(cowplot)
require(ggplot2)
require(viridis)
require(reshape2)
require(biomaRt)
source("Function_Library.R")

# qsub task ID for selecting Seurat cluster ID to run through Monocle
args <- commandArgs(trailingOnly = TRUE)

## Inputs

# Monocle round 2

# How to reorder list of monocle objects for plotting
toOrder <- c("0-1-2-4-12", "0-1-4-12", "0-1-2", "0-1", "3-14", "5-6", "7-9", "0", "1", "2", "3"
  , "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")

# Monocle object
inMo <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inMo <- inMo[grep("Monocle_Round2_monocleO", inMo)]
moL <- lapply(inMo, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(mo_filtered)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_monocleO_cluster", "", inMo)
names <- gsub(".Robj", "", names)
names(moL) <- names
# Reorder list of monocle objects for plotting
moL <- moL[toOrder]

# Pseudotime DE
inPtDE <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inPtDE <- inPtDE[grep("Monocle_Round2_DEpseudotime_cluster", inPtDE)]
ptDeL <- lapply(inPtDE, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(ptDE)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEpseudotime_cluster", "", inPtDE)
names <- gsub(".Robj", "", names)
names(ptDeL) <- names
# Reorder list of DE for plotting
ptDeL <- ptDeL[toOrder]

# State DE
inStDE <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inStDE <- inStDE[grep("DEstate_cluster", inStDE)]
# inStDE <- inStDE[8:15]
stDeLDF <- lapply(inStDE, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  # load("../analysis/Monocle_Round2/Monocle_Round2_Comp1-20_DEstate_cluster5-6.Robj")
  return(stateDE)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEstate_cluster", "", inStDE)
names <- gsub(".Robj", "", names)
names(stDeLDF) <- names
# Reorder list of DE for plotting
stDeLDF <- stDeLDF[toOrder]

# Branch DE
inBeam <- list.files("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/")
inBeam <- inBeam[grep("DEbranch_cluster", inBeam)]
beamLLDF <- lapply(inBeam, function(path) {
  load(paste0("../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/", path))
  return(beamLDF)
})
# List names are Seurat clusters used for Monocle round 2
names <- gsub("Monocle_Round2_DEbranch_cluster", "", inBeam)
names <- gsub(".Robj", "", names)
names(beamLLDF) <- names
# Reorder list of DE for plotting
beamLLDF <- beamLLDF[toOrder]

# Seurat clustering
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM

# Telley transcriptional waves
tyDF <- read.csv("../source/Telley_2015_ST3_TranscriptionalWaves.csv", header = TRUE)

## Variables
graphCodeTitle <- "Monocle_Round2_Analysis_Telley.R"
outGraph <- "../analysis/graphs/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_Analysis_Telley_"
outTable <- "../analysis/tables/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_Analysis_Telley_"
outRdat <- "../analysis/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_Analysis_Telley_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outRdat), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Functions

Convert_MouseGenes_To_HumanGenes <- function(geneList) {
  geneList <- data.frame(geneList)
  mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "aug2017.archive.ensembl.org")
  # mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
  mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
  # Look up human orthologs
  ensemblMmHsDF <- getBM(
    attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"
      , "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_perc_id")
    , filters = "mgi_symbol"
    , values = geneList
    , mart = mart
  )
  # Have to add MGI symbol separately, can't call BM to add both ortholog and MGI
  # at the same time
  ensemblMgiSymDF <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol")
    , filters = "mgi_symbol"
    , values = geneList
    , mart = mart
  )
  # Merge mouse gene symbols and human
  genesDF <- merge(ensemblMmHsDF, ensemblMgiSymDF
    , by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
  return(genesDF)
}

GenesList_Expression_By_Pseudotime <- function(
  gene_group_DF
  , exM
  , mo_filtered
  , center_scale = FALSE) {
  
  # Center and scale expression for each gene across cells
  if (center_scale == TRUE) {
    exM <- t(scale(t(exM)))
  }
  
  # Add pseudotime and calculate mean gene expression for each cell
  df1 <- melt(exM[row.names(exM) %in% gene_group_DF$Gene, ])
  idx <- match(df1$Var2, row.names(mo_filtered@phenoData@data))
  df1$Pseudotime <- mo_filtered@phenoData@data$Pseudotime[idx]
  df1$Group <- gene_group_DF$Group[match(df1$Var1, gene_group_DF$Gene)]
  df1 <- aggregate(value~Var2+Pseudotime+Group, df1, mean)
  # Downsample for testing
  # df1 <- df1[sample(1:nrow(df1), 1000), ]
  return(df1)
}

Expression_By_Pseudotime_FitLine_Plot <- function(df){
  gg <- ggplot(df, aes(x = Pseudotime, y = value)) +
    facet_wrap(~Group) +
    geom_jitter(alpha = 0.1, size = 0.1) +
    geom_smooth(method = "auto", aes(color = Group)) +
    theme(legend.position="none") +
    ylab("Expression")
  # stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point")
  return(gg)
}

Expression_By_Pseudotime_FitLine_Combo_Plot <- function(df){
  gg <- ggplot(df, aes(x = Pseudotime, y = value, group = Group)) +
    geom_smooth(method = "auto", aes(color = Group)) +
    ylab("Expression")
  # stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point")
  return(gg)
}
################################################################################

### Plot Telley genes across pseudotime

# Output Directories
outDir <- paste0(outGraph, "Core-P-Ng-Ngenes_Loess")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_")

lapply(names(ptDeL)[1:4], function(cluster){
  
  # cluster <- "0-1-2"
  print(cluster)
  
  diff_test_res <- ptDeL[[cluster]]
  mo_filtered <- moL[[cluster]]
  
  # Telley Proliferative, Neurogenic, Neuronal genes
  telleyPNgNgenesL <- list(
    P = c("Fabp7", "Nes", "Pax6", "Slc1a3", "Sox2", "Vim")
    , Ng = c("Btg2", "Eomes", "Neurog2", "Tuba1a")
    , N = c("Dcx", "Neurod1", "Map2", "Mapt", "Tbr1", "Tubb3", "Rbfox3")  
  )
  # Human homolog
  telleyPNgNgenesL <- lapply(telleyPNgNgenesL, function(genes) {
    df1 <- Convert_MouseGenes_To_HumanGenes(genes)
    return(df1$hsapiens_homolog_associated_gene_name)
  })
    

  ## Plot Telley Proliferative, Neurogenic, Neuronal genes
  
  # Telley P, Ng, N genes filtered by pseudotime DE FDR < 0.001
  gene_group_DF <- melt(telleyPNgNgenesL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.001]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
    )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.001"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE001_Facet_CenterScaled.png")
    , width = 7, height = 5)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.001"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE001_CenterScaled.png")
    , width = 7, height = 5)
  
  # Telley core Proliferative, Neurogenic, Neuronal genes filtered by pseudotime DE FDR < 0.01
  gene_group_DF <- melt(telleyPNgNgenesL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.01]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.01"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE01_Facet_CenterScaled.png")
    , width = 7, height = 5)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.01"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE01_CenterScaled.png")
    , width = 7, height = 5)
  
  # Telley core Proliferative, Neurogenic, Neuronal genes filtered by pseudotime DE FDR < 0.05
  gene_group_DF <- melt(telleyPNgNgenesL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.05]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.05"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE05_Facet_CenterScaled.png")
    , width = 7, height = 5)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nDE across pseudotime filter: FDR < 0.05"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_PseudotimeDE05_CenterScaled.png")
    , width = 7, height = 5)
  
  # Telley core Proliferative, Neurogenic, Neuronal genes
  gene_group_DF <- melt(telleyPNgNgenesL)
  colnames(gene_group_DF) <- c("Gene", "Group")

  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_Facet_CenterScaled.png")
    , width = 7, height = 5)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley core Proliferative, Neurogenic, Neuronal genes"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_PNgNgenes_CenterScaled.png")
    , width = 7, height = 5)
  
  # png(paste0(outSubGraph, "facet_N.png"))
  # plot_genes_in_pseudotime(mo_filtered[N[1:6], ])
  # dev.off()
  # 
  # png(paste0(outSubGraph, "facet_P.png"))
  # plot_genes_in_pseudotime(mo_filtered[P, ])
  # dev.off()
  # 
  # png(paste0(outSubGraph, "facet_Ng.png"))
  # plot_genes_in_pseudotime(mo_filtered[Ng, ])
  # dev.off()
})
################################################################################

### Telley wave genes

# Output Directories
outDir <- paste0(outGraph, "WaveGenes")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_")

lapply(names(ptDeL)[1:4], function(cluster){
  
  # cluster <- "0-1-2"
  print(cluster)
  
  diff_test_res <- ptDeL[[cluster]]
  mo_filtered <- moL[[cluster]]
  
  # Wave genes human homolog
  telleyWaveL <- lapply(tyDF, function(df) {
    df <- Convert_MouseGenes_To_HumanGenes(df)
    return(df$hsapiens_homolog_associated_gene_name)
  })
  
  
  ## Plot Telley wave genes
  
  # Telley wave genes filtered by pseudotime DE FDR < 0.001
  gene_group_DF <- melt(telleyWaveL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.001]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.001"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE001_Facet_CenterScaled.png")
    , width = 7, height = 7)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.001"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE001_CenterScaled.png")
    , width = 7, height = 7)
  
  # Telley wave genes filtered by pseudotime DE FDR < 0.01
  gene_group_DF <- melt(telleyWaveL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.01]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.01"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE01_Facet_CenterScaled.png")
    , width = 7, height = 7)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.01"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE01_CenterScaled.png")
    , width = 7, height = 7)
  
  # Telley wave genes filtered by pseudotime DE FDR < 0.05
  gene_group_DF <- melt(telleyWaveL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.05]
  gene_group_DF <- gene_group_DF[gene_group_DF$Gene %in% de_genes, ]
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.05"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE05_Facet_CenterScaled.png")
    , width = 7, height = 7)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nDE across pseudotime filter: FDR < 0.05"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_PseudotimeDE05_CenterScaled.png")
    , width = 7, height = 7)
  
  # Telley wave genes
  gene_group_DF <- melt(telleyWaveL)
  colnames(gene_group_DF) <- c("Gene", "Group")
  
  # Mean expression of genes per cell and format for ggplot
  ggDF <- GenesList_Expression_By_Pseudotime(
    gene_group_DF = gene_group_DF
    , exM = noCentExM
    , mo_filtered = mo_filtered
    , center_scale = TRUE
  )
  
  # Faceted graph
  Expression_By_Pseudotime_FitLine_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nCentered and scaled expression"))
  # Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_Facet_CenterScaled.png")
    , width = 7, height = 7)
  
  # Same graph
  Expression_By_Pseudotime_FitLine_Combo_Plot(ggDF) +
    ggtitle(paste0(graphCodeTitle
      , "\n\nExpression of Telley transcriptional wave genes"
      , "\nCentered and scaled expression"))
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-")
    , "_Waves_CenterScaled.png")
    , width = 7, height = 7)
})
################################################################################
    # 
  # 
  # 
  # ## Heatmaps of Telley transcriptional wave genes
  # 
  # # Output Directories
  # outDir <- paste0(outGraph, "Pseudotime_Heatmap_TelleyWaves")
  # dir.create(outDir, recursive = TRUE)
  # outSubGraph <- paste0(outDir, "/", basename(outDir), "_")
  # # Plot
  # lapply(names(ptDeL), function(cluster){
  #   
  #   diff_test_res <- ptDeL[[cluster]]
  #   mo_filtered <- moL[[cluster]]
  #   
  #   Convert_MouseGenes_To_HumanGenes <- function(geneList) {
  #     geneList <- data.frame(geneList)
  #     mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  #     mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
  #     mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
  #     # Look up human orthologs
  #     ensemblMmHsDF <- getBM(
  #       attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"
  #         , "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_perc_id")
  #       , filters = "mgi_symbol"
  #       , values = geneList
  #       , mart = mart
  #     )
  #     # Have to add MGI symbol separately, can't call BM to add both ortholog and MGI
  #     # at the same time
  #     ensemblMgiSymDF <- getBM(
  #       attributes = c("ensembl_gene_id", "mgi_symbol")
  #       , filters = "mgi_symbol"
  #       , values = geneList
  #       , mart = mart
  #     )
  #     # Merge mouse gene symbols and human
  #     genesDF <- merge(ensemblMmHsDF, ensemblMgiSymDF
  #       , by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
  #     return(genesDF)
  #   }
  #   
  #   lapply(colnames(tyDF), function(wave) {
  #     
  #     telleyGenes <- tyDF[ ,wave]
  #     
  #     # Convert mouse to human
  #     df1 <- Convert_MouseGenes_To_HumanGenes(geneList = telleyGenes)
  #     telleyGenes <- df1$hsapiens_homolog_associated_gene_name
  #     
  #     # Subset to TFs, co-factors, chromatin remodelers and by qval
  #     diff_test_res <- diff_test_res[diff_test_res$qval < 0.05
  #       & diff_test_res$gene_short_name %in% telleyGenes, ]
  #     
  #     # Order by qval
  #     diff_test_res <- diff_test_res[order(diff_test_res$qval), ]
  #     
  #     genes <- row.names(diff_test_res)
  #     
  #     png(paste0(outSubGraph, "Cluster", cluster, "_", wave, ".png")
  #       , width = 12, height = 16, units = "in", res = 300)
  #     tryCatch(plot_pseudotime_heatmap(mo_filtered[genes,]
  #       , num_clusters = 3
  #       , cores = 1
  #       , show_rownames = T)
  #       , error = function(cond) {
  #         message(paste0("Error for cluster: ", cluster))
  #         message("Error message:")
  #         message(cond)
  #         # Choose a return value in case of error
  #         return(NULL)
  #       }
  #     )
  #     dev.off()
  #     
  #   })
  # })