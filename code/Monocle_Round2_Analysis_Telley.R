# Damon Polioudakis
# 2017-08-15
# Run monocle
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
toOrder <- c("0-1-4-12", "0-1-2", "0-1", "3-14", "5-6", "7-9", "0", "1", "2", "3"
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

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# Molyneaux markers
mmDF <- read.csv("../source/Molyneaux_LayerMarkers_Format.csv", header = TRUE)

# Telley transcriptional waves
tyDF <- read.csv("../source/Telley_2015_ST3_TranscriptionalWaves.csv", header = TRUE)

# Human TFs
tfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_TF_EnsemblID.txt")
# Human chromatin remodeling factors
crDF <- read.table(
  "../source/AnimalTFDB_Homo_sapiens_chr_remodeling_factor_EnsemblID.txt")
# Human co-factors
cfDF <- read.table("../source/AnimalTFDB_Homo_sapiens_cofactor_EnsemblID.txt")

# biomaRt gene info
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
graphCodeTitle <- "Monocle_Round2_Analysis_Telley.R"
outGraph <- "../analysis/graphs/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_"
outTable <- "../analysis/tables/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_"
outRdat <- "../analysis/Monocle_Round2_Analysis_Telley/Comp1-10/Monocle_Round2_"

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



Convert_MouseGenes_To_HumanGenes <- function(geneList) {
  geneList <- data.frame(geneList)
  mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
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


# Output Directories
outDir <- paste0(outGraph, "Pseudotime_Heatmap_TelleyCore-P-Ng-Ngenes")
dir.create(outDir, recursive = TRUE)
outSubGraph <- paste0(outDir, "/", basename(outDir), "_")

lapply(names(ptDeL)[1:4], function(cluster){
  
  # cluster <- "0-1"
  print(cluster)
  
  diff_test_res <- ptDeL[[cluster]]
  mo_filtered <- moL[[cluster]]
  print(str(mo_filtered))
  
  P <- c("Fabp7", "Nes", "Pax6", "Slc1a3", "Sox2", "Vim")
  Ng <- c("Btg2", "Eomes", "Neurog2", "Tuba1a")
  N <- c("Dcx", "Neurod1", "Map2", "Mapt", "Tbr1", "Tubb3", "Rbfox3")
  
  P <- Convert_MouseGenes_To_HumanGenes(P)
  P <- P$hsapiens_homolog_associated_gene_name
  
  Ng <- Convert_MouseGenes_To_HumanGenes(Ng)
  Ng <- Ng$hsapiens_homolog_associated_gene_name
  
  N <- Convert_MouseGenes_To_HumanGenes(N)
  N <- N$hsapiens_homolog_associated_gene_name
  
  telleyWaveL <- lapply(tyDF, function(df) {
    df <- Convert_MouseGenes_To_HumanGenes(df)
    return(df$hsapiens_homolog_associated_gene_name)
  })
  
  
  
  GenesList_Expression_By_Pseudotime <- function(geneList, exM) {
    # df1 <- melt(as.matrix(centSO@scale.data)[row.names(centSO@scale.data) %in% P, ])
    
    df1 <- melt(exM[row.names(exM) %in% geneList, ])
    idx <- match(df1$Var2, row.names(mo_filtered@phenoData@data))
    df1$Pseudotime <- mo_filtered@phenoData@data$Pseudotime[idx]
    df1 <- aggregate(value~Var2+Pseudotime, df1, mean)
    # Downsample for testing
    # df1 <- df1[sample(1:nrow(df1), 1000), ]
    return(df1)
  }
  
  Expression_By_Pseudotime_FitLine_Plot <- function(df){
    gg <- ggplot(df, aes(x = Pseudotime, y = value)) +
      geom_jitter(alpha = 0.1, size = 0.1) +
      geom_smooth(method = "auto")
    # stat_summary(fun.y = "mean", color = "red", size = 1, geom = "point")
    return(gg)
  }
  
  ggDF <- GenesList_Expression_By_Pseudotime(geneList = P, exM = noCentExM)
  Expression_By_Pseudotime_FitLine_Plot(ggDF)
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_P.png"))
  
  ggDF <- GenesList_Expression_By_Pseudotime(geneList = Ng, exM = noCentExM)
  Expression_By_Pseudotime_FitLine_Plot(ggDF)
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_Ng.png"))
  
  ggDF <- GenesList_Expression_By_Pseudotime(geneList = N, exM = noCentExM)
  Expression_By_Pseudotime_FitLine_Plot(ggDF)
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_N.png"))
  
  ggL <- lapply(telleyWaveL, function(genes) {
    ggDF <- GenesList_Expression_By_Pseudotime(genes, exM = noCentExM)
    gg <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
    return(gg)
  })
  Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_Waves.png")
    , width = 12, height = 7)
  
  ggL <- lapply(telleyWaveL, function(genes) {
    ggDF <- GenesList_Expression_By_Pseudotime(genes, exM = as.matrix(centSO@scale.data))
    gg <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
    return(gg)
  })
  Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_Waves_CenteredScaled.png")
    , width = 12, height = 7)
  
  ggL <- lapply(telleyWaveL, function(genes) {
    de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.05]
    genes <- genes[genes %in% de_genes]
    ggDF <- GenesList_Expression_By_Pseudotime(genes, exM = noCentExM)
    gg <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
    return(gg)
  })
  Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_Waves_PseudotimeDE.png")
    , width = 12, height = 7)
  
  ggL <- lapply(telleyWaveL, function(genes) {
    de_genes <- diff_test_res$gene_short_name[diff_test_res$qval < 0.05]
    genes <- genes[genes %in% de_genes]
    ggDF <- GenesList_Expression_By_Pseudotime(genes, exM = as.matrix(centSO@scale.data))
    gg <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
    return(gg)
  })
  Plot_Grid(ggL, ncol = 3, rel_height = 0.2, title = "Telley Waves")
  ggsave(paste0(outSubGraph, "Cluster", paste0(cluster, collapse = "-"), "_Waves_PseudotimeDE_CenteredScaled.png")
    , width = 12, height = 7)
  
  
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[1]])
  # gg1 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[2]])
  # gg2 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[3]])
  # gg3 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[4]])
  # gg4 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[5]])
  # gg5 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # ggDF <- GenesList_Expression_By_Pseudotime(telleyWaveL[[6]])
  # gg6 <- Expression_By_Pseudotime_FitLine_Plot(ggDF)
  # 
  # 
  png(paste0(outSubGraph, "facet_N.png"))
  plot_genes_in_pseudotime(mo_filtered[N[1:6], ])
  dev.off()
  
  png(paste0(outSubGraph, "facet_P.png"))
  plot_genes_in_pseudotime(mo_filtered[P, ])
  dev.off()
  
  png(paste0(outSubGraph, "facet_Ng.png"))
  plot_genes_in_pseudotime(mo_filtered[Ng, ])
  dev.off()
})
  
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