# Damon Polioudakis
# 2017-11-01
# Separate cycling vRG and oRG

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

# problem is at Format_GenesGroups_Pseudotime_For_GGplot

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(reshape2)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(viridis)
require(monocle)
source("Function_Library.R")
source("Seurat_Cluster_Cycling_vRGoRG_Functions.R")

options(stringsAsFactors = FALSE)

# Keep CC genes from variable gene list used for clustering
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TESTcluster0278910_seuratO.Robj")
centSO <- ssCentSO; rm(ssCentSO)
noCentExM <- ssNoCentExM; rm(ssNoCentExM)

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

### Functions

################################################################################

### DE RG(cluster 7,9) vs IP(2); RG(7,9) vs Neuron(0); IP(2) vs Neuron(0)

DE_Clusters_Vs_Clusters <- function(clusters1, clusters2) {
  ids <- names(centSO@ident)[centSO@ident %in% c(clusters1, clusters2)]
  exM <- as.matrix(centSO@data)
  exM <- exM[ ,colnames(exM) %in% ids]
  # DE Linear model
  termsDF <- centSO@meta.data[
    row.names(centSO@meta.data) %in% ids
    , c("nUMI", "librarylab", "individual", "res.0.6")]
  # Add term TRUE/FALSE cell is in cluster
  termsDF$clusters <- "clusters1"
  termsDF$clusters[termsDF$res.0.6 %in% clusters2] <- "clusters2"
  deLM <- DE_Linear_Model(
    exDatDF = exM
    , termsDF = termsDF
    , mod = "y ~ clusters+nUMI+librarylab+individual")
  # Format LM DE
  deDF <- data.frame(Log_FC_C2vC1 = deLM$coefmat[ ,"clustersclusters2"]
    , Pvalue = deLM$pvalmat[ ,"clustersclusters2"])
  deDF$Gene <- row.names(deDF)
  deDF <- deDF[order(deDF$Log_FC), ]
  deDF$Pvalue[deDF$Pvalue == "NaN"] <- 1
  # FDR correct
  deDF$FDR <- p.adjust(deDF$Pvalue, method="BH")
  # Check
  table(deDF$Pvalue < 0.05)
  table(deDF$FDR < 0.05)
  return(deDF)
}

# de_RG_v_IP_DF <- DE_Clusters_Vs_Clusters(clusters1 = 2, clusters2 = c(7,9))
# # Save as csv
# write.csv(de_RG_v_IP_DF, file = paste0(outTable, "DE_RGvsIP.csv")
#   , quote = FALSE)
#
# de_RG_v_Ne_DF <- DE_Clusters_Vs_Clusters(clusters1 = 0, clusters2 = c(7,9))
# # Save as csv
# write.csv(de_RG_v_Ne_DF, file = paste0(outTable, "DE_RGvsNeuron.csv")
#   , quote = FALSE)
#
# de_IP_v_Ne_DF <- DE_Clusters_Vs_Clusters(clusters1 = 0, clusters2 = c(2))
# # Save as csv
# write.csv(de_IP_v_Ne_DF, file = paste0(outTable, "DE_IPvsNeuron.csv")
#   , quote = FALSE)

# Load DE tables
de_RG_v_IP_DF <- read.csv(paste0(outTable, "DE_RGvsIP.csv")
  , header = TRUE, row.names = 1)
de_RG_v_Ne_DF <- read.csv(paste0(outTable, "DE_RGvsNeuron.csv")
  , header = TRUE, row.names = 1)
de_IP_v_Ne_DF <- read.csv(paste0(outTable, "DE_IPvsNeuron.csv")
  , header = TRUE, row.names = 1)

de_RG_v_IP_DF$Gene <- row.names(de_RG_v_IP_DF)
de_RG_v_Ne_DF$Gene <- row.names(de_RG_v_Ne_DF)
de_IP_v_Ne_DF$Gene <- row.names(de_IP_v_Ne_DF)
de_RG_v_IP_DF$Comparison <- "RG_vs_IP"
de_RG_v_Ne_DF$Comparison <- "RG_vs_Neuron"
de_IP_v_Ne_DF$Comparison <- "IP_vs_Neuron"
deDF <- rbind(de_RG_v_IP_DF, de_RG_v_Ne_DF, de_IP_v_Ne_DF)
################################################################################

### Monocle

Monocle_Run <- function(cellIDs, max_components) {

  print("Monocle: running...")

  # Subset to cells of interest
  exDF <- centSO@raw.data[ ,colnames(centSO@raw.data) %in% cellIDs]
  metDF <- centSO@meta.data[row.names(centSO@meta.data) %in% cellIDs, ]
  # Sort order of cells in exDF to match metadata
  exDF <- exDF[ ,match(row.names(metDF), colnames(exDF))]

  # Format data for monocle
  feature_data = data.frame(gene_short_name = rownames(exDF))
  rownames(feature_data) = feature_data$gene_short_name
  pd <- new("AnnotatedDataFrame", data = metDF)
  fd <- new("AnnotatedDataFrame", data = feature_data)

  # Initialize monocle object
  ssMo <- newCellDataSet(cellData = as(as.matrix(exDF), "sparseMatrix"),
    phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size())

  print("Monocle: Metadata dimensions...")
  print(dim(pData(ssMo)))

  # Filter genes and cells
  # retains all
  ssMo <- detectGenes(ssMo, min_expr = 0.1)
  print(dim(fData(ssMo)))
  # genes expressed in at least 10 cells
  expressed_genes <- row.names(subset(fData(ssMo), num_cells_expressed >= 10))
  ssMo <- ssMo[expressed_genes, ]
  print(dim(fData(ssMo)))
  pData(ssMo)$Total_mRNAs <- Matrix::colSums(exprs(ssMo))

  print("Monocle: Size factors and dispersion")
  # Size factors and dispersion
  ssMo <- estimateSizeFactors(ssMo)
  ssMo <- estimateDispersions(ssMo)

  print("Monocle: Dispersed genes to use for pseudotime ordering")
  # Dispersed genes to use for pseudotime ordering
  disp_table <- dispersionTable(ssMo)
  ordering_genes <- subset(disp_table
    , mean_expression >= 0.01 & dispersion_empirical >= 0.75 * dispersion_fit)$gene_id
  print(length(expressed_genes))
  ssMo_filtered <- setOrderingFilter(ssMo, ordering_genes)
  # Plot
  # png(paste0(outGraph, "Monocle_OrderingGenesDispersion_Cluster.png"))
  #   # , paste0(clid, collapse = '-'), ".png"))
  # plot_ordering_genes(ssMo_filtered)
  # dev.off()

  # print("Variance explained by each PC")
  # # Variance explained by each PC
  # png(paste0(outGraph, "Monocle_PCA_VarianceExplained_Cluster"
  #   , title, ".png"))
  # plot_pc_variance_explained(ssMo_filtered, verbose = TRUE
  #   , use_existing_pc_variance = TRUE, return_all = FALSE)
  # dev.off()

  print("Monocle: Reduce data dimensionality")
  # Reduce data dimensionality
  # Use number of genes expressed or total mRNAs?
  ssMo_filtered <- reduceDimension(ssMo_filtered
    , max_components = max_components
    , residualModelFormulaStr = "~individual + librarylab + Total_mRNAs"
    , verbose = TRUE)

  print("Monocle: Order cells along trajectory")
  # Order cells along trajectory
  ssMo_filtered <- orderCells(ssMo_filtered)

  return(ssMo_filtered)

  ## Differential gene test as function of pseudotime
}

Format_Macosko_Cell_Cycle_Genes_Table <- function(ccDF) {
  print("Format_Macosko_Cell_Cycle_Genes_Table")
  # Subset to genes
  # Format Macosko cell cycle genes table
  genesGroupDF <- melt(ccDF, measure.vars = c("G1.S", "S", "G2.M", "M", "M.G1"))
  colnames(genesGroupDF) <- c("Grouping", "Gene.Symbol")
  genesGroupDF$Gene.Symbol <- gsub(" *", "", genesGroupDF$Gene.Symbol)
  genesGroupDF <- genesGroupDF[! genesGroupDF$Gene.Symbol == "", ]
  genesGroupDF <- genesGroupDF[! is.na(genesGroupDF$Grouping), ]
  genesGroupDF$Grouping <- gsub(" *", "", genesGroupDF$Grouping)
  genesGroupDF$Grouping <- factor(genesGroupDF$Grouping,
     levels = unique(genesGroupDF$Grouping))
  return(genesGroupDF)
}

Format_Genes_Table_To_Plot <- function(ccDF, kmDF, deDF, fold_changes) {

  print("Format_Genes_Table_To_Plot")

  # Faceting label for ggplot
  ccDF$Facet <- "Cell cycle"
  kmDF$Facet <- "Marker genes"
  deDF$Facet <- paste0("DE gene signatures ", deDF$Comparison[1])

  # Add known markers to cell cycle genes table
  genesGroupDF <- rbind(ccDF
    , kmDF[
        kmDF$Grouping %in% c("RG", "IP", "Neuron")
        , c("Grouping","Gene.Symbol","Facet")
      ]
    )

  # Add comparison column for later rbind to deDF
  genesGroupDF$Comparison <- NA

  # Add DE genes
  names(fold_changes) <- fold_changes
  dfl1 <- lapply(fold_changes, function(fold_change){

    print(paste("Fold change:", fold_change))

    df1 <- data.frame(Gene.Symbol = deDF$Gene[
      deDF$Log_FC_C2vC1 > fold_change & deDF$FDR < 0.05])
    df1$Grouping <- paste0("Log2 fold change > ", fold_change)
    df1$Comparison <- deDF$Comparison[1]
    df1$Facet <- deDF$Facet[1]

    df2 <- data.frame(Gene.Symbol = deDF$Gene[
      deDF$Log_FC_C2vC1 < -(fold_change) & deDF$FDR < 0.05])
    df2$Grouping <- paste0("Log2 fold change < -", fold_change)
    df2$Comparison <- deDF$Comparison[1]
    df2$Facet <- deDF$Facet[1]

    return(rbind(df1, df2))
  })
  df1 <- do.call("rbind", dfl1)
  genesGroupDF <- rbind(genesGroupDF, df1)

  return(genesGroupDF)
}

Format_GenesGroups_Pseudotime_For_GGplot <- function(
  genesGroupDF, ssMo_filtered){

  print("Format_GenesGroups_Pseudotime_For_GGplot")
  genes <- genesGroupDF$Gene.Symbol
  df1 <- melt(noCentExM[row.names(noCentExM) %in% genes
    , colnames(noCentExM) %in% pData(ssMo_filtered)$CELL])
  df1$Pseudotime <- pData(ssMo_filtered)$Pseudotime[
    match(df1$Var2, pData(ssMo_filtered)$CELL)]
  df1 <- merge(df1, genesGroupDF
    , by.x = "Var1", by.y = "Gene.Symbol", all = TRUE)
  df1 <- aggregate(value~Var2+Pseudotime+Grouping+Facet, df1, mean)
  df1$Cluster <- pData(ssMo_filtered)$res.0.6[
    match(df1$Var2, pData(ssMo_filtered)$CELL)]
  # df1$Grouping <- genesGroupDF$Grouping[
  #   match(df1$Var1, genesGroupDF$Gene.Symbol)]
  # df1$Facet <- genesGroupDF$Facet[
  #   match(df1$Var1, genesGroupDF$Gene.Symbol)]
  return(df1)
}

GGplot_GenesGroups_Pseudotime <- function(ggDF, title) {
  print("GGplot_GenesGroups_Pseudotime")
  ggL <- lapply(unique(ggDF$Facet), function(facet){
    ggDF <- ggDF[ggDF$Facet == facet, ]
    ggplot(ggDF, aes(x = Pseudotime, y = value, color = Grouping)) +
      geom_smooth() +
      ylab("Normalized expression") +
      ggtitle(facet)
  })
  # Collapse cells down to 1 row per cell
  gg2DF <- ggDF[! duplicated(ggDF$Var2), ]
  gg2 <- ggplot(gg2DF, aes(x = Pseudotime, y = Cluster, color = Cluster)) +
    geom_jitter()
  Plot_Grid(append(ggL, list(gg2)), ncol = 2, title = title, rel_height = 0.1
    , align = 'v', axis = 'l')
  ggsave(paste0(outGraph, "Monocle_MarkerExpr_Pseudotime_", title, ".png")
    , width = 11)
}

Monocle_Plot_Trajectory_By_Pseudotime <- function(mo, ...){
  print("Monocle_Plot_Trajectory_By_Pseudotime")
  gg <- plot_cell_trajectory(mo, color_by = "Pseudotime"
    , cell_size = 1, ...) +
    scale_color_viridis() +
    theme(legend.position = "right") +
    ggtitle("Pseudotime")
  return(gg)
}

Monocle_Plot_Trajectory_By_Seurat_Cluster <- function(mo, ...){
  print("Monocle_Plot_Trajectory_By_Seurat_Cluster")
  gg <- plot_cell_trajectory(mo, color_by = "res.0.6"
    , cell_size = 1, ...) +
    theme(legend.position = "right") +
    ggtitle("Seurat cluster")
  return(gg)
}

Monocle_Plot_Trajectory_By_State <- function(mo, title, ...) {
  print("Monocle_Plot_Trajectory_By_State")
  # Plot faceted by state
  plot_cell_trajectory(mo, color_by = "State", cell_size = 0.5, ...) +
    facet_wrap(~State, ncol = 4) +
    theme(legend.position = "none") +
    ggtitle(paste0(graphCodeTitle
      , "\n\nMonocle trajectory faceted by state"
      , "\n", title
      , "\n"))
  ggsave(paste0(outGraph, "Monocle_Trajectory_State_Comp12_", title, ".png"))
}

Monocle_Plot_Trajectory_By_Expression <- function(mo, genesGroupsDF, ...){
  print("Monocle_Plot_Trajectory_By_Expression")
  # Plot Luis markers - Seurat norm center scaled expression
  ggL <- lapply(unique(genesGroupsDF$Grouping), function(Group) {
    print(Group)
    genes <- genesGroupsDF$Gene.Symbol[genesGroupsDF$Grouping %in% Group]
    gg <- Plot_Trajectory_Gene_Expression(monocleO = mo
      , genes = genes
      # , exprM = as.matrix(exprs(mo_filtered))
      # , exprM = as.matrix(centSO@scale.data)
      , exprM = noCentExM
      , limHigh = 1, limLow = 0, title = Group, cell_size = 1, ...)
    gg <- gg + theme(text = element_text(size = 12))
    return(gg)
  })
  return(ggL)
}

Plot_Monocle_Analysis <- function(mo, ccDF, kmDF, deDF, fold_changes, title){
  print("Plot_Monocle_Analysis")
  # Gene expression by pseudotime
  ccDF <- Format_Macosko_Cell_Cycle_Genes_Table(ccDF)
  genesGroupDF <- Format_Genes_Table_To_Plot(
    ccDF
    , kmDF
    , deDF
    , fold_changes
  )
  ggDF <- Format_GenesGroups_Pseudotime_For_GGplot(
    genesGroupDF, ssMo_filtered = mo
  )
  print("Formatted ggDF:")
  print(head(ggDF))
  GGplot_GenesGroups_Pseudotime(ggDF, title)
  # Monocle trajectory
  # Component 1 and 2
  Monocle_Plot_Trajectory_By_State(mo, title)
  gg1 <- Monocle_Plot_Trajectory_By_Pseudotime(mo)
  gg2 <- Monocle_Plot_Trajectory_By_Seurat_Cluster(mo)
  ggL <- Monocle_Plot_Trajectory_By_Expression(
    mo = mo
    , genesGroupsDF = genesGroupDF
  )
  ggL <- append(list(gg1), ggL)
  ggL <- append(list(gg2), ggL)
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r'
    , title = paste0(title
      , "\n", genesGroupDF$Comparison[! is.na(genesGroupDF$Comparison)][1])
    )
  ggsave(paste0(outGraph, "Monocle_Trajectory_NormCentScale_", title, ".png"), width = 13, height = 16)
  # Component 1 and 3
  gg1 <- Monocle_Plot_Trajectory_By_Pseudotime(mo, x = 1, y = 3)
  gg2 <- Monocle_Plot_Trajectory_By_Seurat_Cluster(mo, x = 1, y = 3)
  ggL <- Monocle_Plot_Trajectory_By_Expression(
    mo = mo
    , genesGroupsDF = genesGroupDF, x = 1, y = 3
  )
  ggL <- append(list(gg1), ggL)
  ggL <- append(list(gg2), ggL)
  Plot_Grid(ggL, ncol = 3, align = 'v', axis = 'r'
    , title = paste0(title
      , "\n", genesGroupDF$Comparison[! is.na(genesGroupDF$Comparison)][1])
    )
  ggsave(paste0(outGraph, "Monocle_Trajectory_NormCentScale_", title, "_1v3.png"), width = 13, height = 16)
}
################################################################################

### Run

# Monocle for cluster 7, 8, 9, 10 cells
# Subset cells
clusters <- c(7,8,9,10)
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- names(centSO@ident[centSO@ident %in% clusters])
# Run monocle testing multiple max_components
lapply(c(6,8,10,20,40), function(max_components){
  ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
  Plot_Monocle_Analysis(
    mo = ssMO
    , ccDF = ccDF
    , kmDF = kmDF
    , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
    , fold_changes = c(0, 0.25, 0.5)
    , title = paste0("Cluster", paste0(clusters, collapse = "")
      , "_RG_IP_N_Comp", max_components
    )
  )
})

# Monocle for cluster 2, 7, 8, 9, 10 cells
# Subset cells
clusters <- c(2,7,8,9,10)
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- names(centSO@ident[centSO@ident %in% clusters])
# Run monocle testing multiple max_components
lapply(c(6,8,10,20,40), function(max_components){
  ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
  Plot_Monocle_Analysis(
    mo = ssMO
    , ccDF = ccDF
    , kmDF = kmDF
    , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
    , fold_changes = c(0, 0.25, 0.5)
    , title = paste0("Cluster", paste0(clusters, collapse = "")
      , "_RG_IP_N_Comp", max_components
    )
  )
})

# RG -> IP: Monocle for RG+ Neuron- cluster 7, 8, 9, 10 cells
# Subset cells
clusters <- c(7,8,9,10)
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
cellIDs <- intersect(cellIDs
  , names(centSO@ident[centSO@ident %in% clusters])
)
# Run monocle testing multiple max_components
lapply(c(6,8,10,20,40), function(max_components){
  ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
  Plot_Monocle_Analysis(
    mo = ssMO
    , ccDF = ccDF
    , kmDF = kmDF
    , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
    , fold_changes = c(0, 0.25, 0.5)
    , title = paste0("Cluster", paste0(clusters, collapse = "")
      , "_RGp_IP_Nn_Comp", max_components
    )
  )
})

# RG -> Neuron: Monocle for RG+ IP- cluster 7, 8, 9, 10 cells
# Subset cells
clusters <- c(7,8,9,10)
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
cellIDs <- intersect(cellIDs
  , names(centSO@ident[centSO@ident %in% clusters])
)
# Run monocle testing multiple max_components
lapply(c(6,8,10,20,40), function(max_components){
  ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
  Plot_Monocle_Analysis(
    mo = ssMO
    , ccDF = ccDF
    , kmDF = kmDF
    , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
    , fold_changes = c(0, 0.25, 0.5)
    , title = paste0("Cluster", paste0(clusters, collapse = "")
      , "_RGp_IPn_N_Comp", max_components
    )
  )
})

# IP -> Neuron: Monocle for IP+ RG- cluster 2, 8, 10 cells
# Subset cells
clusters <- c(2,8,10)
df1 <- Average_MarkersExp_Per_Cell(
  exM = noCentExM, seuratO = centSO)
cellIDs <- row.names(df1)[df1$IP > 0.5 & df1$RG < 0.25]
cellIDs <- intersect(cellIDs
  , names(centSO@ident[centSO@ident %in% clusters])
)
# Run monocle testing multiple max_components
lapply(c(6,8,10,20,40), function(max_components){
  ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
  Plot_Monocle_Analysis(
    mo = ssMO
    , ccDF = ccDF
    , kmDF = kmDF
    , deDF = deDF[deDF$Comparison == "IP_vs_Neuron", ]
    , fold_changes = c(0, 0.25, 0.5)
    , title = paste0("Cluster", paste0(clusters, collapse = "")
      , "_RGn_IPp_N_Comp", max_components
    )
  )
})

#
# # Monocle for RG+ IP- cluster 0, 8, 9, 10 cells
# # Subset cells
# clusters <- c(0,8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_N_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ IP- cluster 8, 9, 10 cells
# # Subset cells
# clusters <- c(8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_N_Comp", max_components
#     )
#   )
# })

# # Monocle for RG+ IP- cluster 0, 8, 10 cells
# # Subset cells
# clusters <- c(0,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8,10,15,20), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_N_Comp", max_components
#     )
#   )
# })

# # Monocle for RG+ Neuron- cluster 7, 8, 10 cells
# # Subset cells
# clusters <- c(7,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IP_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ Neuron- cluster 8, 9, 10 cells
# # Subset cells
# clusters <- c(8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IP_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ Neuron- cluster 7, 8, 9, 10 cells
# # Subset cells
# clusters <- c(7,8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IP_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+/Neuron+ IP- cluster 0, 8, 9, 10 cells
# # Subset cells
# clusters <- c(0,8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 |
#   df1$Neuron > 0.5 & df1$IP < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_Np_Comp", max_components
#     )
#   )
# })
#
# # Monocle for IP+/Neuron+ RG- cluster 0, 2, 8, 10 cells
# # Subset cells
# clusters <- c(0,2,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs <- row.names(df1)[df1$IP > 0.5 & df1$RG < 0.25 |
#   df1$Neuron > 0.5 & df1$RG < 0.25]
# cellIDs <- intersect(cellIDs
#   , names(centSO@ident[centSO@ident %in% clusters])
# )
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "IP_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGn_IPp_Np_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ IP+ Neuron- cluster 7, 8, 10 cells
# # Subset cells
# clusters <- c(7,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(7)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 & df1$IP > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPp_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ IP+ Neuron- cluster 8, 9, 10 cells
# # Subset cells
# clusters <- c(8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(9)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 & df1$IP > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPp_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ IP+ Neuron- cluster 7, 8, 9, 10 cells
# # Subset cells
# clusters <- c(7,8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(7,9)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$Neuron < 0.25 & df1$IP > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_IP", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPp_Nn_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ Neuron+ IP- cluster 7, 8, 10 cells
# # Subset cells
# clusters <- c(7,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(7)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_Np_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ Neuron+ IP- cluster 8, 9, 10 cells
# # Subset cells
# clusters <- c(8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(9)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_Np_Comp", max_components
#     )
#   )
# })
#
# # Monocle for RG+ Neuron+ IP- cluster 7, 8, 9, 10 cells
# # Subset cells
# clusters <- c(7,8,9,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(7,9)])
# )
# cellIDs2 <- row.names(df1)[df1$RG > 0.5 & df1$IP < 0.25 & df1$Neuron > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "RG_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGp_IPn_Np_Comp", max_components
#     )
#   )
# })
#
# # Monocle for IP+ Neuron+ RG- cluster 2, 8, 10 cells
# # Subset cells
# clusters <- c(2,8,10)
# df1 <- Average_MarkersExp_Per_Cell(
#   exM = noCentExM, seuratO = centSO)
# cellIDs1 <- row.names(df1)[df1$IP > 0.5 & df1$RG < 0.25]
# cellIDs1 <- intersect(cellIDs1
#   , names(centSO@ident[centSO@ident %in% c(2)])
# )
# cellIDs2 <- row.names(df1)[df1$IP > 0.5 & df1$RG < 0.25 & df1$Neuron > 0.5]
# cellIDs2 <- intersect(cellIDs2
#   , names(centSO@ident[centSO@ident %in% c(8,10)])
# )
# cellIDs <- c(cellIDs1, cellIDs2)
# # Run monocle testing multiple max_components
# lapply(c(4,6,8), function(max_components){
#   ssMO <- Monocle_Run(cellIDs = cellIDs, max_components = max_components)
#   Plot_Monocle_Analysis(
#     mo = ssMO
#     , ccDF = ccDF
#     , kmDF = kmDF
#     , deDF = deDF[deDF$Comparison == "IP_vs_Neuron", ]
#     , fold_changes = c(0, 0.25, 0.5)
#     , title = paste0("Cluster", paste0(clusters, collapse = "")
#       , "_RGn_IPp_Np_Comp", max_components
#     )
#   )
# })






# print("Differential gene test as function of pseudotime")
#
# ptDE <- differentialGeneTest(ssMo_filtered
#   , fullModelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
#   , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")


#
#
# diff_test_res <- ptDE
#
# # Subset to genes and by qval
# genes <- kmDF$Gene.Symbol[
#   kmDF$Grouping %in% c("RG", "IP", "Neuron", "vRG", "oRG")]
#
# diff_test_res <- diff_test_res[diff_test_res$gene_short_name %in% genes, ]
#
# genes <- row.names(diff_test_res)
#
# png(paste0(outGraph, "Monocle_DEpseudoTime_Heatmap.png")
#   , width = 9, height = 6, units = "in", res = 300)
# tryCatch(plot_pseudotime_heatmap(ssMo_filtered[genes, ]
#   , num_clusters = 3
#   , cores = 1
#   , show_rownames = T)
#   , error = function(cond) {
#     message("Error message:")
#     message(cond)
#     # Choose a return value in case of error
#     return(NULL)
#   }
# )
# dev.off()
################################################################################
