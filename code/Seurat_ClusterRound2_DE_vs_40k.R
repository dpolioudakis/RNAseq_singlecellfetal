# Damon Polioudakis
# 2018-11-20
# DE genes for each cluster

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(tidyverse)
require(reshape2)
source("Function_Library.R")

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
cluster_id <- (as.numeric(args[1]) - 1)
print(paste0("Cluster ID: ", cluster_id))

## Inputs
# Seurat
# Seurat cluster round 2
load("../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_Cluster3_seuratO.Robj")

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
ds_so <- centSO
rm(centSO)
rm(noCentExM)
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# ds_so <- ssCentSO

# Biomart to add ensembl IDs
bm_annot_df <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

## Variables
script_name <- "Seurat_ClusterRound2_DE_vs_40k.R"
date <- format(Sys.Date(), "%Y%m%d")
cluster_annot_tb <- tribble(
    ~cluster_number, ~cluster_annot
    , "9",  "vRG"
    , "7",  "oRG"
    , "8",  "PgS"
    , "10", "PgG2M"
    , "2",  "IP"
    , "0",  "ExN"
    , "1",  "ExM"
    , "4",  "ExCal"
    , "3",  "ExDp1"
    # , "13", "ExDp2"
    , "17", "SP"
    , "5",  "InSST"
    , "6",  "InCALB2"
    , "11", "OPC"
    , "12", "End"
    , "14", "Per"
    , "15", "Mic"
    , "16", "NA"
  )

# outputs
out_table <- paste0(
  "../analysis/tables/Seurat_ClusterRound2_DE_vs_40k/", date, "/")

# make output directories
dir.create(dirname(out_table), recursive = TRUE)
################################################################################

### Functions

################################################################################

### Differentially expressed genes for each cluster versus all other cells

de_cell_clusters <- function(
  expr_m
  , cluster_id
  , cell_cluster_key
  , design_matrix
  , mod
  , min_percent_expr = NULL
  , cluster_annot = NULL
  ){
    print("de_cell_clusters")

    print(paste0("cluster id: ", cluster_id))

    # check that all design matrix observations are found in expression matrix
    if(any(row.names(design_matrix) %in% colnames(expr_m)) == FALSE){
      stop("there are design_matrix rows not found in expr_m")
    }

    # Filter cells
    ss_expr_m <- DE_Filters_ExpMatrix(
      expr_m = expr_m, minPercent = min_percent_expr, clusterID = cluster_id
      , cell_cluster_key = cell_cluster_key
    )

    # DE Linear model
    de_lm <- DE_Linear_Model(
      exDatDF = ss_expr_m, termsDF = design_matrix, mod = mod)

    # Format LM output into data frame
    de_df <- Format_DE(
      deLM = de_lm, expr_m = expr_m, clusterID = cluster_id
      , cell_cluster_key = cell_cluster_key, cluster_annot = cluster_annot)

    # Add ensembl
    de_df$Ensembl <- Convert_Mixed_GeneSym_EnsID_To_EnsID(
        ids = as.character(de_df$Gene)
        , bmDF = bm_annot_df)

    # FDR correct
    de_df$FDR <- p.adjust(de_df$Pvalue, method = "BH")

    return(de_df)
}

make_cell_id_cluster_id_key <- function(){
  print("make_cell_id_cluster_id_key")
  # Gather cell IDs and cluster IDs (including SP subcluster) for heatmap
  cellid_clusterid <- ds_so@ident %>% enframe(
    name = "cell_id", value = "cluster_id") %>%
    left_join(
      (so@ident %>% enframe(name = "cell_id", value = "sub_cluster_id"))
      , by = "cell_id") %>%
    mutate_if(is.factor, as.character) %>%
    # Label subcluster 3:2 cells as 17 (subplate)
    # mutate(cluster_id = replace(cluster_id, sub_cluster_id %in% c(2), 17)) %>%
    # Label subcluster 3:2,0 cells as 17 (subplate)
    mutate(cluster_id = replace(cluster_id, sub_cluster_id %in% c(2,0), 17)) %>%
    # Label cluster 13 cells as 17 (subplate)
    mutate(cluster_id = replace(cluster_id, cluster_id == 13, 17)) %>%
    left_join(cluster_annot_tb
      , by = c("cluster_id" = "cluster_number")) %>%
      select(cell_id, cluster_id) %>% deframe
  return(cellid_clusterid)
}

make_design_matrix <- function(){
  print("make_design_matrix")
  design_matrix <- ds_so@meta.data[c("nUMI", "librarylab", "individual")]
  # Add term TRUE/FALSE cell is in cluster
  design_matrix$cluster <- FALSE
  design_matrix$cluster[make_cell_id_cluster_id_key() == cluster_id] <- TRUE
  return(design_matrix)
}

cluster_annot <- cluster_annot_tb %>%
  filter(cluster_number == cluster_id) %>%
  pull(cluster_annot)

de_df <- de_cell_clusters(
  expr_m = ds_so@data
  , cluster_id = cluster_id
  , cell_cluster_key = make_cell_id_cluster_id_key()
  , design_matrix = make_design_matrix()
  , mod = "y ~ cluster+nUMI+librarylab+individual"
  , min_percent_expr = 10
  , cluster_annot = cluster_annot
)

# Check
table(de_df$Pvalue < 0.05)
table(de_df$FDR < 0.05)

# Write to tab delimited table
print("Writing DE table to text file")
write.csv(x = de_df
  , file = paste0(out_table, "Cluster", cluster_id, "_Vs_All_Clusters.csv")
  , quote = FALSE, row.names = FALSE)
################################################################################
