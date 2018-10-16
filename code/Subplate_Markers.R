# Damon Polioudakis
# 2018-09-18
# Identify interesting subplate markers for follow up

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

options(stringsAsFactors = FALSE)

require(tidyverse)
require(reshape2)
require(cowplot)
require(viridis)
require(Seurat)
require(Matrix)
source("Function_Library.R")
source("GGplot_Theme.R")
################################################################################

## Inputs

# Seurat cluster round 2
load("../analysis/analyzed_data/Seurat_ClusterRound2/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/VarGenes/RegNumiLibBrain/PC1-40/Seurat_ClusterRound2_DS2-11_Cluster3_seuratO.Robj")

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
ds_so <- centSO
ds_ex_m <- noCentExM
rm(centSO)
rm(noCentExM)
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# ds_so <- ssCentSO
# ds_ex_m <- ssNoCentExM

# Subplate markers
sp_markers_tb <- read_csv(
  "../analysis/tables/Subplate/Subplate_Marker_Lists.csv")

# biomaRt gene info
bm_tb <- read_csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv")

# Allen Developmental Macaque human specific genes
hs_tb <- read_csv(
  "../source/Bakken_2016_AllenDevMacaque_ST10_HumanSpecific.csv"
  , skip = 1)

# Human gained enhancers
# Human gained enhancers gene list from Luis de la Torre-Ubieta.  Human gained
# enhancers assigned to genes by correlation of ATAC peaks on enhancer regions
# to ATAC peaks on promoters of genes.
# Human enhancer regions from Riley et al. ChIP-seq for enhancer marks comparing
# human, macaque, rodent at 3 different developmental time points.  Human gained
# have more pull down than other 2 species in ChIP-seq
hge_tb <- read_csv("../source/Gene_Lists/HGEs_Genes.csv")

# Disease risk gene datasets
asd_tada_tb <- load_asd_sanders_genes()
asd_ihart_tb <- load_asd_ihart_genes()
epilepsy_tb <- load_epilepsy_high_conf_genes()
id_genes_tb <- load_id_nejm_lancet_genes()

# Miller LCM from Allen
# Row annotations
miller_lcm_row_annot_df <- read.csv(
  "../allen_brain_data/miller_LCM/lmd_matrix_12840/rows_metadata.csv")
# Processed data from Luis
# [26] "layer.exp"
# [27] "layer.exp.sc"
# [28] "layer.max"
load("../allen_brain_data/miller_LCM/fetal_Ctx_layer_exp.Rdata")
idx <- match(rownames(layer.max), miller_lcm_row_annot_df$gene_symbol)
idx <- idx[! is.na(idx)]
miller_lcm_row_annot_df <- miller_lcm_row_annot_df[idx, ]

# TFs, chromatin remodelers, and co-factors
tf_cr_cf_tb <- load_tf_cofactors_remodelers()

# Known cell type markers
known_marker_tb <- read_csv(
  "../source/MarkersforSingleCell_2017-10-11_Markers.csv")


## Output Directories
out_graph <- "../analysis/graphs/Subplate/Subplate_Markers/Subplate_Markers_"
out_table <- "../analysis/tables/Subplate/Subplate_Markers/Subplate_Markers_"
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)

## Other variables

script_name <- "Subplate_Markers.R"
cluster_annot_key <- c(
  "9" = "vRG"
  , "7" = "oRG"
  , "8" = "PgS"
  , "10" = "PgG2M"
  , "2" = "IP"
  , "0" = "ExN"
  , "1" = "ExM"
  , "4" = "ExCal"
  , "3" = "ExDp1"
  # , "13" = "ExDp2"
  , "17" = "SP"
  , "5" = "InSST"
  , "6" = "InCALB2"
  , "11" = "OPC"
  , "12" = "End"
  , "14" = "Per"
  , "15" = "Mic"
  , "16" = "NA"
)
################################################################################

### Functions

plot_miller_lcm <- function(gene_group_df){
  # browser()
  print("plot_miller_lcm")
  # Allen LCM
  # By area and zone
  gg_dfl <- lapply(split(gene_group_df, gene_group_df$group)
    , function(ss_gene_group_df){
    genes <- ss_gene_group_df$gene
    print(genes)
    gg_df <- format_miller_lcm_by_region_for_ggplot(genes = genes)
    gg_df$gene_group <- ss_gene_group_df$group[1]
    return(gg_df)
  })
  # By zone
  gg_dfl <- lapply(gg_dfl, function(gg_df){
    aggregate(value~Zone+Var2+gene_group, gg_df, mean)
  })
  gg_df <- do.call("rbind", gg_dfl)
  # Plot
  ggplot(gg_df, aes(x = Var2, y = Zone, fill = value)) +
    # facet_wrap(~gene_group, space = "free", nrow = 1) +
    facet_grid(~gene_group, space = "free", scales = "free") +
    geom_tile() +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
        , palette = 5, direction = -1, limits = c(-1.5, 1.5)) +
    ylab("Region") +
    xlab("Gene") +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = element_text(angle = 90))
}

format_miller_lcm_by_region_for_ggplot <- function(genes){
  print("format_miller_lcm_by_region_for_ggplot")
  # browser()
  # Subset to genes of interest
  genes <- genes[! duplicated(genes)]
  ss_row_annot_DF <- miller_lcm_row_annot_df[
    miller_lcm_row_annot_df$gene_symbol %in% genes, ]
  idx <- match(rownames(layer.max), miller_lcm_row_annot_df$gene_symbol)
  miller_lcm_row_annot_df <- miller_lcm_row_annot_df[idx, ]
  ex_m <- layer.max[
    rownames(layer.max) %in% ss_row_annot_DF$gene_symbol, , drop = FALSE]
  # Select probe with highest expression for each gene
  ex_m <- t(ex_m)
  ex_m <- melt(ex_m)
  idx <- match(ex_m$Var2, miller_lcm_row_annot_df$gene_symbol)
  ex_m$gene_symbol <- miller_lcm_row_annot_df$gene_symbol[idx]
  ex_m$gene_symbol <- factor(ex_m$gene_symbol, levels = rev(genes))
  ex_m <- aggregate(value~gene_symbol+Var1, ex_m, max)
  ex_m <- ex_m[! duplicated(ex_m[ ,c("Var1", "gene_symbol")]), ]
  ex_m <- ex_m[! ex_m$Var1 == "", ]
  ex_m <- dcast(ex_m, Var1~gene_symbol, value.var = "value")
  rownames(ex_m) <- ex_m[ ,1]
  ex_m <- ex_m[ ,-1, drop = FALSE]
  # Center and scale
  ex_m <- scale(ex_m)
  ex_m <- melt(ex_m)
  # Set limits to -1.5 1.5
  ex_m$value[ex_m$value < -1.5] <- -1.5
  ex_m$value[ex_m$value > 1.5] <- 1.5
  # Add area
  ex_m$Area <- substr(ex_m$Var1, 1, 1)
  # Add zone
  ex_m$Zone <- substr(ex_m$Var1, 2, 3)
  ex_m$Zone <- factor(ex_m$Zone, levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ"))
  return(ex_m)
}

plot_tsne_expression <- function(genes){
  print("plot_tsne_expression")
  gg_l <- FeaturePlot(
    genes = genes
    , tsneDF = as.data.frame(ds_so@dr$tsne@cell.embeddings)
    , seuratO = ds_so
    , exM = ds_so@scale.data
    , limLow = -1.5
    , limHigh = 1.5
    , centScale = TRUE
  )
  return(gg_l)
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
    mutate(cluster_id = replace(cluster_id, sub_cluster_id %in% c(2), 17)) %>%
    # Label subcluster 3:2,0 cells as 17 (subplate)
    # mutate(cluster_id = replace(cluster_id, sub_cluster_id %in% c(2,0), 17)) %>%
    # Label cluster 13 cells as 17 (subplate)
    mutate(cluster_id = replace(cluster_id, cluster_id == 13, 17)) %>%
    left_join(
      (cluster_annot_key %>% enframe(
        name = "cluster_id", value = "cluster_annot"))
      , by = "cluster_id") %>%
      select(cell_id, cluster_annot) %>% deframe
  return(cellid_clusterid)
}

plot_expression_by_cluster_boxplot <- function(
  genes
  , expr_m
  , cell_cluster_key
  , cluster_order = NULL
  ){
  print("plot_expression_by_cluster_boxplot")

  # make cell_cluster_key a tibble if it's not
  if (any(class(cell_cluster_key) != "tibble")) {
    cell_cluster_key_tb <- cell_cluster_key %>%
      enframe(name = "cell_id", value = "cluster") %>%
      as_tibble
  }

  # gather data
  gg_tb <- expr_m %>%
    as_tibble(rownames = "gene") %>%
    filter(gene %in% genes) %>%
    gather(., key = "cell_id", value = "expression", -gene) %>%
    left_join(., cell_cluster_key_tb) %>%
    # order clusters
    mutate(cluster = fct_relevel(cluster, cluster_annot_key)) %>%
    # order genes
    mutate(genes = fct_relevel(gene, unique(genes)))

  # plot
  ggplot(gg_tb, aes(x = cluster, y = expression, fill = cluster)) +
    facet_wrap(~gene, ncol = 4, scales = "free") +
    geom_boxplot(outlier.size = 0.1) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      , legend.position = "none")
}

determine_number_of_cells_genes_are_both_expressed_in <- function(
  genes, expr_m, title){

  print("determine_number_of_cells_genes_are_both_expressed_in")

  m1 <- expr_m[row.names(expr_m) %in% genes, ] > 0.5
  l1 <- apply(m1, 2, function(col) row.names(m1)[col])

  ldf <- lapply(l1, function(genes) {
    v2 <- genes[genes %in% genes]
    v3 <- genes[genes %in% genes]
    expand.grid(v2, v3)
  })
  df1 <- do.call("rbind", ldf)

  # format for matrix
  df3 <- dcast(df1, Var1 ~ Var2)
  row.names(df3) <- df3$Var1
  df3 <- df3[ ,-1]

  # long format
  df3$gene1 <- row.names(df3)
  df3 <- melt(df3)
  df3$gene1 <- factor(df3$gene, levels = genes)
  df3$gene2 <- factor(df3$variable, levels = genes)
  df3 <- df3[ ,c("gene1", "gene2", "value")]

  return(df3)
}

calculate_jaccard_overlap_of_cells_expressing_two_genes_of_interest <-
  function(gene1, gene2){

  print("calculate_jaccard_overlap_of_cells_expressing_genes")
  print(gene1, gene2)

  # tibble of gene1, gene2, cell IDs, and if gene is expressed in a given cell
  expressed_tb <- ds_ex_m %>% as_tibble(rownames = "gene") %>% filter(gene %in% c(gene1, gene2)) %>% gather(., key = "cell_id", value = "expression", -gene) %>%
    mutate(expressed = expression > 0.5)

  # select cell ids of cells that express gene 1 or 2
  cell_ids_1 <- filter(expressed_tb, gene == gene1 & expressed == TRUE) %>%
    pull(cell_id)
  cell_ids_2 <- filter(expressed_tb, gene == gene2 & expressed == TRUE) %>%
    pull(cell_id)

  # calculate jaccard index using cell ids
  jaccard_index <- length(intersect(cell_ids_1, cell_ids_2)) /
    length(union(cell_ids_1, cell_ids_2))

  return(jaccard_index)
}
################################################################################

### Clean data

# Clean variable names
sp_markers_tb <- clean_variable_names(sp_markers_tb)
bm_tb <- clean_variable_names(bm_tb)
hs_tb <- clean_variable_names(hs_tb)
hge_tb <- clean_variable_names(hge_tb)
asd_tada_tb <- clean_variable_names(asd_tada_tb)
asd_ihart_tb <- clean_variable_names(asd_ihart_tb)
epilepsy_tb <- clean_variable_names(epilepsy_tb)
id_genes_tb <- clean_variable_names(id_genes_tb)
tf_cr_cf_tb <- clean_variable_names(tf_cr_cf_tb)
known_marker_tb <- clean_variable_names(known_marker_tb)

# Subset marker sets
sp_markers_tb <- sp_markers_tb %>% filter(group %in% c(
  "SP_vs_all_DE04", "SP_vs_all_DE04_top80", "SP_vs_cluster3_DE04"
  , "ST18_correlate_top80", "Intersect_SPDE_ST18cor"))
################################################################################

### Heatmaps of expression in Miller LCM of SP markers intersect with gene of
#interest

plot_miller_heatmaps <- function(){
  print("plot_miller_heatmaps")

  # Are SP markers human or primate specific?
  # Allen Developmental Macaque human specific genes
  sp_markers_tb %>%
    left_join(hs_tb, by = "gene") %>%
    filter(set %in% c("Human-specific", "Primate-specific")) %>%
    plot_miller_lcm() +
    ggtitle(paste0(script_name
      , "\nIntersection of subplate markers and human specific expression pattern genes"))
  ggsave(paste0(out_graph, "human_specific_miller_heatmap.pdf")
    , width = 9, height = 5)

  # Do SP markers have human gained enhancers?
  gg_1 <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "GZ>CP") %>%
    plot_miller_lcm() +
      ggtitle("GZ>CP")
  gg_2 <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "CP>GZ") %>%
    plot_miller_lcm() +
      ggtitle("CP>GZ")
  Plot_Grid(list(gg_1, gg_2), ncol = 1
      , title = paste0(script_name, "\nIntersection of subplate markers and human gained enhancers"))
  ggsave(paste0(out_graph, "hge_miller_heatmap.pdf")
    , width = 9, height = 10)

  # Are SP markers disease risk genes?
  # ASD Sanders
  sp_markers_tb %>%
    inner_join(asd_tada_tb, by = "gene") %>%
    plot_miller_lcm() +
      ggtitle(paste0(script_name
        , "\nIntersection of subplate markers and ASD (Sanders)"))
  ggsave(paste0(out_graph, "asd_sanders_miller_heatmap.pdf")
    , width = 9, height = 5)
  # ASD ihart
  sp_markers_tb %>%
    inner_join(asd_ihart_tb, by = c("gene" = "hgnc_gene_symbol")) %>%
    plot_miller_lcm() +
      ggtitle(paste0(script_name
        , "\nIntersection of subplate markers and ASD (ihart)"))
  ggsave(paste0(out_graph, "asd_ihart_miller_heatmap.pdf")
    , width = 9, height = 5)
  # Epilepsy
  sp_markers_tb %>%
    inner_join(epilepsy_tb, by = "gene") %>%
    plot_miller_lcm() +
      ggtitle(paste0(script_name
        , "\nIntersection of subplate markers and epilepsy (Ruzzo)"))
  ggsave(paste0(out_graph, "epilepsy_miller_heatmap.pdf")
    , width = 9, height = 5)
  # ID
  sp_markers_tb %>%
    inner_join(id_genes_tb, by = "gene") %>%
    plot_miller_lcm() +
      ggtitle(paste0(script_name
        , "\nIntersection of subplate markers and ID (deLigt and Rauch)"))
  ggsave(paste0(out_graph, "id_miller_heatmap.pdf")
    , width = 9, height = 5)

    # TFs
    sp_markers_tb %>%
      inner_join(tf_cr_cf_tb, by = c("gene" = "hgnc_symbol")) %>%
      plot_miller_lcm() +
      ggtitle(paste0(script_name
        , "\nIntersection of subplate markers and TFs, co-factors, and chromatin remodelers"))
    ggsave(paste0(out_graph, "tfs_miller_heatmap.pdf")
      , width = 9, height = 5)
}
################################################################################

### tSNE expression of subplate markers intersected with genes of interest

plot_dropseq_expression_tsnes <- function(){
  print("plot_dropseq_expression_tsnes")

  sp_markers_tb %>%
    left_join(hs_tb, by = "gene") %>%
    filter(set %in% c("Human-specific", "Primate-specific")) %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and human specific expression pattern genes"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "human_specific_featureplotind_normcentscale.png")
      , width = 13, height = 18)

  # Do SP markers have human gained enhancers?
  sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and human gained enhancer linked genes"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "hge_featureplotind_normcentscale.png")
      , width = 13, height = 18)

  # Are SP markers disease risk genes?
  # ASD Sanders
  sp_markers_tb %>%
    inner_join(asd_tada_tb, by = "gene") %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and ASD (Sanders)"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "asd_sanders_featureplotind_normcentscale.png")
      , width = 13, height = 15)

  # ASD ihart
  sp_markers_tb %>%
    inner_join(asd_ihart_tb, by = c("gene" = "hgnc_gene_symbol")) %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and ASD (ihart)"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "asd_ihart_featureplotind_normcentscale.png")
      , width = 13, height = 13)

  # Epilepsy
  sp_markers_tb %>%
    inner_join(epilepsy_tb, by = "gene") %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and epilepsy (Ruzzo)"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "epilepsy_featureplotind_normcentscale.png")
      , width = 13, height = 18)

  # ID
  sp_markers_tb %>%
    inner_join(id_genes_tb, by = "gene") %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and ID (deLigt and Rauch)"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "id_featureplotind_normcentscale.png")
      , width = 13, height = 9)

  sp_markers_tb %>%
    inner_join(tf_cr_cf_tb, by = c("gene" = "hgnc_symbol")) %>%
    pull(gene) %>%
    # Feature plot normalized, mean centered scaled
    # Individual expression
    plot_tsne_expression(.) %>%
    Plot_Grid(., ncol = 4, align = 'v', axis = 'r', rel_height = c(0.1)
      , title = paste0(script_name
        , "\n\nIntersection of subplate markers and TFs, co-factors, and chromatin remodelers"
        , "\nNormalized expression, mean centered and variance scaled"
        , "\n")
    )
    ggsave(paste0(out_graph, "tfs_featureplotind_normcentscale.png")
      , width = 13, height = 21)
}
################################################################################

### Heatmaps expression of subplate markers intersected with genes of interest

plot_dropseq_expression_heatmaps <- function(){
  print("plot_dropseq_expression_heatmaps")

  ## Allen human specific
  gene_group_df <- sp_markers_tb %>%
    left_join(hs_tb, by = "gene") %>%
    filter(set %in% c("Human-specific")) %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and human specific expression pattern genes"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "human_specific_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 9
  )

  ## Allen primate specific
  gene_group_df <- sp_markers_tb %>%
    left_join(hs_tb, by = "gene") %>%
    filter(set %in% c("Primate-specific")) %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and primate specific expression pattern genes"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "primate_specific_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 9
  )

  ## HGE
  # GZ>CP
  gene_group_df <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "GZ>CP") %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  gg1 <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
  ggtitle("GZ>CP")
  # CP>GZ
  gene_group_df <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "CP>GZ") %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  gg2 <- Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
  ggtitle("CP>GZ")
  Plot_Grid(list(gg1, gg2), ncol = 1
    , title = paste0(script_name
      , "\n\nIntersection of subplate markers and human gained enhancer linked genes"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "hge_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 15
  )

  ## ASD Sanders
  gene_group_df <- sp_markers_tb %>%
    inner_join(asd_tada_tb, by = "gene") %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and ASD (Sanders)"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "asd_sanders_specific_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 10
  )

  ## ASD ihart
  gene_group_df <- sp_markers_tb %>%
    inner_join(asd_ihart_tb, by = c("gene" = "hgnc_gene_symbol")) %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and ASD (ihart)"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "asd_ihart_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 10
  )

  ## Epilepsy
  gene_group_df <- sp_markers_tb %>%
    inner_join(epilepsy_tb, by = "gene") %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and epilepsy (Ruzzo)"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "epilepsy_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 10
  )

  ## ID
  gene_group_df <- sp_markers_tb %>%
    inner_join(id_genes_tb, by = "gene") %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  # Plot
  # Normalized centered scaled
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and ID (deLigt and Rauch)"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "id_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 6
  )

  ## TFs
  gene_group_df <- sp_markers_tb %>%
    inner_join(tf_cr_cf_tb, by = c("gene" = "hgnc_symbol")) %>%
    select(Gene = gene, Group = group) %>% as.data.frame
  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = ds_so@scale.data
    , cellID_clusterID <- make_cell_id_cluster_id_key()
    , clusters = cluster_annot_key[cluster_annot_key != "NA"]
    , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
  ) +
    ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and TFs, co-factors, and chromatin remodelers"
      , "\nx-axis: Genes"
      , "\ny-axis: Cells ordered by cluster"
      , "\nNormalized expression, mean centered, variance scaled"
      , "\n")
    )
  ggsave(paste0(
      out_graph
      , "tfs_heatmapsetcolwid_normcentscale.png")
    , width = 10, height = 10
  )
}
################################################################################

### Drop-seq expression boxplots

plot_dropseq_expression_boxplots <- function(){

  print("plot_dropseq_expression_boxplots")

  # Allen Developmental Macaque human specific genes
  sp_markers_tb %>%
    left_join(hs_tb, by = "gene") %>%
    filter(set %in% c("Human-specific")) %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
        , "\n\nIntersection of subplate markers and ASD (Sanders)"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "asd_sanders_miller_heatmap.pdf")
    , width = 11, height = 5)

  # Do SP markers have human gained enhancers?
  gg_1 <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "GZ>CP") %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle("GZ>CP")
  gg_2 <- sp_markers_tb %>%
    inner_join(hge_tb, by = c("gene" = "hgnc")) %>%
    filter(gz_or_cp == "CP>GZ") %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
    ggtitle("CP>GZ")
  Plot_Grid(list(gg_1, gg_2), ncol = 1
      , title = paste0(
        script_name
        , "\nIntersection of subplate markers and human gained enhancers"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "hge_miller_heatmap.pdf")
    , width = 11, height = 5)

  # ASD Sanders
  sp_markers_tb %>%
    inner_join(asd_tada_tb, by = "gene") %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
        , "\n\nIntersection of subplate markers and ASD (Sanders)"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "asd_sanders_miller_heatmap.pdf")
    , width = 11, height = 5)

  # ASD ihart
  sp_markers_tb %>%
    inner_join(asd_ihart_tb, by = c("gene" = "hgnc_gene_symbol")) %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
        , "\n\nIntersection of subplate markers and ASD (ihart)"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "asd_ihart_miller_heatmap.pdf")
    , width = 11, height = 5)

  # Epilepsy
  sp_markers_tb %>%
    inner_join(epilepsy_tb, by = "gene") %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
        , "\n\nIntersection of subplate markers and epilepsy (Ruzzo)"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "epilepsy_miller_heatmap.pdf")
    , width = 11, height = 5)

  # ID
  sp_markers_tb %>%
    inner_join(id_genes_tb, by = "gene") %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
        , "\n\nIntersection of subplate markers and ID (deLigt and Rauch)"
        , "\nNormalized expression"))
  ggsave(paste0(out_graph, "id_miller_heatmap.pdf")
    , width = 11, height = 5)

  # TFs, chromatin remodelers, co-factors
  sp_markers_tb %>%
    inner_join(tf_cr_cf_tb, by = c("gene" = "hgnc_symbol")) %>%
    pull(gene) %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
      , "\n\nIntersection of subplate markers and TFs, co-factors, and chromatin remodelers"
      , "\nNormalized expression"))
  ggsave(paste0(out_graph, "tfs_boxplot_normcentscale.png")
    , width = 11, height = 11)

}
################################################################################

### Cells expressing combinations of subplate and / or deep layer markers

# How much do subplate markers overlap with deep layer markers?

# Are subplate markers in the same cell?

plot_cells_expressing_gene_combinations_heatmaps <- function(){

  print("plot_cells_expressing_gene_combinations_heatmaps")

  # Numbers of cells both genes are expressed in
  genes <- c(sp_markers_tb %>%
      filter(group == "Intersect_SPDE_ST18cor") %>% pull(gene)
    , known_marker_tb %>% filter(grouping == "Excitatory Deep Layer Cortical") %>% pull(gene_symbol)
  )
  genes %>% determine_number_of_cells_genes_are_both_expressed_in(
      genes = ., expr_m = ds_ex_m) %>%
    as_tibble() %>%
    ggplot(aes(x = gene1, y = gene2, fill = value)) +
      geom_tile() +
      geom_text(aes(x = gene1, y = gene2, label = value), color = "black") +
      scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "Number of cells") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Genes") +
      ylab("Genes") +
      ggtitle(paste0(script_name
        , "\n\nNumber of cells both genes are expressed in (>0.5 normalized expression)"))
      ggsave(paste0(out_graph, "number_cells_both_genes_expressed_heatmap.png")
        , width = 11, height = 11)


  # jaccard index heatmap of cells expressing both genes

  genes <- c(sp_markers_tb %>%
      filter(group == "Intersect_SPDE_ST18cor") %>% pull(gene)
    , known_marker_tb %>% filter(grouping == "Excitatory Deep Layer Cortical") %>% pull(gene_symbol)
  )

  expand.grid(genes, genes) %>%
    as_tibble %>%
    mutate_if(is.factor, as.character) %>%
    rename(gene1 = Var1, gene2 = Var2) %>%
    mutate(jaccard_index = pmap(.
      , calculate_jaccard_overlap_of_cells_expressing_two_genes_of_interest) %>%
        unlist %>% round(., 2)) %>%
    mutate(gene1 = factor(gene1, levels = unique(gene1))) %>%
    mutate(gene2 = factor(gene2, levels = unique(gene1))) %>%
    ggplot(aes(x = gene1, y = gene2, fill = jaccard_index)) +
      geom_tile() +
      geom_text(aes(x = gene1, y = gene2, label = jaccard_index), color = "black") +
      scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "Number of cells") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Genes") +
      ylab("Genes") +
      ggtitle(paste0(script_name
        , "\n\nJaccard index of cells both genes are expressed in (>0.5 normalized expression)"))
      ggsave(paste0(out_graph
          , "number_cells_both_genes_expressed_jaccard_heatmap.png")
        , width = 14, height = 14)
}
################################################################################

### Heatmap of deep layer markers

plot_deep_layer_marker_expression <- function(){

  print("plot_deep_layer_marker_expression")

  known_marker_tb %>%
    filter(grouping == "Excitatory Deep Layer Cortical") %>%
    select(gene_symbol, grouping) %>%
    rename(Gene = gene_symbol, Group = grouping) %>% as.data.frame %>%
    Plot_Marker_Genes_Heatmap_SetColWidths(
      geneGroupDF = .
      , exprM = ds_so@scale.data
      , cellID_clusterID <- make_cell_id_cluster_id_key()
      , clusters = cluster_annot_key[cluster_annot_key != "NA"]
      , clusterOrder = cluster_annot_key[cluster_annot_key != "NA"]
    ) +
      ggtitle(paste0(script_name
        , "\n\nDeep layer markers"
        , "\nx-axis: Genes"
        , "\ny-axis: Cells ordered by cluster"
        , "\nNormalized expression, mean centered, variance scaled"
        , "\n")
      )
    ggsave(paste0(
        out_graph
        , "deeplayermarkers_heatmapsetcolwid_normcentscale.png")
      , width = 10, height = 5
    )

  # TFs, chromatin remodelers, co-factors
  known_marker_tb %>%
    filter(grouping == "Excitatory Deep Layer Cortical") %>%
    select(gene_symbol, grouping) %>%
    pull(gene_symbol) %>%
    plot_expression_by_cluster_boxplot(
      genes = .
      , expr_m = ds_ex_m
      , cell_cluster_key = make_cell_id_cluster_id_key()
      , cluster_order = cluster_annot_key
    ) +
      ggtitle(paste0(script_name
      , "\n\nDeep layer markers"
      , "\nNormalized expression"))
  ggsave(paste0(out_graph, "deeplayermarkers_boxplot_normcentscale.png")
    , width = 11, height = 11)
}
################################################################################

### Run

main_function <- function(){
  plot_miller_heatmaps()
  plot_dropseq_expression_tsnes()
  plot_dropseq_expression_heatmaps()
  plot_dropseq_expression_boxplots()
  plot_cells_expressing_gene_combinations_heatmaps()
  plot_deep_layer_marker_expression()
}
main_function()
################################################################################
