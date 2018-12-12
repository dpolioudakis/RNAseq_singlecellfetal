# Damon Polioudakis
# 2018-11-14
# DE genes for each cluster

# must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(Seurat)
require(reshape2)
require(tidyverse)
require(cowplot)
require(viridis)
source("Function_Library.R")
source("GGplot_Theme.R")

## command args to input cluster ID
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# cluster_id <- (as.numeric(args[1]) - 1)
# print(paste0("Cluster ID: ", cluster_id))

## other variables
script_name <- "Seurat_ClusterRound2_DE.R"
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
    , "13", "ExDp2"
    , "5",  "InSST"
    , "6",  "InCALB2"
    , "11", "OPC"
    , "12", "End"
    , "14", "Per"
    , "15", "Mic"
  )

## inputs
# sub-clustering DE enrichment tables
# (DE sub-cluster vs other sub-clusters in the cluster)
in_subclust_de_dir <- "../analysis/tables/Seurat_ClusterRound2/DE/20180910/"
# sub-clustering DE enrichment vs full dataset tables
# (DE sub-cluster vs all other cells in dataset)
subclust_de_vs_40k_tb <- read_csv("../analysis/tables/Seurat_ClusterRound2/DE_vs_40k/20181210/Seurat_ClusterRound2_DE_vs_40k_ClusterX_Vs_All_Clusters.csv")
# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST5000_seuratO.Robj")
# centSO <- ssCentSO
# noCentExM <- ssNoCentExM
# cluster DE enrichment table
# Cluster DE table
cluster_de_tb <- read_tsv(
  "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt")

## outputs
out_graph <- paste0(
  "../analysis/graphs/Seurat_ClusterRound2/DE/", date, "/Seurat_ClusterRound2_")
################################################################################

### main function

main_function <- function(){
  run_plot_de_vs_40k_expression_heatmaps()
  run_plot_de_vs_subcluster_expression_heatmaps()
  run_plot_subcluster_expression_heatmaps()
  run_plot_enriched_gene_overlap_odds_ratio_heatmaps()
}
################################################################################

### inputs / outputs

# DE for sub-clusters
subclust_de_tb <- in_subclust_de_dir %>%
  # list.files(pattern = paste0("_Cluster", cluster_id), full.names = TRUE) %>%
  list.files(pattern = paste0("_Cluster"), full.names = TRUE) %>%
  map_df(~read_csv(.))

# make output directories
dir.create(dirname(out_graph), recursive = TRUE)
################################################################################

### differentially expressed genes for each cluster vs other cells in cluster

plot_subcluster_expression_heatmaps <- function(subclust_de_tb, cluster_id){

  print("plot_subcluster_expression_heatmaps")
  print(paste0("cluster: ", cluster_id))

  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", cluster_id, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)

  gene_group_df <-
    subclust_de_tb %>%
    filter(Cluster == cluster_id & FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
    group_by(Subcluster) %>%
    top_n(n = 30, wt = Log2_Fold_Change) %>%
    select(Gene = Gene, Group = Subcluster) %>%
    as.data.frame()

  expr_m <- as.matrix(so@scale.data)

  cluster_annot <- cluster_annot_tb %>%
    filter(cluster_number == cluster_id) %>%
    pull(cluster_annot)

  cluster_order <- so@ident %>% unique %>% sort

  graph_height = 2 + 4*(length(so@ident %>% unique))

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = expr_m
    , cellID_clusterID = so@ident
    , clusterOrder = cluster_order
    , lowerLimit = -1.5
    , upperLimit = 1.5
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered variance scaled by gene"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id, "_subclust_expression_heatmap_zscore.png")
    , height = graph_height, width = 9)

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = noCentExM
    , cellID_clusterID = so@ident
    , clusterOrder = cluster_order
    , lowerLimit = 0
    , upperLimit = 4
    , color_viridis = TRUE
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id, "_subclust_expression_heatmap.png")
      , height = graph_height, width = 9)
}

run_plot_subcluster_expression_heatmaps <- function(){

  print("run_plot_subcluster_expression_heatmaps")

  cluster_annot_tb$cluster_number %>%
    map(
      ~tryCatch({
        plot_subcluster_expression_heatmaps(
        subclust_de_tb = subclust_de_tb
        # , so = so
        , cluster_id = .
          )
      }
        , error = function(cond){
          message(paste0("Error for cluster: ", .))
          message(paste0(cond, "\n"))
          return(NA)
      }
    ))
}
################################################################################

### differentially expressed genes for each sub-cluster vs all other cells

plot_expression_de_vs_40k_heatmaps <- function(
  cluster_id, seurat_obj, noCent40kExM){

  print("plot_expression_de_vs_40k_heatmaps")
  print(paste0("cluster: ", cluster_id))

  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", cluster_id, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)

  # Gather cell IDs and cluster IDs (including SP subcluster) for heatmap
  cellid_clusterid <- seurat_obj@ident %>% enframe(
    name = "cell_id", value = "cluster_ids") %>%
    left_join(
      (so@ident %>% enframe(name = "cell_id", value = "sub_cluster_id"))
      , by = "cell_id") %>%
    mutate_if(is.factor, as.character) %>%
    mutate(cluster_ids = if_else(cluster_ids == cluster_id, paste0(.$cluster_ids, "_", .$sub_cluster_id), cluster_ids))  %>%
    left_join(cluster_annot_tb
      , by = c("cluster_ids" = "cluster_number")) %>%
      select(cell_id, cluster_ids) %>% deframe

  # # fixes column formatting... delete this after rerunning DE script...
  # subclust_de_vs_40k_tb <- subclust_de_vs_40k_tb %>%
  #   mutate(Subcluster = gsub("*._", "", Cluster)) %>%
  #   mutate(Cluster_Subcluster = Cluster) %>%
  #   mutate(Cluster = gsub("_.*", "", Cluster))

  gene_group_df <-
    subclust_de_vs_40k_tb %>%
    filter(Cluster == cluster_id & FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
    group_by(Subcluster) %>%
    top_n(n = 30, wt = Log2_Fold_Change) %>%
    arrange(Subcluster, Log2_Fold_Change) %>%
    select(Gene = Gene, Group = Subcluster) %>%
    as.data.frame()

  expr_m <- as.matrix(seurat_obj@scale.data)

  cluster_annot <- cluster_annot_tb %>%
    filter(cluster_number == cluster_id) %>%
    pull(cluster_annot)

  cluster_order <- cluster_annot_tb %>%
    left_join(.
      , subclust_de_vs_40k_tb %>%
          filter(Cluster == cluster_id) %>%
          select(Cluster, Cluster_Subcluster, Subcluster) %>%
          distinct %>%
          arrange(Subcluster) %>%
          mutate(Cluster = as.character(Cluster))
      , by = c("cluster_number" = "Cluster")) %>%
    mutate(cluster_number = if_else(
      cluster_number == cluster_id, Cluster_Subcluster, cluster_number)) %>%
    pull(cluster_number)

  graph_height = 2 + 4*(length(so@ident %>% unique))

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = expr_m
    , cellID_clusterID = cellid_clusterid
    , clusters = cellid_clusterid %>% unique
    , clusterOrder = cluster_order
    , lowerLimit = -1.5
    , upperLimit = 1.5
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nEnrichment determined by DE of cells in sub-cluster versus other in dataset"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered variance scaled by gene"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id, "_de_vs_40k_expression_heatmap_zscore.png")
    , height = graph_height, width = 9)

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = noCent40kExM
    , cellID_clusterID = cellid_clusterid
    , clusters = cellid_clusterid %>% unique
    , clusterOrder = cluster_order
    , lowerLimit = 0
    , upperLimit = 4
    , color_viridis = TRUE
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nEnrichment determined by DE of cells in sub-cluster versus other in dataset"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id, "_de_vs_40k_expression_heatmap.png")
      , height = graph_height, width = 9)
}

run_plot_de_vs_40k_expression_heatmaps <- function(){

  print("run_plot_de_vs_40k_expression_heatmaps")

  cluster_annot_tb$cluster_number %>%
    map(
      ~tryCatch({
          plot_expression_de_vs_40k_heatmaps(
            seurat_obj = centSO
            , noCent40kExM = noCentExM
            # , so = so
            , cluster_id = .
            )
        }
          , error = function(cond){
            message(paste0("Error for cluster: ", .))
            message(paste0(cond, "\n"))
            return(NA)
        }
      ))
}
################################################################################

### differentially expressed genes for each sub-cluster vs all other cells

plot_de_vs_subclust_expression_heatmaps <- function(
  cluster_id, seurat_obj, noCent40kExM){

  print("plot_de_vs_subclust_expression_heatmaps")
  print(paste0("cluster: ", cluster_id))

  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", cluster_id, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)

  # Gather cell IDs and cluster IDs (including SP subcluster) for heatmap
  cellid_clusterid <- seurat_obj@ident %>% enframe(
    name = "cell_id", value = "cluster_ids") %>%
    left_join(
      (so@ident %>% enframe(name = "cell_id", value = "sub_cluster_id"))
      , by = "cell_id") %>%
    mutate_if(is.factor, as.character) %>%
    mutate(cluster_ids = if_else(cluster_ids == cluster_id, paste0(.$cluster_ids, "_", .$sub_cluster_id), cluster_ids))  %>%
    left_join(cluster_annot_tb
      , by = c("cluster_ids" = "cluster_number")) %>%
      select(cell_id, cluster_ids) %>% deframe

  # standardize column formatting
  subclust_de_tb <- subclust_de_tb %>%
    mutate(Cluster_Subcluster = paste0(Cluster, "_", Subcluster)) %>%
    mutate(Cluster = as.character(Cluster))

  gene_group_df <-
    subclust_de_tb %>%
    filter(Cluster == cluster_id & FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
    group_by(Subcluster) %>%
    top_n(n = 30, wt = Log2_Fold_Change) %>%
    arrange(Subcluster, Log2_Fold_Change) %>%
    select(Gene = Gene, Group = Subcluster) %>%
    as.data.frame()

  expr_m <- as.matrix(seurat_obj@scale.data)

  cluster_annot <- cluster_annot_tb %>%
    filter(cluster_number == cluster_id) %>%
    pull(cluster_annot)

  cluster_order <- cluster_annot_tb %>%
    left_join(.
      , subclust_de_tb %>%
          filter(Cluster == cluster_id) %>%
          select(Cluster, Cluster_Subcluster, Subcluster) %>%
          distinct %>%
          arrange(Subcluster) %>%
          mutate(Cluster = as.character(Cluster))
      , by = c("cluster_number" = "Cluster")) %>%
    mutate(cluster_number = if_else(
      cluster_number == cluster_id, Cluster_Subcluster, cluster_number)) %>%
    pull(cluster_number)

  graph_height = 2 + 4*(length(so@ident %>% unique))

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = expr_m
    , cellID_clusterID = cellid_clusterid
    , clusters = cellid_clusterid %>% unique
    , clusterOrder = cluster_order
    , lowerLimit = -1.5
    , upperLimit = 1.5
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nEnrichment determined by DE of cells in sub-cluster versus other in cluster"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression, mean centered variance scaled by gene"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id
    , "_de_vs_subclust_expression_heatmap_zscore.png")
    , height = graph_height, width = 9)

  Plot_Marker_Genes_Heatmap_SetColWidths(
    geneGroupDF = gene_group_df
    , exprM = noCent40kExM
    , cellID_clusterID = cellid_clusterid
    , clusters = cellid_clusterid %>% unique
    , clusterOrder = cluster_order
    , lowerLimit = 0
    , upperLimit = 4
    , color_viridis = TRUE
  ) + ggtitle(paste0(
    script_name
    , "\n\nExpression of top 30 enriched genes for sub-clusters of cluster: "
    , cluster_annot, " (", cluster_id, ")"
    , "\nEnrichment determined by DE of cells in sub-cluster versus other in cluster"
    , "\nx-axis: Genes"
    , "\ny-axis: Cells ordered by cluster"
    , "\nNormalized expression"
  ))
  ggsave(
    paste0(out_graph, "cluster", cluster_id,
     "_de_vs_subclust_expression_heatmap.png")
      , height = graph_height, width = 9)
}

run_plot_de_vs_subclust_expression_heatmaps <- function(){

  print("run_plot_de_vs_subclust_expression_heatmaps")

  cluster_annot_tb$cluster_number[11] %>%
    map(
      ~tryCatch({
          plot_de_vs_subclust_expression_heatmaps(
            seurat_obj = centSO
            , noCent40kExM = noCentExM
            , cluster_id = .
          )
        }
          , error = function(cond){
            message(paste0("Error for cluster: ", .))
            message(paste0(cond, "\n"))
            return(NA)
        }
    ))
}
################################################################################

### Odds ratio overlaps of cluster enriched genes

calculate_gene_lists_overlap_odds_ratio <- function(
  gene_list_1, gene_list_2, background_gene_list){

  print("calculate_gene_lists_overlap_odds_ratio")

  # binarize gene lists
  binary_gene_list_1 <- as.numeric(background_gene_list %in% gene_list_1)
  binary_gene_list_2 <- as.numeric(background_gene_list %in% gene_list_2)
  # glm
  gene_lists_glm <- glm(
    binary_gene_list_1~binary_gene_list_2, family = binomial)
  enrich_pval <- summary(gene_lists_glm)$coefficients[2,4]
  enrich_or <- summary(gene_lists_glm)$coefficients[2,1]
  # log2 odds ratio
  enrich_log2_or <- log2(exp(enrich_or))
  enrich_log2_or_tb <- tibble(
    enrich_log2_or = enrich_log2_or
    , pvalue = enrich_pval
    , number = sum(binary_gene_list_1 == 1 & binary_gene_list_2 == 1)
  )

  return(enrich_log2_or_tb)
}

calculate_gene_lists_overlap_odds_ratio_wrapper <- function(
  gene_lists_tb
  , background_gene_list
  , gene_col = "gene"
  , group_col = "group"
  , group_combn_tb){

    print("calculate_gene_lists_overlap_odds_ratio_wrapper")

    if (gene_col != "gene") {
      gene_lists_tb <- gene_lists_tb %>% rename(gene = gene_col)
    }
    if (group_col != "group") {
      gene_lists_tb <- gene_lists_tb %>% rename(group = group_col)
    }
    if (is.null(group_combn_tb)) {
      group_combn_tb <-
        expand.grid(
          gene_lists_tb %>% distinct(group) %>% pull()
          , gene_lists_tb %>% distinct(group) %>% pull()
        ) %>%
        as_tibble %>%
        mutate_all(as.character) %>%
        rename(group1 = Var1, group2 = Var2)
    }

    gene_lists_or_tb <- bind_cols(group_combn_tb, pmap(group_combn_tb, function(group1, group2){
        print(paste(group1, group2))
        gene_list_1 <-
          gene_lists_tb %>%
          filter(group == group1) %>%
          pull(gene)
        gene_list_2 <-
          gene_lists_tb %>%
          filter(group == group2) %>%
          pull(gene)
        enrich_log2_or <- calculate_gene_lists_overlap_odds_ratio(
          gene_list_1 = gene_list_1
          , gene_list_2 = gene_list_2
          , background_gene_list = background_gene_list)
        return(enrich_log2_or)
      }) %>%
      bind_rows)

    return(gene_lists_or_tb)
}

plot_gene_lists_overlaps_odds_ratios <- function(
  gene_lists_tb, dataset1, dataset2, background_gene_list
  , group2_levels = NULL){

  print("plot_gene_lists_overlaps_odds_ratios")

  # groups to calculate OR for
  group_combn_tb <-
    expand.grid(
      gene_lists_tb %>%
        filter(dataset == dataset1) %>%
        distinct(group) %>%
        pull
      , gene_lists_tb %>%
        filter(dataset == dataset2) %>%
        distinct(group) %>%
        pull
    ) %>%
    as_tibble %>%
    mutate_all(as.character) %>%
    rename(group1 = Var1, group2 = Var2)

  gene_lists_or_tb <- calculate_gene_lists_overlap_odds_ratio_wrapper(
    gene_lists_tb = gene_lists_tb
    , background_gene_list = background_gene_list
    , gene_col = "gene"
    , group_col = "group"
    , group_combn_tb = group_combn_tb
  )

  gg <- gene_lists_or_tb %>%
    mutate(group1 = gsub(paste0("_", dataset1), "", group1)) %>%
    mutate(group2 = gsub(paste0("_", dataset2), "", group2)) %>%
    mutate(cluster = gsub("_.*", "", group1) %>%
      factor(., levels = cluster_annot_tb$cluster_annot)
      ) %>%
    arrange(cluster) %>%
    mutate(group1 = factor(group1, levels = .$group1 %>% unique)) %>%
    {if(! is.null(group2_levels)){
      mutate(., group2 = factor(group2, levels = group2_levels))} else .} %>%
    mutate(pvalue = if_else(pvalue < 4.195783e-296, 4.195783e-296, pvalue)) %>%
    ggplot(aes(x = group1, y = group2, fill = -log(pvalue, 10))) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "-log10 p-value") +
      geom_text(aes(x = group1, y = group2, label = round(enrich_log2_or, 1))
        , color = "black", size = 3) +
      ggplot_set_theme_publication +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(gg)
}

run_plot_enriched_gene_overlap_odds_ratio_heatmaps <- function(){

  print("run_plot_enriched_gene_overlap_odds_ratio_heatmaps")

  ## Nowakowski

  # Nowakowski cluster DE table
  nowakowski_de_tb <- read_csv(
    "../data/nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv"
  )

  # make cluster enriched gene lists
  gene_lists_tb <- bind_rows(
    # nowakowski
    nowakowski_de_tb %>%
      # nowakowski DE not very cell type specifc, set DE filter higher
      filter(avg_diff > 0.75) %>%
      mutate(dataset = "nowakowski") %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene, enrichment = avg_diff, cluster, dataset, group)
    # geschwind
    , subclust_de_tb %>%
      filter(FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
      mutate(dataset = "geschwind") %>%
      mutate(Cluster = (as.character(Cluster))) %>%
      left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
      mutate(cluster = paste0(cluster_annot, "_", Subcluster)) %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene = Gene, enrichment = Log2_Fold_Change, cluster, dataset, group)
    )

  # background gene list
  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", 9, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)
  background_gene_list <- rownames(so@scale.data)

  plot_gene_lists_overlaps_odds_ratios(
    gene_lists_tb = gene_lists_tb
    , dataset1 = "geschwind"
    , dataset2 = "nowakowski"
    , background_gene_list = background_gene_list
  ) + ylab("Nowakowski et al. 2017") +
      xlab("Geschwind fetal drop-seq") +
      ggtitle(paste0(script_name
        , "\n\nOverlap of Nowakowski et al. 2017 and Geschwind fetal drop-seq cluster enriched genes"
        , "\nLog2 odds ratio (numbers)"))
  ggsave(paste0(out_graph, "overlap_geschwind_nowakowski_heatmap.png")
      , height = 12, width = 16)


  ## Lake

  # Lake 2017 cluster DE tables
  lake_de_tb <- read_csv("../source/Lake_2018_TS3_Cluster_DE_Genes.csv", skip = 4)
  # format
  lake_de_tb <- lake_de_tb %>% rename_all(clean_strings)

  # make cluster enriched gene lists
  gene_lists_tb <- bind_rows(
    # lake
    lake_de_tb %>%
      mutate(dataset = "lake") %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene, enrichment = average_difference_log_fold_change, cluster
        , dataset, group)
    # geschwind
    , subclust_de_tb %>%
      filter(FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
      mutate(dataset = "geschwind") %>%
      mutate(Cluster = (as.character(Cluster))) %>%
      left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
      mutate(cluster = paste0(cluster_annot, "_", Subcluster)) %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene = Gene, enrichment = Log2_Fold_Change, cluster, dataset, group)
    )

  # background gene list
  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", 9, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)
  background_gene_list <- rownames(so@scale.data)

  plot_gene_lists_overlaps_odds_ratios(
    gene_lists_tb = gene_lists_tb
    , dataset1 = "geschwind"
    , dataset2 = "lake"
    , background_gene_list = background_gene_list
  ) + ylab("Lake et al. 2017") +
      xlab("Geschwind fetal drop-seq") +
      ggtitle(paste0(script_name
        , "\n\nOverlap of Lake et al. 2017 and Geschwind fetal drop-seq cluster enriched genes"
        , "\nLog2 odds ratio (numbers)"))
  ggsave(paste0(out_graph, "overlap_geschwind_lake_heatmap.png")
      , height = 12, width = 16)


  ## dropseq sub-clusters vs dropseq clusters

  # make cluster enriched gene lists
  gene_lists_tb <- bind_rows(
    subclust_de_tb %>%
      filter(FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
      mutate(dataset = "geschwind_subclusters") %>%
      mutate(Cluster = (as.character(Cluster))) %>%
      left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
      mutate(cluster = paste0(cluster_annot, "_", Subcluster)) %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene = Gene, enrichment = Log2_Fold_Change, cluster, dataset, group)
    , cluster_de_tb %>%
      filter(FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
      mutate(dataset = "geschwind_clusters") %>%
      mutate(Cluster = (as.character(Cluster))) %>%
      left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
      mutate(cluster = cluster_annot) %>%
      mutate(group = paste0(cluster, "_", dataset)) %>%
      select(gene = Gene, enrichment = Log2_Fold_Change, cluster, dataset, group)
    )

  # background gene list
  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", 9, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)
  background_gene_list <- rownames(so@scale.data)

  plot_gene_lists_overlaps_odds_ratios(
    gene_lists_tb = gene_lists_tb
    , dataset1 = "geschwind_subclusters"
    , dataset2 = "geschwind_clusters"
    , background_gene_list = background_gene_list
    , group2_levels = cluster_annot_tb$cluster_annot
  ) + ylab("Geschwind fetal drop-seq") +
      xlab("Geschwind fetal drop-seq") +
      ggtitle(paste0(script_name
        , "\n\nOverlap of Geschwind fetal drop-seq sub-cluster enriched genes"
        , "\nLog2 odds ratio (numbers)"))
  ggsave(paste0(out_graph, "overlap_geschwind_clusters_heatmap.png")
      , height = 12, width = 16)


  ## dropseq sub-clusters vs dropseq sub-clusters

  # make cluster enriched gene lists
  gene_lists_tb <-
    subclust_de_tb %>%
    filter(FDR < 0.05 & Log2_Fold_Change > 0.2) %>%
    mutate(dataset = "geschwind") %>%
    mutate(Cluster = (as.character(Cluster))) %>%
    left_join(., cluster_annot_tb, by = c("Cluster" = "cluster_number")) %>%
    mutate(cluster = paste0(cluster_annot, "_", Subcluster)) %>%
    mutate(group = paste0(cluster, "_", dataset)) %>%
    select(gene = Gene, enrichment = Log2_Fold_Change, cluster, dataset, group)

  # background gene list
  # load seurat clustering round 2 object
  in_subclust_so <- paste0(
    "../analysis/analyzed_data/Seurat_ClusterRound2/ClusterRound2/20180907/"
    , "Seurat_ClusterRound2_Cluster", 9, "_seuratO.Robj")
  # so, noCentExM, rd1CentExM
  load(in_subclust_so)
  background_gene_list <- rownames(so@scale.data)

  plot_gene_lists_overlaps_odds_ratios(
    gene_lists_tb = gene_lists_tb
    , dataset1 = "geschwind"
    , dataset2 = "geschwind"
    , background_gene_list = background_gene_list
  ) + ylab("Geschwind fetal drop-seq") +
      xlab("Geschwind fetal drop-seq") +
      ggtitle(paste0(script_name
        , "\n\nOverlap of Geschwind fetal drop-seq sub-cluster enriched genes"
        , "\nLog2 odds ratio (numbers)"))
  ggsave(paste0(out_graph, "overlap_geschwind_heatmap.png")
      , height = 12, width = 16)

}
################################################################################

### Run

main_function()
################################################################################
