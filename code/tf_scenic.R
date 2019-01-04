# Damon Polioudakis
# 2018-11-25
# analysis of SCENIC TF motif enrichment

# must load modules:
#  module load R/3.51
#  module load gcc/7.2.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(AUCell)
require(tidyverse)
require(cowplot)
source("Function_Library.R")
source("GGplot_Theme.R")

### variables

## inputs
# seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

# SCENIC analysis done by Xu from Gerstein lab
load("../analysis/analyzed_data/Xu_TF_analysis/20181129/results_withJaspar/cells_AUC.RData")

# enrichment of TF regulons by Seurat cluster
scenic_enr_tb <- read_csv("../analysis/analyzed_data/Xu_TF_analysis/20181129/results_withJaspar/cells_enrich_TF_final.csv")

## other variables
script_name <- "tf_scenic.R"
date <- format(Sys.Date(), "%Y%m%d")
cluster_annot_tb <- tribble(
    # xu had old cluster annotations
    ~cluster_number, ~cluster_annot, ~cluster_annot_xu
    , "9",  "vRG", "vRG"
    , "7",  "oRG", "oRG"
    , "8",  "PgS", "Cycling progenitor S phase"
    , "10", "PgG2M", "Cycling progenitor G2/M phase"
    , "2",  "IP", "IPC"
    , "0",  "ExN", "Excitatory Neuron new born migrating"
    , "1",  "ExM", "Excitatory Neuron"
    # callosal is misspelled
    , "4",  "ExCal", "Excitatory Neuron (collosal)"
    , "3",  "ExDp1", "Deep layer excitatory neuron 1"
    , "13", "ExDp2", "Deep layer excitatory neuron 2"
    , "5",  "InSST", "Interneuron (SST)"
    , "6",  "InCALB2", "Interneuron (CALB2)"
    , "11", "OPC", "Oligodendrocyte precursor"
    , "12", "End", "Endothelial"
    , "14", "Per", "Pericyte"
    , "15", "Mic", "Microglia"
  )

## outputs
out_graph <- paste0(
  "../analysis/graphs/tf_scenic/", date, "/tf_scenic_")
out_scratch <- paste0(
  "/u/flashscratch/d/dpolioud/tf_scenic/", date, "/tf_scenic_")

# make output directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_scratch), recursive = TRUE)
################################################################################

### main function

main_function <- function(){
  plot_tsne_colored_by_auc()
  plot_heatmap_tf_regulon_cluster_enrichment()
}

################################################################################

### functions

plot_tsne_colored_by_variable <- function(
  tsne_1, tsne_2, variable_value, facet_variable
  , title = NULL, guide_size = 4, legend_title = NULL
  , alpha = 0.5, size = 0.1, ncol = 4
  , expression_color_gradient = FALSE){

  print("plot_tsne_colored_by_variable")

  gg_tb <- tibble(tsne_1 = tsne_1, tsne_2 = tsne_2
    , value = variable_value)
  if (! missing(facet_variable)){
    gg_tb <- mutate(gg_tb, facet_variable = facet_variable)}
  # Plot
  if (class(gg_tb$value) %in% c("character", "factor")){
    gg <- ggplot(gg_tb, aes(
        x = tsne_1, y = tsne_2, shape = value, col = value, group = value)) +
      geom_point(size = size, alpha = alpha, fill = NA) +
      scale_shape_manual(name = legend_title, values = rep(0:6,100)) +
      # option to facet
      { if (! missing(facet_variable)){
        facet_wrap(~facet_variable, scales = "free", ncol = ncol)
      }} +
      guides(color = guide_legend(title = legend_title
        , override.aes = list(size = guide_size, alpha = 1)))
  } else {
    gg <- ggplot(gg_tb, aes(x = tsne_1, y = tsne_2)) +
        geom_point(size = size, alpha = alpha, aes(col = value)) +
      scale_shape_manual(values = rep(0:6,100)) +
      # option to facet
      { if (! missing(facet_variable)){
        facet_wrap(~facet_variable, scales = "free", ncol = ncol)
      }} +
      # color scale options
      { if (expression_color_gradient == TRUE){
        # scale_color_gradientn(name = legend_title
        #   , colours = c(
        #     # desert
        #     # "lightgrey", "#fee090", "#fdae61", "#f46d43", "#ca0020")
        #     # colortest1
        #     # "lightgrey", "#d1e5f0", "#67a9cf", "#b2182b")
        #     # colortest2
        #     # "lightgrey", "#e0f3f8", "#74add1", "#d73027", "#a50026")
        #     # colortest3
        #     "#d9d9d9", "#bdbdbd", "#6baed6", "#4292c6", "#cb181d", "#99000d")
        #   )
        scale_colour_distiller(
          name = legend_title, type = "seq", palette = "PuRd", direction = 1)
      } else {
        scale_color_viridis(name = legend_title)
      }}
  }
    gg <- gg +
      ggplot_set_theme_publication +
      ggtitle(title)

  return(gg)
}
################################################################################

### tsne colored by AUC

plot_tsne_colored_by_auc <- function(){

  print("plot_tsne_colored_by_auc")

  getAUC(cells_AUC) %>%
    # gather auc scores
    t() %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(key = "gene_set", value = "auc", -cell) %>%
    # set max AUC to 0.5 for plotting
    mutate(auc = if_else(auc > 0.5, 0.5, auc)) %>%
    # add tsne coordinates
    left_join(.
      , as.data.frame(centSO@dr$tsne@cell.embeddings) %>%
        rownames_to_column("cell")
    ) %>%
    clean_variable_names(.) %>%
    # plot
    with(.,
      plot_tsne_colored_by_variable(
      tsne_1 = tsne_1, tsne_2 = tsne_2
      , variable_value = auc, facet_variable = gene_set
      , guide_size = 4, legend_title = "AUCell"
      , alpha = 0.5, size = 0.1, ncol = 4
      , expression_color_gradient = TRUE
      , title = paste0(script_name
        , "\n\ntSNE colored by SCENIC AUC"
        , "\nHuman fetal brain drop-seq dataset"))) +
        theme(
          , text = element_text(size = 10, colour = "black")
          , plot.title = element_text(size = 10)
          , axis.line = element_line(colour = "black")
          , axis.title = element_blank()
          , axis.text = element_blank()
          , axis.ticks = element_blank()
        )
  ggsave(paste0(out_graph, "tsne_auc.png")
    , width = 13, height = 90, limitsize = FALSE)
}



# gg_l <- map(rownames(getAUC(cells_AUC)), function(tf){
#   getAUC(cells_AUC)[rownames(getAUC(cells_AUC)) == tf, 1:1000, drop = FALSE] %>%
#     t() %>%
#     data.frame() %>%
#     rownames_to_column("cell") %>%
#     gather(key = "gene_set", value = "auc", -cell) %>%
#     mutate(gene_set = tf) %>%
#     # add tsne coordinates
#     left_join(.
#       , as.data.frame(centSO@dr$tsne@cell.embeddings) %>%
#         rownames_to_column("cell")
#       , by = c("cell" = "cell")
#     ) %>%
#     clean_variable_names(.) %>%
#     with(.,
#       plot_tsne_colored_by_variable(
#       tsne_1 = tsne_1, tsne_2 = tsne_2
#       , variable_value = auc, facet_variable = gene_set
#       , title = NULL, guide_size = 4, legend_title = NULL
#       , alpha = 0.5, size = 0.1
#       , expression_color_gradient = TRUE))
# })
# # extract legend - make sure it exists (NA genes plot with no legend)
# plot_legend <- lapply(gg_l, function(x) tryCatch(
#   get_legend(x), error = function(e) NA)) %>%
#     .[! is.na(.)] %>% .[[1]]
# # remove legends
# gg_l <- lapply(gg_l, function(gg){
#   gg + theme(legend.position = "none")})
# # plot
# Plot_Grid(gg_l, ncol = 4, rel_height = 0.1
#   , align = 'v', axis = 'r'
#   , title = paste0(script_name
#     , "\n\ntSNE colored by SCENIC AUC"
#     , "\nHuman fetal brain drop-seq dataset")
#   ) %>%
#   # add legend
#   plot_grid(., plot_legend, rel_widths = c(1, .1))
################################################################################

### TF regulon enrichment by Seurat cluster

plot_heatmap_tf_regulon_cluster_enrichment <- function(){

  print("plot_heatmap_tf_regulon_cluster_enrichment")

  # order TFs by cluster enrichment
  tf_order <-
    scenic_enr_tb %>%
    gather(cluster, pvalue, -TF) %>%
    clean_variable_names %>%
    # change NA pvalues to 1
    mutate(pvalue = if_else(is.na(pvalue), 1, pvalue)) %>%
    arrange(pvalue) %>%
    filter({! duplicated(.$tf)}) %>%
    mutate(cluster = factor(
      cluster, levels = cluster_annot_tb$cluster_annot_xu)) %>%
    arrange(cluster, tf) %>%
    # remove TFs with no cluster enrichment
    filter(pvalue < 0.05) %>%
    pull(tf)

  # use TF order to order by cluster enrichment and plot
  scenic_enr_tb %>%
    gather(cluster, pvalue, -TF) %>%
    clean_variable_names %>%
    mutate(pvalue = if_else(is.na(pvalue), 1, pvalue)) %>%
    # convert 0 p-values to smallest floating point
    mutate(pvalue = if_else(pvalue == 0, 3.978217e-320, pvalue)) %>%
    # -log10 transform p-values
    mutate(pvalue = -log(pvalue, 10)) %>%
    # set ordering
    mutate(tf = factor(tf, levels = tf_order)) %>%
    mutate(cluster = factor(
      cluster, levels = cluster_annot_tb$cluster_annot_xu)) %>%
    # plot
    ggplot(aes(x = cluster, y = tf, fill = pvalue)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red"
        , space = "Lab", name = "-log10 p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste0(script_name
        , "\nTF regulon cluster enrichment from Xu"
        , "\nFishers test using AUC score cutoff for each cell"))
      ggsave(paste0(out_graph, "enrichment_heatmap.png"), height = 20)
      ggsave(paste0(out_graph, "enrichment_heatmap_paper.png")
        , height = 20, width = 25)
}
################################################################################

### run

main_function()
################################################################################
