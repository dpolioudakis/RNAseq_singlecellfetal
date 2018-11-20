# Damon Polioudakis
# Function library

################################################################################

### Differential expression

# Filter expression matrix by:
# percent of cells in cluster gene is expressed in
# fold change of gene in cluster versus all other cells
DE_Filters_ClustersAvsB_ExpMatrix <- function(
  so
  , minPercent = NULL
  , foldChange = NULL
  , clusterIDs1 = NULL
  , clusterIDs2 = NULL
) {

  if (! is.null(minPercent)) {
    # Expressed > 0 counts in > X% of cells in cluster
    if (! is.null(clusterIDs1)) {
      # Subset expression matrix to cluster
      cdf <- as.matrix(so@data)[ ,so@ident %in% clusterIDs1]
    }
    # Expressed > 0 counts in > X% of cells in cluster
    idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (minPercent / 100)
    print(paste0("Genes expressed in > ", minPercent, "% of cells in cluster"))
    print(table(idxp))
  } else {
    idxp <- rep(TRUE, nrow(so@data))
  }

  if (! is.null(foldChange)) {
    # Fold change > Y of gene in cluster versus all other cells
    if (! is.null(clusterIDs1)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,so@ident %in% clusterIDs1]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ ,so@ident %in% clusterIDs2]
    }
    # Fold change
    v1 <- rowMeans(cdf) - rowMeans(ndf)
    idxf <- v1 > foldChange
    print(paste0("Genes > ", foldChange, " fold change in cluster versus all other cells"))
    print(table(idxf))
  } else {
    idxf <- rep(TRUE, nrow(so@data))
  }

  # Filter exDF
  exDF <- as.matrix(so@data[idxp & idxf, so@ident %in% c(clusterIDs1, clusterIDs2)])
  return(exDF)
}

# Filter expression matrix by:
# percent of cells in cluster gene is expressed in
# fold change of gene in cluster versus all other cells
DE_Filters_ExpMatrix <- function(
  expr_m
  , minPercent = NULL
  , foldChange = NULL
  , clusterID = NULL
  , cell_cluster_key
  , cellID = NULL) {

  if (! is.null(minPercent)) {
    # Expressed > 0 counts in > X% of cells in cluster
    if (! is.null(clusterID)) {
      # Subset expression matrix to cluster
      cdf <- as.matrix(expr_m)[ ,cell_cluster_key == clusterID]
    }
    if (! is.null(cellID)) {
      # Subset expression matrix to cluster
      cdf <- as.matrix(expr_m)[ ,colnames(expr_m) %in% cellID]
    }
    # Expressed > 0 counts in > X% of cells in cluster
    idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (minPercent / 100)
    print(paste0("Genes expressed in > ", minPercent, "% of cells in cluster"))
    print(table(idxp))
  } else {
    idxp <- rep(TRUE, nrow(expr_m))
  }

  ### not currently generalized outside of Seurat
  if (! is.null(foldChange)) {
    # Fold change > Y of gene in cluster versus all other cells
    if (! is.null(clusterID)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,cell_cluster_key == clusterID]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! cell_cluster_key == clusterID]
    }
    if (! is.null(cellID)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,colnames(expr_m) %in% cellID]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! colnames(expr_m) %in% cellID]
    }
    # Fold change
    v1 <- rowMeans(cdf) - rowMeans(ndf)
    idxf <- v1 > foldChange
    print(paste0("Genes > ", foldChange, " fold change in cluster versus all other cells"))
    print(table(idxf))
  } else {
    idxf <- rep(TRUE, nrow(expr_m))
  }

  # Filter exDF
  exDF <- as.matrix(expr_m[idxp & idxf, ])
  return(exDF)
}

## Function: DE Linear model
# termsDF:
# ExpCondition RIN.y     Seq.PC1      Seq.PC2
# 1            CP   8.4  0.04792498 -0.090448567
# 2            CP   8.1  0.53502697 -0.287629654
# 3            CP   8.0  0.18824922 -0.155651102
# 4            VZ   8.4  0.02529722 -0.100858264
# 5            VZ   8.7  0.45139297  0.856908177
# 6            VZ   9.1  0.27861748 -0.248868277
# mod: "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2"
DE_Linear_Model <- function(exDatDF, termsDF, mod) {
  print("DE_Linear_Model")
  lmmod <- apply(as.matrix(exDatDF), 1
    , function(y) {
      mod <- as.formula(mod)
      lm(mod, data = termsDF)})
  coefmat <- matrix(NA, nrow = nrow(exDatDF)
    , ncol = length(coef(lmmod[[1]])))
  pvalmat <- matrix(NA, nrow = nrow(exDatDF)
    , ncol = length(summary(lmmod[[1]])[[4]][ ,4]))
  colnames(coefmat) <- names(coef(lmmod[[1]]))
  rownames(coefmat) <- rownames(exDatDF)
  colnames(pvalmat) <- names(summary(lmmod[[1]])[[4]][ ,4])
  rownames(pvalmat) <- rownames(exDatDF)
  for (i in 1:nrow(exDatDF)) {
    if (i%%100 == 0) {cat(".")}
    coefmat[i, ] <- coef(lmmod[[i]])
    pvalmat[i, ] <- summary(lmmod[[i]])[[4]][ ,4]
  }
  deCoefPvalLM <- list(coefmat = coefmat, pvalmat = pvalmat)
  return(deCoefPvalLM)
}

# Format output of linear model into data frame
Format_DE <- function(
  deLM, expr_m, clusterID, cell_cluster_key, cluster_annot = NULL) {
  print("Format_DE")
  # Combine log2 fold changes, p-values
  deDF <- data.frame(Gene = row.names(deLM$coefmat)
    , Log2_Fold_Change = deLM$coefmat[ ,2]
    , Pvalue = deLM$pvalmat[ ,2])
  # Add cluster ID
  deDF$Cluster <- clusterID
  # Add cluster annotation
  if (! is.null(cluster_annot_tb)){
    deDF$Cluster_Annot <- cluster_annot
  }
  # Percent of cells in cluster expressing gene > 0 counts
  cdf <- expr_m[
    row.names(expr_m) %in% deDF$Gene, cell_cluster_key == clusterID]
  deDF$Percent_Cluster <- (rowSums(cdf > 0) / ncol(cdf)) * 100
  deDF$Percent_Cluster <- round(deDF$Percent_Cluster, 1)
  # Percent of all cells expressing gene > 0 counts
  deDF$Percent_All <- (rowSums(expr_m[
    row.names(expr_m) %in% deDF$Gene, ] > 0)
    / ncol(expr_m[row.names(expr_m) %in% deDF$Gene, ])) * 100
  deDF$Percent_All <- round(deDF$Percent_All, 1)
  # Order by log fold change
  deDF <- deDF[order(-deDF$Log2_Fold_Change), ]
  return(deDF)
}
################################################################################

### Expression heatmaps

# Format data for expression heatmap, cells ordered by cluster
Heatmap_By_Cluster_Format_Data <- function(
  geneGroupDF, exprM, cellID_clusterID, clusters, lowerLimit, upperLimit
  , geneOrder = FALSE, clusterOrder = NULL) {

  print("Heatmap_By_Cluster_Format_Data")

  # Subset expression matrix to genes of interest by merging
  ggDF <- merge(geneGroupDF[c("Gene", "Group")], exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)

  if (geneOrder == TRUE) {
    # levels <- expand.grid(geneOrder, unique(ggDF$Group))
    # levels <- paste0(levels$Var1, "   ", levels$Var2)
    gene_levels <- paste0(geneGroupDF$Gene, "   ", geneGroupDF$Group)
  }

  # Set group factor levels
  ggDF$Group <- factor(ggDF$Group
    , levels = unique(as.character(geneGroupDF$Group))[unique(as.character(geneGroupDF$Group)) %in% as.character(ggDF$Group)])
  # colnames(ggDF)[1:2] <- c("Gene", "Group")
  # Remove blanks
  ggDF <- ggDF[! ggDF$Gene == "", ]
  # Format for ggplot
  ggDF <- melt(ggDF)
  # Add seurat clusters
  idx <- match(as.character(ggDF$variable), names(cellID_clusterID))
  ggDF$Cluster <- cellID_clusterID[idx]
  # Set cluster levels
  if (! is.null(clusterOrder)) {
    ggDF$Cluster <- factor(ggDF$Cluster, levels = clusterOrder)
  }
  # Subset clusters
  ggDF <- ggDF[as.character(ggDF$Cluster) %in% as.character(clusters), ]
  # Add group to gene name
  ggDF$Gene_Group <- paste0(ggDF$Gene, "   ", ggDF$Group)

  # Order genes by setting levels
  if (! geneOrder == FALSE) {
    ggDF$Gene_Group <- factor(ggDF$Gene_Group, levels = gene_levels)
  }
  # Set limits
  ggDF$value[ggDF$value < lowerLimit] <- lowerLimit
  ggDF$value[ggDF$value > upperLimit] <- upperLimit

  return(ggDF)
}

# Expression heatmap, cells ordered by cluster
Heatmap_By_Cluster <- function(
  geneGroupDF, exprM, cellID_clusterID, clusters, lowerLimit, upperLimit
  , geneOrder = FALSE, clusterOrder = NULL) {

  print("Heatmap_By_Cluster")

  ggDF <- Heatmap_By_Cluster_Format_Data(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , cellID_clusterID = cellID_clusterID
    , clusters = clusters
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneOrder
    , clusterOrder = clusterOrder
  )

  print("Heatmap_By_Cluster: plotting...")
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = Gene_Group, fill = value)) +
    geom_tile() +
    facet_grid(Group~Cluster, space = "free", scales = "free") +
    # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(lowerLimit, upperLimit)) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(strip.background = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10)) +
    ylab("Genes") +
    xlab("Cells ordered by cluster")
  # gg <- gg + ...
  return(gg)
}

# Wrapper function for Heatmap_By_Cluster
# Uses plot_grid to combine 3 heatmaps to deal with cell number scaling
Heatmaps_By_Cluster_Combined <- function(geneGroupDF, exprM, cellID_clusterID
  , clusters1, clusters2, clusters3, lowerLimit, upperLimit, geneOrder = NULL) {

  p1 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , cellID_clusterID = cellID_clusterID, clusters = clusters1
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder = geneOrder
  )
  p1 <- p1 + theme(
    # axis.title.x = element_blank()
    , strip.text.y = element_blank()
    , legend.position = "none"
    )
  p2 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , cellID_clusterID = cellID_clusterID, clusters = clusters2
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder = geneOrder
  )
  p2 <- p2 + theme(
    strip.text.y = element_blank()
    , legend.position = "none"
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank()
  )
  p3 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , cellID_clusterID = cellID_clusterID, clusters = clusters3
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder = geneOrder
  )
  p3 <- p3 + theme(
    # axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank()
    )
  ggL <- list(p1, p2, p3)
  # # plot_grid combine
  # pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # # now add the title
  # title <- ggdraw() + draw_label(title)
  # # rel_heights values control title margins
  # pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  return(ggL)
}

# Expression heatmap of list of genes, cells ordered by 1st gene
Heatmap_By_Gene_Expression <- function (genes, exprM, limLow, limHigh) {
  ggDF <- merge(genes, exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF <- melt(ggDF)
  ggDF$x <- factor(ggDF$x, levels = genes)
  ggDF <- ggDF[with(ggDF, order(x, -value)), ]
  ggDF$variable <- factor(ggDF$variable, levels = unique(ggDF$variable))
  ggDF$value[ggDF$value > limHigh] <- limHigh
  ggDF$value[ggDF$value < limLow] <- limLow
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = x, fill = value)) +
    geom_tile() +
    # facet_grid(Group~Cluster, space = "free", scales = "free") +
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10))
  return(gg)
}

# Expression heatmap, cells ordered by cluster
Plot_Marker_Genes_Heatmap_SetColWidths <- function(
  geneGroupDF
  , exprM
  , cellID_clusterID
  , clusters = c(0:15)
  , lowerLimit = -1.5
  , upperLimit = 1.5
  , geneOrder = TRUE
  , clusterOrder = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
  , color_viridis = FALSE
  ) {

  print("Plot_Marker_Genes_Heatmap_SetColWidths")
  print(str(exprM))

  # Heatmap plot
  # Normalized, mean centered and scaled
  ggDF <- Heatmap_By_Cluster_Format_Data(
    geneGroupDF = geneGroupDF
    , exprM = exprM
    , cellID_clusterID = cellID_clusterID
    , clusters = clusters
    , lowerLimit = lowerLimit
    , upperLimit = upperLimit
    , geneOrder = geneOrder
    , clusterOrder = clusterOrder
  )

  print("Heatmap_By_Cluster: plotting...")
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = Gene_Group, fill = value)) +
    geom_tile() +
    facet_grid(Group~Cluster, space = "free_y", scales = "free"
      , drop = TRUE) +
    { if(color_viridis == TRUE) {
      scale_fill_viridis(name = "Normalized expression"
      , limits = c(lowerLimit, upperLimit))
    } else {
      scale_fill_distiller(name = "Normalized\nexpression\nzscore"
        , type = "div", palette = 5, direction = -1
        , limits = c(lowerLimit, upperLimit)
        , na.value = "grey90")
    }} +
    # scale_fill_distiller(name = "Normalized\nexpression", type = "div"
    #   , palette = 5, direction = -1, limits = c(lowerLimit, upperLimit)
    #   , na.value = "grey90") +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(strip.background = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10)) +
    ylab("Genes") +
    xlab("Cells ordered by cluster")
  # gg <- gg + ...
  return(gg)
}

# # Example:
# # Normalized, no mean centering scaling
# p1 <- Heatmap_By_Gene_Expression(genes = genes, exprM = centSO@scale.data
#   , limLow = -1.5, limHigh = 1.5)
# p1 <- p1 + ylab("Genes") +
#   # scale_fill_distiller(name = "Normalized\nexpression\nz-score") +
#   ylab("Genes") +
#   xlab("Cells ordered by PAX6 expression") +
#   ggtitle("Normalized mean centered and scaled expression")
#
# # Normalized
# p2 <- Heatmap_By_Gene_Expression(genes = genes, exprM = noCentExM
#   , limLow = -1, limHigh = 3)
# p2 <- p2 + ylab("Genes") +
#   xlab("Cells ordered by PAX6 expression") +
#   ggtitle("Normalized expression")
#
# # plot_grid
# pg <- plot_grid(p1, p2, ncol = 2)
# # now add the title
# title <- paste0(graphCodeTitle
#   , "\n\nExpression of neuronal differentiation marker genes"
#   , "\nx-axis: Genes"
#   , "\ny-axis: Cells ordered by PAX6 expression"
#   , "\n")
# title <- ggdraw() + draw_label(title)
# # rel_heights values control title margins
# pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.25, 1))
# ggsave(paste0(outGraph, "Markers_Heatmap.png"), width = 12, height = 7)
################################################################################

### Feature plot

# Calculate mean expression of group of genes for each cell using seurat scaled
# expression values
Mean_Expression <- function(tsneDF, genes, exM) {
  genesExDF <- exM[which(row.names(exM) %in% genes), ]
  # genesExDF <- as.matrix(seuratO@data)[which(row.names(as.matrix(seuratO@data)) %in% genes), ]
  # Only calculate column means if there are multiple genes
  print("Length genes:")
  print(length(genes))
  if (is.matrix(genesExDF)) {
    mnExDF <- colMeans(genesExDF)
  } else {
    mnExDF <- genesExDF
  }
  # Add to ggplot data frame
  tsneDF$EXPRESSION <- mnExDF[match(row.names(tsneDF), names(mnExDF))]
  return(tsneDF)
}

# Transform data to desired limits for ggplot2
Set_Limits <- function(tsneDF, limHigh, limLow) {
  tsneDF$EXPRESSION[tsneDF$EXPRESSION < limLow] <- limLow
  tsneDF$EXPRESSION[tsneDF$EXPRESSION > limHigh] <- limHigh
  return(tsneDF)
}

# Color tSNE plot by expression from Mean_Expression()
FeaturePlot_Graph <- function(tsneDF, title, limLow, limHigh
  , alpha = 0.5, size = 0.05) {
  ggFp <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    # geom_point(size = size, alpha = alpha) +
    geom_point(size = size, alpha = alpha, stroke = 0.5, shape = 21) +
    # scale_color_distiller(name = "Normalized\nexpression", type = "div"
    #   , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    scale_color_distiller(name = "Normalized\nexpression", type = "seq"
      , palette = "BuPu", direction = 1, limits = c(limLow, limHigh)) +
    # scale_color_viridis(option = "magma", limits = c(limLow, limHigh)) +
    ggtitle(title)
  return(ggFp)
}

# Color tSNE plot by expression from Mean_Expression()
FeaturePlot_Graph_CentScale <- function(tsneDF, title, limLow, limHigh
  , alpha = 0.5, size = 0.02) {
  ggFp <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = size, alpha = alpha, stroke = 0.5, shape = 21) +
    scale_color_distiller(name = "Normalized\nexpression\nz-score", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    ggtitle(title)
  return(ggFp)
}

# tSNE plot colored by Seurat clustering
TSNE_Plot <- function(seuratO, do.label = TRUE) {
  # tSNE graph colored by cluster
  gg_tsne <- TSNEPlot(seuratO, do.label = do.label, pt.size = 0.02
    , do.return = TRUE, no.legend = FALSE)
  gg_tsne <- gg_tsne + ggtitle(paste0("Color indicates cluster assignment"))
  gg_tsne <- gg_tsne + theme_set(theme_bw()) +
    theme_set(theme_get() + theme(text = element_text(size = 16))) +
    theme_update(plot.title = element_text(size = 12)) +
    theme_update(axis.line = element_line(colour = "black")
      , plot.background = element_blank()
      , panel.border = element_blank()
    )
  return(gg_tsne)
}
# ggL <- append(list(gg_tsne), ggL)

# Wrapper function for feature plots
FeaturePlot_CentScale <- function(genes, tsneDF, seuratO, limLow, limHigh) {
  # tsneDF: as.data.frame(centSO@dr$tsne@cell.embeddings)
  # Collect tSNE values for ggplot
  tsneDF <- as.data.frame(seuratO@dr$tsne@cell.embeddings)
  # Loop through and plot each group of genes
  ggL <- lapply(genes, function(gene) {
    print(gene)
    tsneDF <- Mean_Expression(tsneDF, gene, seuratO@scale.data)
    tsneDF <- Set_Limits(tsneDF, limLow = limLow, limHigh = limHigh)
    ggFp <- FeaturePlot_Graph_CentScale(tsneDF, title = paste0("\n", gene)
      , limLow = limLow, limHigh = limHigh)
    return(ggFp)
  })
  gg_tsne <- TSNE_Plot(seuratO) +
    guides(color = guide_legend(title = "Cluster")) +
    theme(legend.position = "none")
  gg_tsne_nolabel <- TSNE_Plot(seuratO, do.label = FALSE) +
    guides(color = guide_legend(title = "Cluster")) +
    theme(legend.position = "none")
  ggL <- append(list(gg_tsne_nolabel), ggL)
  ggL <- append(list(gg_tsne), ggL)
  return(ggL)
}
# Wrapper function for feature plots
FeaturePlot <- function(genes, tsneDF, seuratO, exM, limLow, limHigh
  , size = 0.02, alpha = 0.5, geneGrouping = NULL, centScale = FALSE) {
  # Collect tSNE values for ggplot
  # tsneDF <- as.data.frame(seuratO@dr$tsne@cell.embeddings)
  # tSNE_1    tSNE_2
  # ACCTAAGGATTA  4.872695 13.282958
  # CCGTTTGTGATA  6.102554  3.158785
  # GGCACAAGTGGC  3.886729  7.046055
  # Loop through and plot each group of genes

  # Mean expression of gene group if grouping is provided
  if (! is.null(geneGrouping)) {

    genesGroupingDF <- data.frame(Gene = genes, Grouping = geneGrouping)

    ggL <- lapply(unique(genesGroupingDF$Grouping), function(grouping) {
      print(grouping)
      genes <- genesGroupingDF$Gene[genesGroupingDF$Grouping == grouping]
      print(head(genes))
      tsneDF <- Mean_Expression(tsneDF, genes, exM)
      tsneDF <- Set_Limits(tsneDF, limLow = limLow, limHigh = limHigh)
      if (centScale == FALSE) {
        ggFp <- FeaturePlot_Graph(tsneDF, title = paste0("\n", grouping)
          , limLow = limLow, limHigh = limHigh, size = size, alpha = alpha)
      }
      if (centScale == TRUE) {
        ggFp <- FeaturePlot_Graph_CentScale(tsneDF, title = paste0("\n", grouping)
          , limLow = limLow, limHigh = limHigh, size = size, alpha = alpha)
      }
      return(ggFp)
    })
  }

  # If grouping is not provided plot each gene individually
  else {
    ggL <- lapply(genes, function(gene) {
      print(gene)
      tsneDF <- Mean_Expression(tsneDF, gene, exM)
      tsneDF <- Set_Limits(tsneDF, limLow = limLow, limHigh = limHigh)
      if (centScale == FALSE) {
        ggFp <- FeaturePlot_Graph(tsneDF, title = paste0("\n", gene)
          , limLow = limLow, limHigh = limHigh, size = size, alpha = alpha)
      }
      if (centScale == TRUE) {
        ggFp <- FeaturePlot_Graph_CentScale(tsneDF, title = paste0("\n", gene)
          , limLow = limLow, limHigh = limHigh, size = size, alpha = alpha)
      }
      return(ggFp)
    })
  }

  gg_tsne <- TSNE_Plot(seuratO) + theme(legend.position = "none")
  gg_tsne_nolabel <- TSNE_Plot(seuratO, do.label = FALSE) +
    theme(legend.position = "none")
  ggL <- append(list(gg_tsne_nolabel), ggL)
  ggL <- append(list(gg_tsne), ggL)
  # Format
  ggL <- lapply(ggL, function(gg){gg + ggplot_set_theme_publication})
  ggL[1:2] <- lapply(ggL[1:2], function(gg){
      gg + guides(color = guide_legend(
        ncol = 2, title = "Cluster", override.aes = list(size = 3)))
  })
  # Output
  return(ggL)
}

# # Genes to plot
# genes <- c("PAX6", "EOMES")
# # FeaturePlot outputs list of ggplots
# ggL <- FeaturePlot_CentScale(genes, tsneDF, centSO, -1.5, 1.5)
# # plot_grid combine tSNE graphs
# pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'v', axis = 'r')
# # now add the title
# title <- ggdraw() + draw_label(paste0(graphCodeTitle
#   , "\n\nExpression of genes of interest"
#   , "\nRemove cells < 200 genes detected"
#   , "\nRemove cells > 3192 (3 SD) genes detected"
#   , "\nRemove genes detected in < 3 cells"
#   , "\nRegress out nUMI, donor, library lab"
#   , "\nNormalized expression, mean centered and variance scaled"
#   , "\ntSNE PC 1-40"
#   , "\n"))
# # rel_heights values control title margins
# plot_grid(title, pg, ncol = 1, rel_heights = c(0.5, 1))
# ggsave(paste0(outGraph, "FeaturePlot_PAX6_EOMES_NormalizedCenteredScaled.png")
#   , width = 14, height = 3+length(ggL)*1.25)
################################################################################

### General functions

# Add gene list as column to 0 1 matrix with genes as row names
Add_Gene_List_To_Binary_Matrix <- function(
  gene_binary_M, gene_list, gene_list_name){
  print("Add_Gene_List_To_Binary_Matrix")
  binary_gene_list <- as.numeric(rownames(gene_binary_M) %in% gene_list)
  gene_binary_M <- cbind(gene_binary_M, binary_gene_list)
  colnames(gene_binary_M)[dim(gene_binary_M)[2]] <- gene_list_name
  return(gene_binary_M)
}

Convert_Mixed_GeneSym_EnsID_To_EnsID <- function(ids, bmDF = bmDF){
  # Converts vector of gene symbols and ensembl IDs to ensembl IDs
  # ids must be character format
  # Need to load:
  # bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  # , header = TRUE)
  print("Convert_Mixed_GeneSym_EnsID_To_EnsID")
  idx <- match(ids, bmDF$hgnc_symbol)
  ens <- bmDF$ensembl_gene_id[idx]
  ids[! is.na(ens)] <- as.character(ens[! is.na(ens)])
  return(ids)
}

clean_variable_names <- function(data){
  cleaned <- data %>%
    rename_all(
      funs(clean_strings)
    )
  return(cleaned)
}

clean_strings <- function(string_vector){
  print("clean_strings")
  cleaned <- string_vector %>%
    gsub("* ", "_", .) %>%
    gsub("\\.", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "_", .) %>%
    gsub("\\+", "and", .) %>%
    gsub("#", "_number", .) %>%
    gsub("_$", "", .) %>%
    gsub("__", "_", .) %>%
    tolower
  return(cleaned)
}

# Cowplot plot_grid and add title
Plot_Grid <- function(ggPlotsL, ncol = 2, title = "", rel_height = 0.1, ...) {
  # cowplot plot_grid ...: align = 'v', axis = 'l'
  # Plot grid
  pg <- plot_grid(plotlist = ggPlotsL, ncol = ncol, ...)
  # now add the title
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(rel_height, 1))
}

# Remove outliers from a column based on SD
Remove_Outliers_By_SD <- function(x, nStdev, na.rm = TRUE, ...) {
  stdev <- sd(x)
  mn <- mean(x)
  x[x > (mn + nStdev*stdev) | x < (mn - nStdev*stdev)] <- NA
  # x <- as.numeric(x)
  return(x)
}
################################################################################

# ggplot

# Publication theme
ggplot_set_theme_publication <-
  theme_bw() +
  theme(
    , text = element_text(size = 10, colour = "black")
    , plot.title = element_text(size = 10)
    , axis.line = element_line(colour = "black")
    , panel.border = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  )
ggplot_set_theme_publication_nolabels <-
  theme_bw() +
  theme(
    , text = element_text(size = 10, colour = "black")
    , plot.title = element_text(size = 10)
    , axis.line = element_line(colour = "black")
    , panel.border = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , axis.title = element_blank()
    , axis.text = element_blank()
    , axis.ticks = element_blank()
  )
################################################################################

### Load datasets

# ASD
# TADA Sanders 2015 = from Luis metaMat
load_asd_sanders_genes <- function(){
  read_csv(
    "../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv") %>%
    filter(tadaFdrAscSscExomeSscAgpSmallDel < 0.1) %>%
    select(gene = RefSeqGeneName) %>%
    clean_variable_names()
}

# ihart
load_asd_ihart_genes <- function(){
  read_csv(
    "../source/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv") %>%
    clean_variable_names() %>%
    mutate(hgnc_gene_symbol = gsub("\"", "", hgnc_gene_symbol)) %>%
    filter(ihart_69 == 1) %>%
    select(hgnc_gene_symbol)
}


# Epilepsy high confidence risk genes from Elizabeth Ruzzo
load_epilepsy_high_conf_genes <- function(){
  read_tsv(
    "../source/Gene_Lists/High-Confidence_Epilepsy_Risk_Genes_Ruzzo_2018-05-11.txt") %>%
    clean_variable_names() %>%
    filter(classification == "High-confidence") %>%
    select(gene)
}

# ID risk genes
# de novo ID: NEJM + Lancet
load_id_nejm_lancet_genes <- function(){
  id_genes_tb <- bind_rows(
    read_tsv("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_deLigt_NEJM.txt") %>%
      clean_variable_names() %>%
      filter(nature_of_mutation == "D") %>%
      select(gene)
    , read_tsv("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_Rauch_Lancet.txt") %>%
      clean_variable_names() %>%
      filter(type %in% c("frameshift", "nonsense", "splice")) %>%
      select(gene = gene_symbol) %>%
      rename()
  ) %>% distinct()
  return(id_genes_tb)
}

## TFs, chromatin remodelers, and co-factors
load_tf_cofactors_remodelers <- function(){
  print("load_tf_cofactors_remodelers")

  # Load
  factors_tb <- bind_rows(
    # Human TFs
    read_tsv("../source/AnimalTFDB_Homo_sapiens_TF_EnsemblID.txt"
      , col_names = FALSE) %>%
      add_column(type = "TF")

    # Human chromatin remodeling factors
    , read_tsv("../source/AnimalTFDB_Homo_sapiens_chr_remodeling_factor_EnsemblID.txt"
      , col_names = FALSE) %>%
      add_column(type = "chromatin remodeler")

    # Human co-factors
    , read_tsv("../source/AnimalTFDB_Homo_sapiens_cofactor_EnsemblID.txt"
      , col_names = FALSE) %>%
      add_column(type = "cofactor")
  )
  # Add hgnc symbols
  bm_tb <- read_csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
    , col_names = TRUE)
  factors_tb <- bm_tb %>% select(ensembl_gene_id, hgnc_symbol) %>%
    right_join(factors_tb, by = c("ensembl_gene_id" = "X1"))

  return(factors_tb)
}
################################################################################

### Monocle functions

# Format expression of gene list for plotting trajectory
# For Plot_Trajectory_Gene_Expression
Add_Expr_To_MonocleObj <- function (monocleO, genes, exprM, limHigh, limLow) {

  # Subset expression matrix to genes of interest
  m <- exprM
  m <- m[row.names(m) %in% genes, ]

  # Mean expression
  if (class(m) != "numeric") {
    m <- colMeans(m)
  }

  # Set expression limits
  m[m > limHigh] <- limHigh
  m[m < limLow] <- limLow

  # Merge with monocle cells to filter any cells not in monocle object
  m <- data.frame(CELLS = names(m), EXPRESSION = m)
  m <- merge(pData(monocleO), m, by.x = "row.names"
    , by.y = "CELLS")

  # Order cells like monocle object
  idx <- match(row.names(pData(monocleO)), m$CELL)
  m <- m[idx, ]

  # Add column with expression values to monocle object pData slot
  pData(monocleO)$EXPRESSION = m$EXPRESSION
  return(monocleO)
}

# Color monocle trajectory by expression of gene list
Plot_Trajectory_Gene_Expression <- function (monocleO, genes, exprM
  , limHigh, limLow, title, cell_size = 0.01, ...) {

  mo <- Add_Expr_To_MonocleObj(monocleO = monocleO, genes = genes
    , exprM = exprM, limHigh, limLow)

  gg <- plot_cell_trajectory(mo, color_by = "EXPRESSION", cell_size = cell_size, ...)
  gg <- gg +
    scale_color_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    theme(legend.position = "right") +
    ggtitle(title)
  return(gg)
}

# Expression heatmap, cells ordered by state
Heatmap_By_State <- function(
  geneGroupDF, exprM, monocleO, clusters, lowerLimit, upperLimit
  , geneOrder = NULL, centScale = NULL) {

  # Subset expression matrix to genes of interest by merging
  ggDF <- merge(geneGroupDF[c("Gene", "Group")], exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  # Set group factor levels
  ggDF$Group <- as.factor(ggDF$Group)
  # colnames(ggDF)[1:2] <- c("Gene", "Group")
  # Remove blanks
  ggDF <- ggDF[! ggDF$Gene == "", ]
  # Save order to set levels later
  if (is.null(geneOrder)) {
    levels <- paste0(ggDF$Gene, "   ", ggDF$Group)
  } else {
    levels <- rev(paste0(geneOrder, "   ", ggDF$Group))
  }
  if (centScale == TRUE) {
    Gene <- ggDF$Gene
    Group <- ggDF$Group
    ggDF <- as.data.frame(t(scale(t(ggDF[ ,-c(1:2)]))))
    ggDF$Gene <- Gene
    ggDF$Group <- Group
  }
  ggDF <- melt(ggDF)
  # Add seurat clusters
  idx <- match(ggDF$variable, pData(monocleO)$CELL)
  ggDF$MONOCLE_STATE <- pData(monocleO)$State[idx]
  # Subset clusters
  ggDF <- ggDF[ggDF$MONOCLE_STATE %in% clusters, ]
  # Add group to gene name and order by levels
  ggDF$Gene_Group <- paste0(ggDF$Gene, "   ", ggDF$Group)
  ggDF$Gene_Group <- factor(ggDF$Gene_Group, levels = levels)
  # Set limits
  ggDF$value[ggDF$value < lowerLimit] <- lowerLimit
  ggDF$value[ggDF$value > upperLimit] <- upperLimit
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = Gene_Group, fill = value)) +
    geom_tile() +
    facet_grid(Group~MONOCLE_STATE, space = "free", scales = "free") +
    # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1, limits = c(lowerLimit, upperLimit)) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(strip.background = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 10)) +
    ylab("Genes") +
    xlab("Cells ordered by cluster")
  # gg <- gg + ...
  return(gg)
}

# Wrapper function for Heatmap_By_State
# Uses plot_grid to combine 3 heatmaps to deal with cell number scaling
Heatmaps_By_State_Combined <- function(geneGroupDF, exprM, monocleO
  , clusters1, clusters2, clusters3, lowerLimit, upperLimit, geneOrder = NULL) {

  p1 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , monocleO = monocleO, clusters = clusters1
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder
  )
  p1 <- p1 + theme(
    axis.title.x = element_blank()
    , strip.text.y = element_blank()
    , legend.position = "none"
  )
  p2 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , monocleO = monocleO, clusters = clusters2
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder
  )
  p2 <- p2 + theme(
    strip.text.y = element_blank()
    , legend.position = "none"
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank()
  )
  p3 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , monocleO = monocleO, clusters = clusters3
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder
  )
  p3 <- p3 + theme(
    axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , axis.text.y = element_blank()
    , axis.ticks.y = element_blank())
  ggL <- list(p1, p2, p3)
  # # plot_grid combine
  # pg <- plot_grid(plotlist = ggL, ncol = 3, align = 'h', axis = 'b')
  # # now add the title
  # title <- ggdraw() + draw_label(title)
  # # rel_heights values control title margins
  # pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
  return(ggL)
}
################################################################################

ScatterPlot_Expression <- function (genes, exprM, limLow, limHigh){
  # Normalized, no mean centering scaling
  ggDF <- merge(genes, exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  row.names(ggDF) <- ggDF$x
  ggDF <- ggDF[ ,-1]
  ggDF <- t(ggDF)
  ggDF <- as.data.frame(ggDF)
  # Order by genes
  ggDF <- ggDF[ ,match(genes, colnames(ggDF))]
  # # Limits
  # ggDF[ggDF > limHigh] <- limHigh
  # ggDF[ggDF < limLow] <- limLow
  print(head(ggDF))
  # PAX6 STMN2 EOMES
  gg <- ggplot(ggDF, aes_string(x = colnames(ggDF)[1], y = colnames(ggDF)[2])) +
    geom_point(size = 0.5, alpha = 0.5, aes_string(color = colnames(ggDF)[3])) +
    # scale_color_distiller(name = "Normalized\nexpression", type = "div"
    #   , palette = 5, direction = -1) +
    scale_color_viridis()
  # # Limits
  # ggDF[ggDF > limHigh] <- limHigh
  # ggDF[ggDF < limLow] <- limLow
  # print(head(ggDF))
  # # PAX6 STMN2 EOMES
  # p1 <- ggplot(ggDF, aes_string(x = colnames(ggDF)[1], y = colnames(ggDF)[2])) +
  #   geom_point(size = 0.5, alpha = 0.5, aes_string(color = colnames(ggDF)[3])) +
  #   coord_cartesian(xlim = c(limLow, limHigh), ylim = c(limLow, limHigh)) +
  #   # scale_color_distiller(name = "Normalized\nexpression", type = "div"
  #   #   , palette = 5, direction = -1) +
  #   scale_color_viridis(limits = c(limLow, limHigh))
  return(gg)
}
################################################################################

### Violin plots

## Violin plots of expression by cluster
Gene_Expression_By_Cluster_ViolinPlot <- function(genes, exprM, clusterIDs
  , geneOrder = NULL, grouping = NULL){

  # Normalized, no mean centering scaling
  ggDF <- merge(genes, exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF <- melt(ggDF)
  # Add cluster ID
  idx <- match(ggDF$variable, names(clusterIDs))
  ggDF$CLUSTER <- clusterIDs[idx]
  # Set gene order if provided
  if (! is.null(geneOrder)) {
    ggDF$x <- factor(ggDF$x, levels = geneOrder)
  }
  print(str(ggDF))
  # Mean expression of gene group if grouping is provided
  if (! is.null(grouping)) {
    idx <- match(ggDF$x, genes)
    ggDF$Grouping <- grouping[idx]
    ggDF <- aggregate(
      ggDF$value, list(ggDF$Grouping, ggDF$variable, ggDF$CLUSTER)
      , mean, na.rm = TRUE)
    colnames(ggDF) <- c("x", "variable", "CLUSTER", "value")
  }

  # Violin plots of expression by cluster
  gg <- ggplot(ggDF, aes(x = CLUSTER, y = value)) +
    geom_violin(aes(fill = CLUSTER)) +
    geom_jitter(aes(x = CLUSTER, y = value)
      , size = 0.02, height = 0, alpha = 0.1) +
    facet_wrap(~x, scales = "free", ncol = 3) +
    theme(legend.position = "none") +
    ylab("Normalized expression") +
    xlab("Clusters")

  return(gg)
}

## Violin plots of expression faceted by cluster
Gene_Expression_Facet_By_Cluster_ViolinPlot <- function(genes, exprM, clusterIDs
  , geneOrder = NULL, grouping = NULL, ggtitle = NULL, ncol){

  # Normalized, no mean centering scaling
  ggDF <- merge(genes, exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  ggDF <- melt(ggDF)
  # Add cluster ID
  idx <- match(ggDF$variable, names(clusterIDs))
  ggDF$CLUSTER <- clusterIDs[idx]
  # Set gene order if provided
  if (! is.null(geneOrder)) {
    ggDF$x <- factor(ggDF$x, levels = geneOrder)
  }
  print(str(ggDF))
  # Mean expression of gene group if grouping is provided
  if (! is.null(grouping)) {
    idx <- match(ggDF$x, genes)
    ggDF$Grouping <- grouping[idx]
    ggDF <- aggregate(
      ggDF$value, list(ggDF$Grouping, ggDF$variable, ggDF$CLUSTER)
      , mean, na.rm = TRUE)
    colnames(ggDF) <- c("x", "variable", "CLUSTER", "value")
  }

  # Violin plots of expression by cluster
  gg <- ggplot(ggDF, aes(x = x, y = value)) +
    geom_violin(aes(fill = x)) +
    geom_jitter(aes(x = x, y = value)
      , size = 0.02, height = 0, alpha = 0.1) +
    facet_wrap(~CLUSTER, scales = "free", ncol = ncol) +
    theme(legend.position = "none") +
    ylab("Normalized expression") +
    xlab("Genes") +
    ggtitle(ggtitle)

  return(gg)
}
################################################################################
