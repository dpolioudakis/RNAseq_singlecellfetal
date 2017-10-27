# Damon Polioudakis
# Function library

################################################################################

### Differential expression

# Filter expression matrix by:
# percent of cells in cluster gene is expressed in
# fold change of gene in cluster versus all other cells
DE_Filters_ExpMatrix <- function(
  so
  , minPercent = NULL
  , foldChange = NULL
  , clusterID = NULL
  , cellID = NULL) {
  
  if (! is.null(minPercent)) {
    # Expressed > 0 counts in > X% of cells in cluster
    if (! is.null(clusterID)) {
      # Subset expression matrix to cluster
      cdf <- as.matrix(so@data)[ ,so@ident == clusterID]  
    }
    if (! is.null(cellID)) {
      # Subset expression matrix to cluster
      cdf <- as.matrix(so@data)[ ,colnames(so@data) %in% cellID]  
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
    if (! is.null(clusterID)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,so@ident == clusterID]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! so@ident == clusterID]
    }
    if (! is.null(cellID)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,colnames(so@data) %in% cellID]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! colnames(so@data) %in% cellID]
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
  exDF <- as.matrix(so@data[idxp & idxf, ])
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
DE_Linear_Model <- function (exDatDF, termsDF, mod) {
  
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
Format_DE <- function (deLM, so, clusterID) {
  # Combine log2 fold changes, p-values
  deDF <- data.frame(GENE = row.names(deLM$coefmat)
    , LOG_FC = deLM$coefmat[ ,2]
    , PVALUE = deLM$pvalmat[ ,2])
  # Order by pvalue
  deDF <- deDF[order(deDF$PVALUE), ]
  # Add cluster ID
  deDF$CLUSTER <- clusterID
  # Percent of cells in cluster expressing gene > 0 counts
  cdf <- as.matrix(so@data)[row.names(so@data) %in% deDF$GENE, so@ident == clusterID]
  deDF$PERCENT_CLUSTER <- (rowSums(cdf > 0) / ncol(cdf)) * 100
  # Percent of all cells expressing gene > 0 counts
  deDF$PERCENT_ALL <- (rowSums(as.matrix(so@data)[row.names(so@data) %in% deDF$GENE, ] > 0)
    / ncol(so@data[row.names(so@data) %in% deDF$GENE, ])) * 100
  # Order by log fold change
  deDF <- deDF[order(-deDF$LOG_FC), ]
  return(deDF)
}

################################################################################

### Expression heatmaps

# Expression heatmap, cells ordered by cluster
Heatmap_By_Cluster <- function(
  geneGroupDF, exprM, seuratO, clusters, lowerLimit, upperLimit, geneOrder = NULL) {
  
  # Subset expression matrix to genes of interest by merging
  ggDF <- merge(geneGroupDF[c("GENE", "GROUP")], exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  # Set group factor levels
  ggDF$GROUP <- as.factor(ggDF$GROUP)
  # colnames(ggDF)[1:2] <- c("GENE", "GROUP")
  # Remove blanks
  ggDF <- ggDF[! ggDF$GENE == "", ]
  # Save order to set levels later
  if (is.null(geneOrder)) {
    levels <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  } else {
    levels <- rev(paste0(geneOrder, "   ", ggDF$GROUP))
  }
  ggDF <- melt(ggDF)
  # Add seurat clusters
  idx <- match(ggDF$variable, names(seuratO@ident))
  ggDF$SEURAT_CLUSTERS <- seuratO@ident[idx]
  # Subset clusters
  ggDF <- ggDF[ggDF$SEURAT_CLUSTERS %in% clusters, ]
  # Add group to gene name and order by levels
  ggDF$GENE_GROUP <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  ggDF$GENE_GROUP <- factor(ggDF$GENE_GROUP, levels = levels)
  # Set limits
  ggDF$value[ggDF$value < lowerLimit] <- lowerLimit
  ggDF$value[ggDF$value > upperLimit] <- upperLimit
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = GENE_GROUP, fill = value)) +
    geom_tile() +
    facet_grid(GROUP~SEURAT_CLUSTERS, space = "free", scales = "free") +
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
Heatmaps_By_Cluster_Combined <- function(geneGroupDF, exprM, seuratO
  , clusters1, clusters2, clusters3, lowerLimit, upperLimit, geneOrder = NULL) {
  
  p1 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , seuratO = seuratO, clusters = clusters1
    , lowerLimit = lowerLimit, upperLimit = upperLimit
    , geneOrder
  )
  p1 <- p1 + theme(
    axis.title.x = element_blank()
    , strip.text.y = element_blank()
    , legend.position = "none"
  )
  p2 <- Heatmap_By_Cluster(geneGroupDF = geneGroupDF, exprM = exprM
    , seuratO = seuratO, clusters = clusters2
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
    , seuratO = seuratO, clusters = clusters3
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
    # facet_grid(GROUP~SEURAT_CLUSTERS, space = "free", scales = "free") +
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
FeaturePlot_Graph <- function(tsneDF, title, limLow, limHigh) {
  ggFp <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = 0.05) +
    # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
    #   , high = "#e31a1c", limits = c(0, 2)) +
    scale_color_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1) +
    # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
    #   , high = "red", limits = c(limLow, limHigh)) +
    ggtitle(title)
  return(ggFp)
}

# Color tSNE plot by expression from Mean_Expression()
FeaturePlot_Graph_CentScale <- function(tsneDF, title, limLow, limHigh) {
  ggFp <- ggplot(tsneDF, aes(x = tSNE_1, y = tSNE_2, col = EXPRESSION)) +
    geom_point(size = 0.05) +
    # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
    #   , high = "#e31a1c", limits = c(0, 2)) +
    scale_color_distiller(name = "Normalized\nexpression\nz-score", type = "div"
      , palette = 5, direction = -1, limits = c(limLow, limHigh)) +
    # scale_colour_gradient(name = "Normalized\nExpression", low = "#a6cee3"
    #   , high = "red", limits = c(limLow, limHigh)) +
    ggtitle(title)
  return(ggFp)
}

# tSNE plot colored by Seurat clustering
TSNE_Plot <- function(seuratO) {
  # tSNE graph colored by cluster
  ggTsne <- TSNEPlot(seuratO, do.label = TRUE, pt.size = 0.1, do.return = TRUE
    , no.legend = FALSE)
  ggTsne <- ggTsne + ggtitle(paste0(
    "tSNE plot, each point is a cell"
    , "\nColor indicates cluster assignment"
  ))
  ggTsne <- ggTsne + theme_set(theme_bw()) +
    theme_set(theme_get() + theme(text = element_text(size = 16))) +
    theme_update(plot.title = element_text(size = 12)) +
    theme_update(axis.line = element_line(colour = "black")
      , plot.background = element_blank() 
      , panel.border = element_blank()
    )
  return(ggTsne)
}
# ggL <- append(list(ggTsne), ggL)

# Wrapper function for feature plots
FeaturePlot_CentScale <- function(genes, tsneDF, seuratO, limLow, limHigh) {
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
  ggTsne <- TSNE_Plot(seuratO) + theme(legend.position = "none")
  ggL <- append(list(ggTsne), ggL)
  return(ggL)
}
# Wrapper function for feature plots
FeaturePlot <- function(genes, tsneDF, seuratO, exM, limLow, limHigh) {
  # Collect tSNE values for ggplot
  tsneDF <- as.data.frame(seuratO@dr$tsne@cell.embeddings)
  # Loop through and plot each group of genes
  ggL <- lapply(genes, function(gene) {
    print(gene)
    tsneDF <- Mean_Expression(tsneDF, gene, exM)
    tsneDF <- Set_Limits(tsneDF, limLow = limLow, limHigh = limHigh)
    ggFp <- FeaturePlot_Graph(tsneDF, title = paste0("\n", gene)
      , limLow = limLow, limHigh = limHigh)
    return(ggFp)
  })
  ggTsne <- TSNE_Plot(seuratO) + theme(legend.position = "none")
  ggL <- append(list(ggTsne), ggL)
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

# Cowplot plot_grid and add title
Plot_Grid <- function(ggPlotsL, ncol, title, rel_height, ...) {
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
    , by.y = "CELLS", all.x = TRUE)
  
  # Order cells like monocle object
  idx <- match(row.names(pData(monocleO)), m$CELL)
  m <- m[idx, ]
  
  # Add column with expression values to monocle object pData slot
  pData(monocleO)$EXPRESSION = m$EXPRESSION
  return(monocleO)
}

# Color monocle trajectory by expression of gene list
Plot_Trajectory_Gene_Expression <- function (monocleO, genes, exprM
  , limHigh, limLow, title) {
  
  mo <- Add_Expr_To_MonocleObj(monocleO = monocleO, genes = genes
    , exprM = exprM, limHigh, limLow)
  
  gg <- plot_cell_trajectory(mo, 1, 2
    , color_by = "EXPRESSION", cell_size = 0.01)
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
  ggDF <- merge(geneGroupDF[c("GENE", "GROUP")], exprM
    , by.x = 1, by.y = "row.names", all.x = TRUE)
  # Set group factor levels
  ggDF$GROUP <- as.factor(ggDF$GROUP)
  # colnames(ggDF)[1:2] <- c("GENE", "GROUP")
  # Remove blanks
  ggDF <- ggDF[! ggDF$GENE == "", ]
  # Save order to set levels later
  if (is.null(geneOrder)) {
    levels <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  } else {
    levels <- rev(paste0(geneOrder, "   ", ggDF$GROUP))
  }
  if (centScale == TRUE) {
    GENE <- ggDF$GENE
    GROUP <- ggDF$GROUP
    ggDF <- as.data.frame(t(scale(t(ggDF[ ,-c(1:2)]))))
    ggDF$GENE <- GENE
    ggDF$GROUP <- GROUP
  }
  ggDF <- melt(ggDF)
  # Add seurat clusters
  idx <- match(ggDF$variable, pData(monocleO)$CELL)
  ggDF$MONOCLE_STATE <- pData(monocleO)$State[idx]
  # Subset clusters
  ggDF <- ggDF[ggDF$MONOCLE_STATE %in% clusters, ]
  # Add group to gene name and order by levels
  ggDF$GENE_GROUP <- paste0(ggDF$GENE, "   ", ggDF$GROUP)
  ggDF$GENE_GROUP <- factor(ggDF$GENE_GROUP, levels = levels)
  # Set limits
  ggDF$value[ggDF$value < lowerLimit] <- lowerLimit
  ggDF$value[ggDF$value > upperLimit] <- upperLimit
  # ggplot
  gg <- ggplot(ggDF, aes(x = variable, y = GENE_GROUP, fill = value)) +
    geom_tile() +
    facet_grid(GROUP~MONOCLE_STATE, space = "free", scales = "free") +
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