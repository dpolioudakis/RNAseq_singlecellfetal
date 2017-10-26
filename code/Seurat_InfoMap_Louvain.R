# Damon Polioudakis
# 2017-03-27
# Compare InfoMap and Louvain clustering to Seurat clustering
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(igraph)
require(RANN)
require(reshape2)
require(Seurat)
require(Matrix)
require(gridExtra)
# require(fpc)

# source("GSE81904_class.R")

load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb
rm(fetb)

## Variables
graphCodeTitle <- "Seurat_InfoMap_Louvain.R"
outGraph <- "../analysis/graphs/Seurat_InfoMap_Louvain_"
outAnalysis <- "../analysis/Seurat_InfoMap_Louvain_"
################################################################################

## Function adopted from Shekhar et al. GSE81904_class.R
doGraph_clustering <- function(object,cells.use=NULL,pcs.use=1:10,num.nn=30
  , do.jaccard=FALSE, method="Louvain") {
  
  if (do.jaccard){
    weights=TRUE;
    method_print = paste0(method,"-","Jaccard")
  } else {
    weights=NULL;
    method_print = method
  }
  
  print(paste0("Performing ", method_print, " clustering. Using ", num.nn
    , " nearest neighbors, and ", max(pcs.use), " PCs"))
  
  if (is.null(cells.use)){
    data.use=object@pca.obj[[1]]$x[,pcs.use]
  } else {
    data.use=object@pca.obj[[1]]$x[cells.use,pcs.use]
  } 
  
  print(dim(data.use))
  
  print("Get edges...")
  Adj = get_edges(data.use,nn=num.nn,do.jaccard=do.jaccard)
  
  print("Graph adjacency...")
  g=graph.adjacency(Adj, mode = "undirected", weighted=weights)
  if (method=="Louvain") graph.out = cluster_louvain(g)
  if (method=="Infomap") graph.out = cluster_infomap(g)
  
  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
  print("Outputting clusters ..")
  # object@meta$clust = NULL
  # object@meta[names(clust.assign),"clust"]=clust.assign
  # object@group=clust.assign; names(object@group)=names(clust.assign);               
  # 
  # 
  # return(object) 
  return(clust.assign)
}

## Function adopted from Shekhar et al. GSE81904_class.R
get_edges=function(X,nn=30,do.jaccard=TRUE) {
  nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
  
  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1
  
  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  if (do.jaccard){
    
    NN = nearest$nn.idx
    jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges$C = jaccard_dist
    edges = subset(edges, C != 0)
    edges$C = edges$C/max(edges$C)
  }
  
  Adj = matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
  Adj[cbind(edges$A,edges$B)] = edges$C
  Adj[cbind(edges$B,edges$A)] = edges$C
  return(Adj)
}




#First lets stash our identities for later
seuratO <- StashIdent(seuratO, save.name = "SeuratCluster")


# Alternatively, we can also cluster the cells using Infomap (note that the
# do.jaccard option has been set to FALSE),
dsq.bip = doGraph_clustering(seuratO, pcs.use = 1:10, num.nn = 350
  , do.jaccard = FALSE, method = "Infomap")

# Add cluster IDs to Seurat object
seuratO@ident <- dsq.bip

## Plot
# tSNE plot - Seurat
plot1 <- TSNEPlot(seuratO, do.label = TRUE, group.by = "SeuratCluster"
  , pt.size = 0.2, do.return = TRUE, no.legend = TRUE)
plot1 <- plot1 + ggtitle(paste0(graphCodeTitle
  , "\n"
  , "\ntSNE plot, each point is a cell"
  , "\nColor indicates cluster assignment"
  , "\nSeurat clustering"
  , "\n"))
# tSNE plot - InfoMap
plot2 <- TSNEPlot(seuratO, do.label = TRUE, pt.size = 0.2, do.return = TRUE
  , no.legend = TRUE)
plot2 <- plot2 + ggtitle(paste0(graphCodeTitle
  , "\n"
  , "\ntSNE plot, each point is a cell"
  , "\nColor indicates cluster assignment"
  , "\nInfomap clustering"
  , "\n"))
# Save
png(paste0(outGraph, "Infomap.png"), width = 8, height = 6, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

## Jacard index

Jaccard_Index <- function(v1, v2) {
  sum(v1 %in% v2) / (length(v1) + length(v2) - sum(v1 %in% v2))
}

# Empty matrix
jiM <- matrix(NA
  , length(unique(seuratO@data.info$SeuratCluster))
  , length(unique(seuratO@ident)))
# Fill with Jaccard index
for (i in 1:length(unique(seuratO@data.info$SeuratCluster))){
  for (j in 1:length(unique(seuratO@ident))){
    v1 <- row.names(seuratO@data.info)[seuratO@data.info$SeuratCluster == i]
    v2 <- row.names(seuratO@data.info)[seuratO@ident == j]
    jiM[i,j] <- Jaccard_Index(v1, v2)
  }
}

# Distance
d <- dist(t(seuratO@scale.data))

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))

cs <- cluster.stats(d = d, clustering = seuratO@data.info$SeuratCluster)

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))


