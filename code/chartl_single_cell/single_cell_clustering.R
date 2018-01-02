# Rscript single_cell_clustering.R ExprMatrix.txt out --input_tpm INPUT_TPM --no_log NO_LOG

suppressPackageStartupMessages(library('scater'))
suppressPackageStartupMessages(library('SC3')) # bioconductor
library(argparse)
library(Rtsne)
library(fpc)
library(FNN)
library(cluster)
suppressPackageStartupMessages(library(dendextend))
library(mvnfast)
library(igraph)

options(stringsAsFactors=F)

getArgs <- function() {
  parser = ArgumentParser()
  parser$add_argument('expression', help='The (single-cell) expression matrix')
  parser$add_argument('out_dir', help='The output directory', default='.')
  parser$add_argument('--marker_list', help=paste0('The known marker list',
    '(these will be included regardless of expression levels, and used to',
    'annotate clusters. Format: gene(tab)cell.type'), default=NULL)
  parser$add_argument('--clust_inflate',
    help='Inflate the pamk cluster estimates by this amount',
    type='integer', default=2)
  parser$add_argument('--min_prop',
    help=paste0('Genes are included if # counts > min_count in at least',
      '`min_prop` proportion of cells'), type='double',
    default=0.05)
  parser$add_argument('--max_prop',
    help=paste0('Genes are included if # counts > min_count in at least',
      '`min_prop` proportion of cells'), type='double',
    default=NULL)
  parser$add_argument('--min_count', help=paste0('Genes are included if # counts > min_count in at least',
    '`min_prop` proportion of cells'), type='double',
    default=2)
  parser$add_argument('--impute', help='impute counts with scImpute',
    type='logical', default=FALSE)
  parser$add_argument('--ncores', help='Number of cores for imputation', default=1,
    type='integer')
  parser$add_argument('--out_prefix', help='File prefix for output (do NOT include directory)',
    default='')
  parser$add_argument('--tsne_iter', help='The tSNE iterations', default=1000, type='integer')
  parser$add_argument('--use_scimpute_clusters', help=paste0('Use scImpute to determine',
    'imputation clusters, rather than kNN+walktrap'), default=F, type='logical')
  parser$add_argument('--input_tpm', help='The input is TPM/FPKM. Do not normalize.', type='logical', default=F)
  parser$add_argument('--no_log', help='INput already logged. Do not do so again.', type='logical', default=F)
  parser$add_argument('--log_offset', help='The log-offset (log(THIS + expr))', type='double', default=1.01)
  parser$add_argument('--force_markers', help='Use the markers in the marker list for clustering, even if they fail filters', type='logical', default=F)
  parser$add_argument('--max_cluster_markers', help='If # of markers > this value, they will be ranked by CV and subset to this number', default=1200,
    type='integer')

  parser$parse_args()
}

dropOutliers <- function(expr.dat) {
  any.outlier <- T
  iter <- 1
  while ( any.outlier ) {
    print(sprintf('outlier iter %d...', iter))
    iter <- 1 + iter
    dat.pc <- svd(expr.dat, nu=15, nv=15)
    outliers <- rep(F, ncol(expr.dat))
    for ( pc in 1:15 ) {
      z.scores <- scale(dat.pc$v[,pc])
      outliers <- outliers | ( abs(z.scores) > 7 )
    }
    expr.dat <- expr.dat[, ! outliers]
    any.outlier <- any(outliers)
  }
  expr.dat
}

main <- function(args) {
  # read in the expression
  dat.counts <- read.table(args$expression, header=T, row.names=1)
  rn <- rownames(dat.counts)
  cn <- colnames(dat.counts)
  dat.counts <- as.matrix(dat.counts)
  rownames(dat.counts) <- rn
  colnames(dat.counts) <- cn
  print(sprintf('Initial dimension: %d x %d', dim(dat.counts)[1], dim(dat.counts)[2]))
  dat.counts <- dat.counts[apply(dat.counts == 0, 1, sum) < ncol(dat.counts),]
  print(sprintf('After removing unobserved genes: %d x %d', dim(dat.counts)[1], dim(dat.counts)[2]))
  cell.counts <- apply(dat.counts, 2, sum)
  if ( any(cell.counts == 0) ) {
    stop('Cells with no counts?')
  }
  if ( args$out_prefix != "" && ! grepl('\\.$', args$out_prefix) && ! grepl('_$', args$out_prefix)) {
    out.prefix <- sprintf("%s.", args$out_prefix)
  } else {
    out.prefix <- args$out_prefix
  }
  above.min.count <- apply(dat.counts > args$min_count, 1, sum)
  if ( is.null(args$max_prop) ) {
    max_prop <- 1 - args$min_prop
  } else {
    max_prop <- args$max_prop
  }
  keep.genes <- rownames(dat.counts)[above.min.count > (args$min_prop * ncol(dat.counts))]
  dropout.rate <- apply(dat.counts == 0, 1, sum)/ncol(dat.counts)
  keep.genes.dropout <- rownames(dat.counts)[dropout.rate < args$max_prop]
  keep.genes <- union(keep.genes, keep.genes.dropout)

  if ( ! is.null(args$marker_list) ) {
    marker.genes <- read.table(args$marker_list, header=T)
    present.genes <- subset(marker.genes, gene %in% rownames(dat.counts))
    if ( args$force_markers ) {
      keep.genes <- union(keep.genes, present.genes$gene)
    }
    marker.gene.list <- marker.genes$gene
  }

  if ( args$input_tpm ) {
    cell.rates <- dat.counts
  } else {
    cell.rates <- sweep(as.matrix(dat.counts), MARGIN=2, cell.counts/10^6, FUN='/')
  }
  if ( args$no_log ) {
    cell.rates <- exp(cell.rates) - args$log_offset
    cell.rates.log <- cell.rates
  } else {
    cell.rates.log <- log10(cell.rates + args$log_offset)
  }
  if ( length(keep.genes) > args$max_cluster_markers ) {
    print(sprintf('Ranking %d genes by CV and subsetting to %d', length(keep.genes), args$max_cluster_markers))
    means <- apply(cell.rates.log[keep.genes,], 1, mean)
    vars <- apply(cell.rates.log[keep.genes,], 1, var)
    cvs <- sqrt(vars)/means
    keep.genes <- keep.genes[order(cvs, decreasing=T)[1:args$max_cluster_markers]]
  }
  cell.rates.log.cluster <- cell.rates.log[keep.genes,]
  print('Pre-filtering:')
  print(dim(cell.rates.log.cluster))
  cell.rates.log.cluster <- dropOutliers(cell.rates.log.cluster)
  print('Initial dimensions:')
  print(dim(cell.rates.log.cluster))
  if ( ! all(is.finite(cell.rates.log.cluster)) ) {
    save(list=ls(), file='debug.Rdata')
    stop('NA values')
  }

  dat.cluster <- dat.counts[rownames(cell.rates.log.cluster), colnames(cell.rates.log.cluster)]
  rate.cluster <- cell.rates[rownames(cell.rates.log.cluster), colnames(cell.rates.log.cluster)]

  sc3.res <- SC3.cluster(dat.cluster, rate.cluster, cell.rates.log.cluster)
  k.guess <- sc3.res$k.est
  print(table(sc3.res$clusters))
  while ( min(table(sc3.res$clusters)) < 5 ) {
    print('Small cluster identified, decreasing k')
    sc3.res <- adjust.sc3.k(sc3.res, sc3.res$k.est-1)
    print(table(sc3.res$clusters))
  }

  n.clust <- length(unique(sc3.res$clusters))
  clust.colors <- rainbow(n.clust, alpha=0.6)[factor(sc3.res$clusters)]
  cluster.df <- data.frame('sample'=colnames(cell.rates.log.cluster), 'cluster'=sc3.res$clusters, 'color'=clust.colors)
  orig.sample.df <- data.frame('sample'=colnames(dat.counts))
  cluster.df <- merge(orig.sample.df, cluster.df, by='sample', all.x=T)
  rownames(cluster.df) <- cluster.df$sample
  cluster.df.sub <- cluster.df[colnames(cell.rates.log.cluster),]
  cluster.df <- cluster.df[colnames(cell.rates.log),]

  tsne.highdim <- Rtsne(t(cell.rates.log.cluster), dims=5, perplexity=100/sqrt(k.guess),
    theta=0.1, max_iter=args$tsne_iter, verbose=TRUE)
  tsne.lowdim <- Rtsne(t(cell.rates.log.cluster), dims=2, perplexity=90/sqrt(k.guess), theta=0.1,
    max_iter=args$tsne_iter, verbose=TRUE)

  pdf(sprintf('%s/%stsne_high_dim.pre_impute.pdf', args$out_dir, out.prefix))
  pairs(tsne.highdim$Y, pch=16, col=cluster.df.sub$color)
  dev.off()

  pdf(sprintf('%s/%stsne_low_dim.pre_impute.pdf', args$out_dir, out.prefix))
  plot(tsne.lowdim$Y, pch=16, xlab='TSNE.1', ylab='TSNE.2', col=cluster.df.sub$color)
  dev.off()

  if ( args$impute ) {
    if ( args$use_scimpute_clusters ) {
      sc.res <- scImpute::imputation_model8(count=cell.rates.log, labeled=F, point=log10(args$log_offset),
        drop_thre=0.0, Kcluster=length(unique(sc3.res$clusters)), ncores=args$ncores,
        out_dir=args$out_dir)
    } else {
      full.clusters <- as.character(cluster.df$cluster)
      full.clusters[is.na(full.clusters)] <- 'NA'
      sc.res <- scImpute::imputation_wlabel_model8(count=cell.rates.log, labeled=T, point=log10(args$log_offset),
        cell_labels=full.clusters, Kcluster=NULL,
        drop_thre=0.0, ncores=args$ncores, out_dir=args$out_dir)
    }

    cell.rates.imp <- sc.res$count_imp
    print(sprintf('N.outliers: %d', length(sc.res$outlier)))
    rownames(cell.rates.imp) <- rownames(cell.rates.log)
    colnames(cell.rates.imp) <- colnames(cell.rates.log)
    cell.rates.log <- cell.rates.imp
    cell.rates.log.cluster <- cell.rates.imp[keep.genes,]  # note: now full sample set
  }

  set.seed(0xABCDEF)
  tsne.highdim <- Rtsne(t(cell.rates.log.cluster), dims=5, perplexity=100/sqrt(k.guess),
    theta=0.1, max_iter=args$tsne_iter, verbose=TRUE)
  tsne.lowdim <- Rtsne(t(cell.rates.log.cluster), dims=2, perplexity=90/sqrt(k.guess), theta=0.1,
    max_iter=args$tsne_iter, verbose=TRUE)

  pdf(sprintf('%s/%stsne_high_dim.pdf', args$out_dir, out.prefix))
  pairs(tsne.highdim$Y, pch=16, col=cluster.df$color)
  dev.off()

  pdf(sprintf('%s/%stsne_low_dim.pdf', args$out_dir, out.prefix))
  plot(tsne.lowdim$Y, pch=16, xlab='TSNE.1', ylab='TSNE.2', col=cluster.df$color)
  dev.off()

  if ( ! is.null(args$marker_list) ) {
    colByGene <- function(gene.name) {
      cr <- colorRamp(c('grey60', 'red'), alpha=0.6)
      gvec <- uniscale(cell.rates.log[gene.name,])
      gvec[is.na(gvec)] <- 0
      cols <- cr(gvec)
      sapply(1:nrow(cols), function(i) {rgb2hex(cols[i,])})
    }
    pdf(sprintf('%s/%stsne.markers.pdf', args$out_dir, out.prefix))
    for ( row in 1:length(present.genes$gene) ) {
      marker = present.genes$gene[row]
      if ( marker %in% rownames(cell.rates.log) ) {
        type = present.genes$cell.type[row]
        plot(tsne.lowdim$Y, pch=16, xlab='TSNE.1', ylab='TSNE.2', col=colByGene(marker),
          main=sprintf('%s (%s)', marker, type))
      }
    }
    dev.off()
  }

  dat.svd <- svd(cell.rates.log.cluster, nu=n.clust + 1, nv=n.clust + 1)
  # (k x ngene) x (ngene x sample)
  proj.svd <- t(dat.svd$u) %*% cell.rates.log.cluster
  dat.hclust <- hclust(dist(t(proj.svd), 'euclidean'))  # hierarchy just on the MDS

  pdf(sprintf('%s/%s_hierarchichal.pdf', args$out_dir, out.prefix))
  hcd <- as.dendrogram(dat.hclust)
  labels_colors(hcd) <- cluster.df$color[order.dendrogram(hcd)]
  #labels(hcd) <- sc3.res$clusters[order.dendrogram(hcd)]
  labels(hcd) <- rep('.', length(order.dendrogram(hcd)))

  plot(hcd)
  for ( pc1 in 1:(nrow(proj.svd)-1) ) {
    pc2 <- 1 + pc1
    plot(proj.svd[pc1,], proj.svd[pc2,], col=cluster.df$color, pch=16,
      xlab=sprintf('Cluster matrix PC%d', pc1), ylab=sprintf('Cluster matrix PC%d', pc2))
  }
  dev.off()
  clust.assign <- cluster.df$cluster
  if ( ! is.null(args$marker_list) ) {
    save(list=ls(), file='scpipe.annot.debug.Rda')
    annot.res <- annotateClusters(clust.assign, present.genes, cell.rates.log)
    clust.annot <- annot.res[[1]]
    marker.expr <- annot.res[[2]]
    write.table(clust.annot, file=sprintf('%s/%scluster_annotations.txt', args$out_dir, out.prefix), sep='\t', quote=F)
    pdf(sprintf('%s/%stsne.markers.annot.pdf', args$out_dir, out.prefix))
    for ( annot.name in colnames(marker.expr) ) {
      marker <- marker.expr[,annot.name]
      scaleColor <- function(zz) {
        cr <- colorRamp(c('grey60', 'red'), alpha=0.6)
        gvec <- uniscale(zz)
        gvec[is.na(gvec)] <- 0
        cols <- cr(gvec)
        sapply(1:nrow(cols), function(i) {rgb2hex(cols[i,])})
      }
      plot(tsne.lowdim$Y, pch=16, xlab='TSNE1', ylab='TSNE.2', col=scaleColor(marker),
        main=annot.name)
    }
    dev.off()
  }
  clust.df <- data.frame(cell=colnames(cell.rates.log), clust.id=clust.assign, dendro.rank=order(dat.hclust$order))
  tsne.df <- as.data.frame(cbind(tsne.highdim$Y, tsne.lowdim$Y))
  colnames(tsne.df) <- c(sprintf('TSNE.HD.%d', 1:5), sprintf('TSNE.LD.%d', 1:2))
  clust.df <- cbind(clust.df, tsne.df)

  write.table(clust.df, file=sprintf('%s/%scluster_definitions.txt', args$out_dir, out.prefix), sep='\t', quote=F)

  write.table(cell.rates.log, file=sprintf('%s/%sexpression_rate_matrix.txt', args$out_dir, out.prefix), sep='\t', quote=F)

}

uniscale <- function(x) {
  z <- x - min(x)
  z/max(z)
}

rgb2hex <- function(x) {
  sprintf("%s66", rgb(x[1], x[2], x[3], maxColorValue=255))
}


p.extreme.sim <- function(mean.vec, covar.mat, n.sims = 10^5, b.size=10^4, group=NULL, group.name='group') {
  # calculate marginal probabilities that each dimension is
  # the "largest" or "smallest" dimension.
  # if `group` is specified, it is the index of a set of variables to test
  # for *ALL* being greater than/less than the others, i.e.
  #   P[X1 > X3 & X1 > X4 & X2 > X3 & X2 > X4]  --> group=c(1,2)
  ## inputs:
  ##  mean.vec -- the beta values -- obtained from summary(lm)$coefficients[,1]
  ##  covar.mat -- the beta covariance -- obtained from summary(lm)$cov.unscaled
  ##  n.sums  -- number of simulations (for monte carlo integration)
  ##  b.size -- the batch size for monte carlo
  ##  group -- for testing whether coefficients 1, 2, 5 > 3, 4; specifiy c(1,2,5)
  ##  group.name -- a name for the group (e.g. 'CORTEX')
  sims <- 0
  if ( is.null(group) ) {
    lg.counts <- sm.counts <- 0 * mean.vec + 1
    while ( sims < n.sims ) {
      dat <- rmvn(b.size, mean.vec, covar.mat, ncores=1)
      lrg.sum <- rowSums(apply(dat, 1, function(x) { x == max(x) }))
      sml.sum <- rowSums(apply(dat, 1, function(x) { x == min(x) }))
      lg.counts <- lrg.sum + lg.counts
      sm.counts <- sml.sum + sm.counts
      sims <- sims + b.size
    }
    probs <- cbind(lg.counts, sm.counts)
    probs <- probs/(sims + length(mean.vec))  # adjust for the pseudocounts
    beta.diff <- sapply(1:length(mean.vec), function(i) { mean.vec[i] - mean(mean.vec[-i])})
    probs <- cbind(probs, beta.diff)
    colnames(probs) <- c('largest', 'smallest', 'delta_beta')
    rownames(probs) <- names(mean.vec)
  } else {
    l.count <- s.count <- 1
    while ( sims < n.sims ) {
      dat <- rmvn(b.size, mean.vec, covar.mat, ncores=1)
      lg <- sum(sapply(1:b.size, function(i) { all(dat[i,group] > max(dat[i, -group])) }))
      sm <- sum(sapply(1:b.size, function(i) { all(dat[i,group] < min(dat[i, -group])) }))
      l.count <- l.count + lg
      s.count <- s.count + sm
      sims <- sims + b.size
    }
    beta.diff <- mean(mean.vec[group]) - mean(mean.vec[-group])
    probs <- data.frame(largest=l.count/(2 + sims), smallest=s.count/(2 + sims), delta_beta=beta.diff)
    rownames(probs) <- group.name
  }
  as.data.frame(probs)
}

p.extreme <- function(expr.dat, cluster.assign) {
  # expr.dat        - matrix; (n.cell x n.markers)
  # cluster.assign  - vector; (n.cell x 1)
  expr.dat <- as.matrix(expr.dat)  # weird...
  if ( ncol(expr.dat) < 1 ) {
    print('No markers? Check the gene symbols match...')
    return(data.frame(cluster.assign))
  }
  expr.pv <- lapply(1:ncol(expr.dat), function(i) {
    mkr <- expr.dat[,i]
    my.df <- data.frame(expr=mkr, cluster=cluster.assign)
    lm.res <- lm(expr ~ factor(cluster) - 1, data=my.df)
    coefs <- summary(lm.res)$coefficients[,1]
    covars <- summary(lm.res)$cov.unscaled
    prob.res <- p.extreme.sim(coefs, covars)
    prob.res$marker <- colnames(expr.dat)[i]
    prob.res$group <- names(coefs)
    prob.res
  })
  do.call(rbind, expr.pv)
}

SC3.cluster <- function(raw.counts, rates, log.rates) {
  if ( ! all(dim(raw.counts) == dim(rates)) ) {
    stop('raw.counts dim does not match rates')
  }
  if ( ! all(dim(rates) == dim(log.rates)) ) {
    stop('rates dim does not match log rates dim')
  }
  if ( ! all(rownames(raw.counts) == rownames(log.rates)) ) {
    stop('gene names differ in input')
  }
  sce <- SingleCellExperiment(assays = list(counts=as.matrix(raw.counts)))
  normcounts(sce) <- rates
  logcounts(sce) <- log.rates
  exprs(sce) <- logcounts(sce)
  rowData(sce)$feature_symbol <- rownames(raw.counts)
  isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
  print('Consensus clustering with SC3')
  sc3.res <- sc3_prepare(sce, gene_filter=FALSE, svm_max=1e9, kmeans_nstart=100)
  # sc3.res <- sc3_prepare(sce, gene_filter=FALSE, ks=5:30, svm_max=1e9, kmeans_nstart=100)
  sc3.res <- sc3_estimate_k(sc3.res)
  sc3.k <- metadata(sc3.res)$sc3$k_estimation
  sc3.res <- sc3_calc_dists(sc3.res)
  sc3.res <- sc3_calc_transfs(sc3.res)
  sc3.res <- sc3_kmeans(sc3.res, ks=5:30)
  # sc3.res <- sc3_kmeans(sc3.res)
  sc3.res <- sc3_calc_consens(sc3.res)
  print('Done!')
  save(sc3.res, sce, file='test.Rdata')
  sc3.k.char <- as.character(sc3.k)
  list(sc3obj=sc3.res, clusters=metadata(sc3.res)$sc3$consensus[[sc3.k.char]]$silhouette[,'cluster'], k.est=sc3.k)
}

adjust.sc3.k <- function(sc3, new.k) {
  new.k.char <- as.character(new.k)
  list(sc3obj=sc3[['sc3obj']], clusters=metadata(sc3[['sc3obj']])$sc3$consensus[[new.k.char]]$silhouette[,'cluster'], k.est=new.k)
}

annotateClusters <- function(labels, marker.df, dat.expr) {
  # marker.df - data frame of (gene, cell.type) -- these should
  #             be up-regulated in a given cell type (and
  #             presumably only in that cell type)
  # labels    - labels discovered by the clustering
  # dat.expr  - the (log-rate normalized) expression matrix
  print('Annotating clusters')
  type.list <- unique(marker.df$cell.type)
  print(sprintf('Labels: %d (%d unique); expr: %d x %d', length(marker.df$cell.type), length(type.list), dim(dat.expr)[1], dim(dat.expr)[2]))
  type.mat <- matrix(0, nrow=ncol(dat.expr), ncol=length(type.list))
  i <- 1
  for ( ctype in type.list ) {
    marker.sub <- subset(marker.df, cell.type == ctype)
    tot.expr <- apply(dat.expr[marker.sub$gene,,drop=F], 2, sum)
    type.mat[,i] <- tot.expr
    i <- i + 1
  }
  rownames(type.mat) <- colnames(dat.expr)
  colnames(type.mat) <- type.list
  annot <- p.extreme(type.mat, labels)
  list(annot, type.mat)
}

if ( ! interactive() ) {
  args <- getArgs()
  main(args)
}
