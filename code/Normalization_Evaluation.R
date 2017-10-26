# Damon Polioudakis
# 2017-06-22
# Evaluate single-cell normalization methods by comparing to bulk CP/GZ DE
################################################################################

rm(list = ls())

require(NODES)
require(scran)
require(scater)
require(Seurat)
require(ggplot2)
require(Matrix)
require(NODES)
require(monocle)
require(R.matlab)
require(data.table)
require(fdrtool)
require(reshape2)
require(GGally)

load("../analysis/Normalization_Evaluation_Workspace.RData")

# Luis RNAseq fetal brain VZ CP
# CQN normalized (depth, length, GC)
load("../bulk_VZ_CP_from_ATAC/analysis/Expression_CQN_Ft0_Luis_RNAseq_VZCP.RData")
# Luis RNAseq fetal brain VZ CP
buExDF <- read.csv(
  "../bulk_VZ_CP_from_ATAC/data/htseq/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics - Luis RNAseq - Seq PC1-5
load("../bulk_VZ_CP_from_ATAC/analysis/SeqPC_Luis_RNAseq_VZCP.RData")

# Luis RNAseq
buMtDF <- read.csv(
  "../bulk_VZ_CP_from_ATAC/metadata/VZCP_sampleinfo.csv", header = TRUE)

# Drop-seq
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")

# MAGIC imputation
# Katherine - no read depth normalization
dir <- "/u/project/eeskin/geschwind/ksheu/Rotation3/magic/"
mgExDF <- fread(paste0(dir, "data_imputed_k2_t6.txt"))
mgExDF = t(mgExDF)
colnames = readMat(paste0(dir, "cells.mat"))
rownames = readMat(paste0(dir, "gene_names.mat"))
colnames(mgExDF) <- unlist(colnames)
row.names(mgExDF) <- unlist(rownames)
# Damon - read depth normalization
mgnExDF <- read.csv(
  "../analysis/Expression_Matrix_Compile_dge_FtMm250_MagicNorm_DS-2-3-4.csv"
  , header = TRUE, row.names = 1)
# Damon - no read depth normalization
mgdExDF <- read.csv(
  "../analysis/Expression_Matrix_Compile_dge_FtMm250_Magic_DS-2-3-4.csv"
  , header = TRUE, row.names = 1)
mgdExDF <- t(mgdExDF)

# BiomaRt
bmDF <- read.csv(
  "../source/BiomaRt_Compile_GeneInfo_GRCh37_Ensembl75.csv", header = TRUE)

outGraph <- "../analysis/graphs/Normalization_Evaluation_"
outAnalysis <- "../analysis/Normalization_Evaluation_"
################################################################################

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
  deCoefPvalLM
}
################################################################################

### Format

exDF <- fetb@raw.data
metDF <- fetb@data.info
metDF$Region <- "GZ"
metDF$Region[grep("CP", row.names(metDF))] <- "CP"
metDF$Cell <- row.names(metDF)
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

# Cleanup ENSEMBL IDs
row.names(vzcpCqnDatDF) <- gsub("\\..*", "", row.names(vzcpCqnDatDF))
row.names(buExDF) <- gsub("\\..*", "", row.names(buExDF))

# Remove ERCCs and STAR stats from tail
# Bulk
tail(buExDF, 10)
buExDF <- head(buExDF, -97)
tail(buExDF, 5)
################################################################################

### Comparing Seurat regression to linear model

noCentSO <- CreateSeuratObject(raw.data = fetb@raw.data, min.cells = 0
  , min.genes = 0
  , normalization.method = "LogNormalize", scale.factor = 10000
  , project = "DS-2-3-4", do.scale = FALSE, do.center = FALSE)

mito.genes <- grep("^MT-", rownames(noCentSO@data), value = TRUE)
percent.mito <- apply(
  expm1(noCentSO@data[mito.genes, ]), 2, sum)/apply(expm1(noCentSO@data), 2, sum)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC
# stats
noCentSO <- AddMetaData(noCentSO, percent.mito, "percent.mito")

noCentSO <- ScaleData(noCentSO, vars.to.regress = c("nUMI", "percent.mito")
  , do.scale = FALSE, do.center = FALSE)

# Make dataframe of covariates to regress out
covDF <- data.frame(nUMI = noCentSO@meta.data$nUMI)

# Regress out confounding variables
# Colnames show indexes to use to subset matrix X
RegressCovariates <- function (exM, covDF) {
  X = model.matrix(~ ., data = covDF)
  Y = exM
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  b = as.data.frame(t(beta))
  # Colnames show indexes to use to subset matrix X
  print(colnames(X))
  to_regress = (as.matrix(X[,3]) %*% (as.matrix(beta[3,])))
  exRegCovM = exM - t(to_regress)
  return(exRegCovM)
}

lm <- RegressCovariates(as.matrix(noCentSO@data), covDF)
################################################################################

### Luis RNAseq

## CPM

# All cells
rDep <- (apply(buExDF, 2, sum) / 10^6)
buCpmDF <- buExDF / rDep
buCpmDF <- log(buCpmDF + 0.01, 2)

## CPM DE

# Terms:
# RIN, SeqPC1, SeqPC2
# termsDF <- data.frame(buMtDF[ ,c("ExpCondition", "RIN.y")], topPCdatSeq[ ,1:2])
termsDF <- data.frame(ExpCondition = buMtDF$ExpCondition)

# DE
buCpmDeLM <- DE_Linear_Model(
  # buCpmDF, termsDF, "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2")
  buCpmDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
buCpmDeDF <- data.frame(ENSEMBL = row.names(buCpmDeLM$coefmat)
  , BULK_CPM_GZCP_LOG2_FC = buCpmDeLM$coefmat[ ,2]
  , BULK_CPM_GZCP_PVAL = buCpmDeLM$pvalmat[ ,2])

## FDR
df <- buCpmDeDF[! is.na(buCpmDeDF$BULK_CPM_GZCP_PVAL), ]
corrected <- fdrtool(df$BULK_CPM_GZCP_PVAL, statistic = "pvalue", plot = FALSE)
df$BULK_CPM_FDR <- corrected$lfdr
buCpmDeDF <- merge(buCpmDeDF, df[ ,c("ENSEMBL", "BULK_CPM_FDR")]
  , by = "ENSEMBL", all = TRUE)

# (CQN outputs log2 transformed values, these are input into linear model, and
#  then beta values output from linear mode are the log2 fold change)

## DE cqn

# Terms:
# RIN, SeqPC1, SeqPC2
# termsDF <- data.frame(buMtDF[ ,c("ExpCondition", "RIN.y")], topPCdatSeq[ ,1:2])
termsDF <- data.frame(ExpCondition = buMtDF$ExpCondition)

# DE
buDeLM <- DE_Linear_Model(
  # vzcpCqnDatDF, termsDF, "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2")
  vzcpCqnDatDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
buDeDF <- data.frame(ENSEMBL = row.names(buDeLM$coefmat)
  , BULK_CQN_GZCP_LOG2_FC = buDeLM$coefmat[ ,2]
  , BULK_CQN_GZCP_PVAL = buDeLM$pvalmat[ ,2])

## FDR
corrected <- fdrtool(buDeDF$BULK_CQN_GZCP_PVAL, statistic = "pvalue", plot = FALSE)
buDeDF$BULK_CQN_FDR <- corrected$lfdr
################################################################################

### Scater

# The counts are then combined into a single matrix for constructing a SCESet
# object. For convenience, metadata for all cells are stored in the same object
# for later access.
pd <- new("AnnotatedDataFrame", data = metDF)
rownames(pd) <- pd$Cell
sce <- newSCESet(countData = exDF, phenoData = pd)

# Normalization of cell-specific biases is performed using the deconvolution
# method in the computeSumFactors function. Here, we cluster similar cells
# together and normalize the cells in each cluster using the deconvolution
# method. This improves normalization accuracy by reducing the number of DE
# genes between cells in the same cluster. Scaling is then performed to ensure
# that size factors of cells in different clusters are comparable.
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
summary(sizeFactors(sce))

# size factors are estimated based on median ratios and are more robust to the
# presence of DE between cells
pdf(paste0(outGraph, "QC_TotalCountsVsSizeFactors.pdf"))
plot(sizeFactors(sce), sce$nUMI, log="xy",
  ylab="nUMI", xlab="Size factor")
dev.off()

# Normalized log-expression values are computed for each endogenous gene or
# spike-in transcript using the appropriate size factors.
sce <- normalize(sce)

v1 <- rowMeans(as.matrix(seuratO@data))
v2 <- rowMeans(seuratO@scale.data)
df1 <- data.frame(data = v1)
df2 <- data.frame(scale.data = v2)
df3 <- merge(df1, df2, by = "row.names")
ggplot(df, aes(x = data, y = scale.data)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 10), ylim = c(-10e-16, 10e-16))
ggsave(paste0(outGraph, "SeuratScaleVsSeuratData.pdf"))

v1 <- rowMeans(as.matrix(fetb@data))
v2 <- rowMeans(fetb@scale.data)
df1 <- data.frame(data = v1)
df2 <- data.frame(scale.data = v2)
df3 <- merge(df1, df2, by = "row.names")
ggplot(df, aes(x = data, y = scale.data)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 10), ylim = c(-10e-16, 10e-16))
ggsave(paste0(outGraph, "SeuratScaleVsSeuratData.pdf"))

v1 <- rowSums(norm_exprs(sce)[ ,sce$Region == "CP"])
v2 <- rowSums(norm_exprs(sce)[ ,sce$Region == "GZ"])
length(v1)
length(v2)
tail(v1)
tail(v2)
df1 <- data.frame(scater = v1/v2)
v1 <- rowSums(as.matrix(seuratO@data)[ ,metDF$Region == "CP"])
v2 <- rowSums(as.matrix(seuratO@data)[ ,metDF$Region == "GZ"])
length(v1)
length(v2)
tail(v1)
tail(v2)
df2 <- data.frame(seurat = v1/v2)
df3 <- merge(df1, df2, by = "row.names")
ggplot(df3, aes(x = scater, y = seurat)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))
ggsave(paste0(outGraph, "SeuratVsScater.pdf"))

## DE via LM for scater

# Terms: CP/GZ
termsDF <- data.frame(ExpCondition = sce$Region)

# DE
sceDeLM <- DE_Linear_Model(norm_exprs(sce), termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
sceDeDF <- data.frame(GENE = row.names(sceDeLM$coefmat)
  , SCE_GZCP_LOG2_FC = sceDeLM$coefmat[ ,2]
  , SCE_GZCP_PVAL = sceDeLM$pvalmat[ ,2])
# Remove NaNs
sceDeDF <- sceDeDF[! is.na(sceDeDF$SCE_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(sceDeDF$SCE_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
sceDeDF$SCE_FDR <- corrected$lfdr

# Add ensembl
sceDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(sceDeDF$GENE, bmDF$hgnc_symbol)]

## DE via LM for scater + regress out nUMI

# Terms: CP/GZ
termsDF <- data.frame(ExpCondition = sce$Region, nUMI = metDF$nUMI)

# DE
sceuDeLM <- DE_Linear_Model(norm_exprs(sce), termsDF, "y~ExpCondition+nUMI")

# Combine log2 fold changes, p-values
sceuDeDF <- data.frame(GENE = row.names(sceuDeLM$coefmat)
  , SCEnUMI_GZCP_LOG2_FC = sceuDeLM$coefmat[ ,2]
  , SCEnUMI_GZCP_PVAL = sceuDeLM$pvalmat[ ,2])
# Remove NaNs
sceuDeDF <- sceuDeDF[! is.na(sceuDeDF$SCEnUMI_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(sceuDeDF$SCEnUMI_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
sceuDeDF$SCEnUMI_FDR <- corrected$lfdr

# Add ensembl
sceuDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(sceuDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### Prabhakar pQ

exDF <- as.matrix(fetb@raw.data)
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]
pqDF <- pQ(exDF, throw_sd = 0, hard_outlier = 200)

# Terms: CP/GZ
termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
pqDeLM <- DE_Linear_Model(pqDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
pqDeDF <- data.frame(GENE = row.names(pqDeLM$coefmat)
  , PQ_GZCP_LOG2_FC = pqDeLM$coefmat[ ,2]
  , PQ_GZCP_PVAL = pqDeLM$pvalmat[ ,2])
# Remove NaNs
pqDeDF <- pqDeDF[! is.na(pqDeDF$PQ_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(pqDeDF$PQ_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
pqDeDF$PQ_FDR <- corrected$lfdr

# Add ensembl
pqDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(pqDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### DE via LM for Seurat CPM

exDF <- as.matrix(fetb@data)
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

# Terms: CP/GZ
termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
cpmDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
cpmDeDF <- data.frame(GENE = row.names(cpmDeLM$coefmat)
  , CPM_GZCP_LOG2_FC = cpmDeLM$coefmat[ ,2]
  , CPM_GZCP_PVAL = cpmDeLM$pvalmat[ ,2])
# Remove NaNs
cpmDeDF <- cpmDeDF[! is.na(cpmDeDF$CPM_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(cpmDeDF$CPM_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
cpmDeDF$CPM_FDR <- corrected$lfdr

# Add ensembl
cpmDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(cpmDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### DE via LM for raw counts

exDF <- fetb@raw.data
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

# Terms: CP/GZ
termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
rawDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
rawDeDF <- data.frame(GENE = row.names(rawDeLM$coefmat)
  , RAW_GZCP_LOG2_FC = rawDeLM$coefmat[ ,2]
  , RAW_GZCP_PVAL = rawDeLM$pvalmat[ ,2])
# Remove NaNs
rawDeDF <- rawDeDF[! is.na(rawDeDF$RAW_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(rawDeDF$RAW_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
rawDeDF$RAW_FDR <- corrected$lfdr

# Add ensembl
rawDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(rawDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### DE via LM nUMI pctMT on CPM like Seurat

exDF <- as.matrix(noCentSO@data)
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

termsDF <- data.frame(ExpCondition = metDF$Region, nUMI = metDF$nUMI
  , pctMT = metDF$percent.mito)

# DE
srtDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition+nUMI+pctMT")

# Combine log2 fold changes, p-values
srtDeDF <- data.frame(GENE = row.names(srtDeLM$coefmat)
  , SRT_GZCP_LOG2_FC = srtDeLM$coefmat[ ,2]
  , SRT_GZCP_PVAL = srtDeLM$pvalmat[ ,2])
# Remove NaNs
srtDeDF <- srtDeDF[! is.na(srtDeDF$SRT_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(srtDeDF$SRT_GZCP_PVAL, statistic = "pvalue", plot = TRUE)
srtDeDF$SRT_FDR <- corrected$lfdr

# Add ensembl
srtDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(srtDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### DE via LM nUMI on CPM like Seurat

exDF <- as.matrix(noCentSO@data)
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

termsDF <- data.frame(ExpCondition = metDF$Region, nUMI = metDF$nUMI)

# DE
srtuDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition+nUMI")

# Combine log2 fold changes, p-values
srtuDeDF <- data.frame(GENE = row.names(srtuDeLM$coefmat)
  , SRTnUMI_GZCP_LOG2_FC = srtuDeLM$coefmat[ ,2]
  , SRTnUMI_GZCP_PVAL = srtuDeLM$pvalmat[ ,2])
# Remove NaNs
srtuDeDF <- srtuDeDF[! is.na(srtuDeDF$SRTnUMI_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(srtuDeDF$SRTnUMI_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
srtuDeDF$SRTnUMI_FDR <- corrected$lfdr

# Add ensembl
srtuDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(srtuDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### MAGIC imputed

## Katherine (no read depth norm)

exDF <- mgExDF[ ,colnames(mgExDF) %in% metDF$Cell]
exDF <- log(exDF + 0.01, 2)

termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
mgcDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
mgcDeDF <- data.frame(GENE = row.names(mgcuDeLM$coefmat)
  , MAGIC_GZCP_LOG2_FC = mgcuDeLM$coefmat[ ,2]
  , MAGIC_GZCP_PVAL = mgcuDeLM$pvalmat[ ,2])
# Remove NaNs
mgcDeDF <- mgcDeDF[! is.na(mgcDeDF$MAGIC_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(mgcDeDF$MAGIC_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
mgcDeDF$MAGIC_FDR <- corrected$lfdr

# Add ensembl
mgcDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(mgcDeDF$GENE, bmDF$hgnc_symbol)]

## Damon (no read depth norm)

exDF <- mgdExDF[ ,colnames(mgdExDF) %in% metDF$Cell]
exDF <- log(exDF + 0.01, 2)

termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
mgdDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
mgdDeDF <- data.frame(GENE = row.names(mgdDeLM$coefmat)
  , MAGIC_Damon_GZCP_LOG2_FC = mgdDeLM$coefmat[ ,2]
  , MAGIC_Damon_GZCP_PVAL = mgdDeLM$pvalmat[ ,2])
# Remove NaNs
mgdDeDF <- mgdDeDF[! is.na(mgdDeDF$MAGIC_Damon_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(mgdDeDF$MAGIC_Damon_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
mgdDeDF$MAGIC_Damon_FDR <- corrected$lfdr

# Add ensembl
mgdDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(mgdDeDF$GENE, bmDF$hgnc_symbol)]

## Damon (read depth norm)

exDF <- mgnExDF[ ,colnames(mgnExDF) %in% metDF$Cell]
exDF <- log(exDF + 0.01, 2)

termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
mgnDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
mgnDeDF <- data.frame(GENE = row.names(mgnDeLM$coefmat)
  , MAGIC_RDN_GZCP_LOG2_FC = mgnDeLM$coefmat[ ,2]
  , MAGIC_RDN_GZCP_PVAL = mgnDeLM$pvalmat[ ,2])
# Remove NaNs
mgnDeDF <- mgnDeDF[! is.na(mgnDeDF$MAGIC_RDN_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(mgnDeDF$MAGIC_RDN_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
mgnDeDF$MAGIC_RDN_FDR <- corrected$lfdr

# Add ensembl
mgnDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(mgnDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### MAGIC imputed, regress nUMI

## No read depth normalization

exDF <- mgExDF[ ,colnames(mgExDF) %in% metDF$Cell]
exDF <- log(exDF + 0.01, 2)

termsDF <- data.frame(ExpCondition = metDF$Region, nUMI = metDF$nUMI)

# DE
mgcuDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition+nUMI")

# Combine log2 fold changes, p-values
mgcuDeDF <- data.frame(GENE = row.names(mgcuDeLM$coefmat)
  , MAGICnUMI_GZCP_LOG2_FC = mgcuDeLM$coefmat[ ,2]
  , MAGICnUMI_GZCP_PVAL = mgcuDeLM$pvalmat[ ,2])
# Remove NaNs
mgcuDeDF <- mgcuDeDF[! is.na(mgcuDeDF$MAGICnUMI_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(mgcuDeDF$MAGICnUMI_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
mgcuDeDF$MAGICnUMI_FDR <- corrected$lfdr

# Add ensembl
mgcuDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(mgcuDeDF$GENE, bmDF$hgnc_symbol)]

## Read depth normalization

exDF <- mgnExDF[ ,colnames(mgnExDF) %in% metDF$Cell]
exDF <- log(exDF + 0.01, 2)

termsDF <- data.frame(ExpCondition = metDF$Region, nUMI = metDF$nUMI)

# DE
mgnuDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition+nUMI")

# Combine log2 fold changes, p-values
mgnuDeDF <- data.frame(GENE = row.names(mgnuDeLM$coefmat)
  , MAGIC_RDNnUMI_GZCP_LOG2_FC = mgnuDeLM$coefmat[ ,2]
  , MAGIC_RDNnUMI_GZCP_PVAL = mgnuDeLM$pvalmat[ ,2])
# Remove NaNs
mgnuDeDF <- mgnuDeDF[! is.na(mgnuDeDF$MAGIC_RDNnUMI_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(mgnuDeDF$MAGIC_RDNnUMI_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
mgnuDeDF$MAGIC_RDNnUMI_FDR <- corrected$lfdr

# Add ensembl
mgnuDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(mgnuDeDF$GENE, bmDF$hgnc_symbol)]
################################################################################

### Monocle 2.0 / Census

# Uses size factor normalization:
# The size factor is the median ratio of the sample over a "pseudosample": for
# each gene, the geometric mean of all samples.

exDF <- fetb@raw.data
exDF <- exDF[ ,colnames(exDF) %in% metDF$Cell]

feature_data = data.frame(gene_short_name = rownames(exDF))
rownames(feature_data) = feature_data$gene_short_name

pd <- new("AnnotatedDataFrame", data = metDF) 
fd <- new("AnnotatedDataFrame", data = feature_data) 

mo <- newCellDataSet(as(as.matrix(exDF), "sparseMatrix"), 
  phenoData = pd, featureData = fd, lowerDetectionLimit = 0, 
  expressionFamily = negbinomial.size())

mo <- estimateSizeFactors(mo)

df <- as.matrix(exprs(mo)) / mo@phenoData$Size_Factor

# Filter to same cells as Seurat

exDF <- df[ ,colnames(df) %in% metDF$Cell]
exDF <- log(exDF + 1, 2)

termsDF <- data.frame(ExpCondition = metDF$Region)

# DE
moDeLM <- DE_Linear_Model(exDF, termsDF, "y~ExpCondition")

# Combine log2 fold changes, p-values
moDeDF <- data.frame(GENE = row.names(moDeLM$coefmat)
  , MONOCLE_GZCP_LOG2_FC = moDeLM$coefmat[ ,2]
  , MONOCLE_GZCP_PVAL = moDeLM$pvalmat[ ,2])
# Remove NaNs
moDeDF <- moDeDF[! is.na(moDeDF$MONOCLE_GZCP_PVAL), ]

# FDR
corrected <- fdrtool(moDeDF$MONOCLE_GZCP_PVAL, statistic= "pvalue", plot = TRUE)
moDeDF$MONOCLE_FDR <- corrected$lfdr

# Add ensembl
moDeDF$ENSEMBL <- bmDF$ensembl_gene_id[match(moDeDF$GENE, bmDF$hgnc_symbol)]

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))
################################################################################

### Correlation of DE

# Compile DE data frames
ldf <- list(buCpmDeDF, rawDeDF, cpmDeDF, srtDeDF, srtuDeDF, sceDeDF, sceuDeDF
  , pqDeDF, mgcDeDF, mgnDeDF, mgdDeDF, mgcuDeDF, mgnuDeDF, moDeDF)
deDF <- buDeDF
for (df in ldf) {
  deDF <- merge(deDF, df, by = "ENSEMBL")
}
df <- deDF[c("BULK_CQN_GZCP_LOG2_FC", "BULK_CPM_GZCP_LOG2_FC", "RAW_GZCP_LOG2_FC"
  , "CPM_GZCP_LOG2_FC", "SRT_GZCP_LOG2_FC", "SRTnUMI_GZCP_LOG2_FC"
  , "SCE_GZCP_LOG2_FC", "SCEnUMI_GZCP_LOG2_FC", "PQ_GZCP_LOG2_FC"
  , "MONOCLE_GZCP_LOG2_FC", "MAGIC_GZCP_LOG2_FC")]
colnames(df) <- c("BULK_CQN", "BULK_CPM", "RAW", "READ_DEPTH", 
  "SEURAT_nUMI_pMT", "SEURAT_nUMI", "SCE", "SCE_nUMI","pQ", "MONOCLE", "MAGIC")
cor(df, method = "spearman")

# GGally plot
# Spearman
png(paste0(outGraph, "CorMatrix_Heatmap_DE.png"), width = 12, height = 12
  , units = "in", res = 300)
ggcorr(df, method = c("pairwise.complete.obs", "spearman"), label = TRUE
  , label_round = 2)
dev.off()

# GGally plot
# Spearman
png(paste0(outGraph, "CorMatrix_DE_Spearman.png"), width = 18, height = 18
  , units = "in", res = 300)
ggpairs(df
  , upper = list(continuous = wrap('cor', method = "spearman", size = 10))
  , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
dev.off()
# Pearson
png(paste0(outGraph, "CorMatrix_DE_Pearson.png"), width = 18, height = 18
  , units = "in", res = 300)
ggpairs(df
  , upper = list(continuous = wrap('cor', method = "pearson", size = 10))
  , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
dev.off()

## MAGIC

df <- deDF[c("MAGIC_GZCP_LOG2_FC", "MAGIC_RDN_GZCP_LOG2_FC"
  , "MAGIC_Damon_GZCP_LOG2_FC", "MAGICnUMI_GZCP_LOG2_FC"
  , "MAGIC_RDNnUMI_GZCP_LOG2_FC")]
colnames(df) <- c("MAGIC", "MAGIC_RDN", "MAGIC_DAMON", "MAGIC_nUMI",
  "MAGIC_RDNnUMI")
cor(df, method = "spearman")

# GGally plot
# Spearman
png(paste0(outGraph, "CorMatrix_MAGIC_Heatmap_DE.png"), width = 12, height = 12
  , units = "in", res = 300)
ggcorr(df, method = c("pairwise.complete.obs", "spearman"), label = TRUE
  , label_round = 2)
dev.off()

# GGally plot
# Spearman
png(paste0(outGraph, "CorMatrix_MAGIC_DE_Spearman.png"), width = 12, height = 12
  , units = "in", res = 300)
ggpairs(df
  , upper = list(continuous = wrap('cor', method = "spearman", size = 10))
  , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
dev.off()
# Pearson
png(paste0(outGraph, "CorMatrix_MAGIC_DE_Pearson.png"), width = 18, height = 18
  , units = "in", res = 300)
ggpairs(df
  , upper = list(continuous = wrap('cor', method = "pearson", size = 10))
  , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
dev.off()



# ## what is R^2 of linear model?  how well does it fit?
# # linear model for bulk read depth normalization?
# 
# 
# ldf <- list(data.frame(sce = rowMeans(norm_exprs(sce)))
#   , data.frame(cpm = rowMeans(as.matrix(fetb@data)))
#   , data.frame(ser = rowMeans(noCentSO@scale.data))
#   , data.frame(raw = rowMeans(fetb@raw.data)))
# # Convert to ensembl
# ldf <- lapply(ldf, function(df) {
#   df$ensembl_gene_id <- bmDF$ensembl_gene_id[match(row.names(df), bmDF$hgnc_symbol)]
#   df <- df[! is.na(df$ensembl_gene_id), ]
# })
# ldf <- c(ldf, list(data.frame(ensembl_gene_id = row.names(buCpmDF)
#   , bulkcpm = rowMeans(buCpmDF))))
# mnExDF <- data.frame(bulk = rowMeans(vzcpCqnDatDF))
# mnExDF$ensembl_gene_id <- row.names(mnExDF)
# for (df in ldf) {
#   mnExDF <- merge(mnExDF, df, by = "ensembl_gene_id")
# }
# cor(mnExDF[-1], method = "spearman", use = "complete.obs")
# cor(mnExDF[-1], method = "pearson", use = "complete.obs")
# ggplot(mnExDF, aes(x = cpm, y = ser)) +
#   geom_point()
# ggsave(paste0(outGraph, "CPMvsSeurat.pdf"))
# 
# 
# pdf(paste0(outGraph, "CorMatrix_MeanExpression.pdf"))
# ggpairs(mnExDF, columns = 2:6
#   , upper = list(continuous = "cor", corMethod = "spearman")
#   , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
# dev.off()
# 
# df1 <- head(mnExDF[order(-mnExDF$bulk), ], 10000)
# cor(df1[-1], method = "spearman", use = "complete.obs")
# pdf(paste0(outGraph, "CorMatrix_MeanExpressionTop1e4.pdf"))
# ggpairs(df1, columns = 2:6
#   , upper = list(continuous = "cor", corMethod = "spearman")
#   , lower = list(continuous = wrap("points", alpha = 0.3, size = 0.3)))
# dev.off()


################################################################################

