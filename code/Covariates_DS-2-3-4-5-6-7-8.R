# Damon Polioudakis
# 2017-02-13
# Normalize Luis' bulk RNAseq from experiments for ATAC work

### TODO
# Add bootstrapping regression for regressing out covariates
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(reshape2)
require(WGCNA)
require(Seurat)
require(irlba)
require(biomaRt)
# require(boot)

### Load data and assign variables

## Load data

# Digital gene expression and metadata: exDF and metDF
load("../analysis/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7-8.Rdata")
# # Subsetting for testing
# idx <- sample(1:ncol(exDF), 500)
# exDF <- exDF[idx]
# metDF <- metDF[idx, ]
# Seurat object for variable genes
load("../analysis/Seurat_Cluster_DS-2-3-4-5-6-7-8/Seurat_Cluster_DS-2-3-4-5-6-7-8_exon_FtMm250_RegNumiPMtLibBrain_seuratO.Robj")

# CDS length and GC
geneInfoDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh37_Ensembl75.csv"
  , header = TRUE)

# Gene lengths and GC content for Union Exon model
# load("../source/ENSEMBLhg19_UnionAnno.rda")

# Output
# Graphs
outGraphs <- "../analysis/graphs/Covariates_DS-2-3-4-5-6-7-8_"
dir.create(dirname(outGraphs), recursive = TRUE)
graphsTitle <- "Covariates_DS-2-3-4-5-6-7-8.R"
# Data
outData <- "../analysis/Covariates_DS-2-3-4-5-6-7-8_"

## ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 16)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , panel.border = element_blank()
)
################################################################################

### Functions

## Function: Query biomaRt
QueryBiomaRt <- function (genesList, filters, attributes) {
  genesDF <- data.frame(genesList)
  # Attribute: values to retrieve - c("hgnc_symbol", "ensembl_gene_id")
  # Filters: input query type
  # Values: input query
  ensembl= useMart(biomart="ENSEMBL_MART_ENSEMBL"
    , host="feb2014.archive.ensembl.org", path="/biomart/martservice"
    , dataset="hsapiens_gene_ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Outputs data frame of attributes parameter from biomart
  bmDF <- getBM(
    attributes = attributes
    , filters = filters
    , values = genesDF
    , mart = ensembl
  )
  bmDF
}

# Useful function for comparing multivariate data
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
################################################################################

### Format and compile covariates

print("Formatting and compiling covariates")

## Add UMI and genes detected to metadata
# Number UMI
nUMI <- colSums(exDF)
metDF$nUMI <- nUMI[match(metDF$CELL, names(nUMI))]
# Number genes detected (counts > 0)
nGene <- apply(exDF, 2, function(col) {sum(col  > 0)})
metDF$nGene <- nGene[match(metDF$CELL, names(nGene))]

## Length and GC
# Convert ensembl ID to gene symbol and subset to ensembl IDs found in biomaRt
df <- QueryBiomaRt(row.names(exDF), "hgnc_symbol"
  , c("ensembl_gene_id", "hgnc_symbol"))
geneInfoDF <- df
# Add gene information
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "gene_biotype", "percentage_gc_content"))
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "cds_length"))
# Take longest CDS
df <- aggregate(cds_length ~ ensembl_gene_id, data = df, max)
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
str(geneInfoDF)
head(geneInfoDF)
tail(geneInfoDF)

# Length and GC
# Length
# Take longest CDS
df <- aggregate(cds_length ~ hgnc_symbol + percentage_gc_content
  , data = geneInfoDF, max)
df <- merge(geneInfoDF, exDF, by.x = "hgnc_symbol", by.y = "row.names")
df <- na.omit(df)
# Length calculation per cell
lengthDF <- data.frame(apply(df[6:ncol(df)], 2, function(col) {
  sum((col * df$cds_length)) / sum(col)
}))
colnames(lengthDF) <- "CDS_LENGTH"
metDF$CDS_LENGTH <- lengthDF[match(metDF$CELL, row.names(lengthDF)), ]
# GC
# Take GC for longest CDS
df <- aggregate(cds_length ~ hgnc_symbol + percentage_gc_content
  , data = geneInfoDF, max)
df <- merge(df, exDF, by.x = "hgnc_symbol", by.y = "row.names")
df <- na.omit(df)
# GC calculation per cell
gcDF <- data.frame(apply(df[6:ncol(df)], 2, function(col) {
  sum((col * df$percentage_gc_content)) / sum(col)
}))
colnames(gcDF) <- "GC_CONTENT"
metDF$GC_CONTENT <- gcDF[match(metDF$CELL, row.names(gcDF)), ]

# Percentage of mitochondrial genes
mito.genes <- grep("^MT-", rownames(exDF), value = TRUE)
percent.mito <- apply(exDF[row.names(exDF) %in% mito.genes, ], 2, sum) / apply(exDF, 2, sum)
metDF$PERCENT_MT <- percent.mito

## Covariates
cvDF <- metDF[c("BRAIN", "REGION", "LIBRARY", "SEQ_RUN", "nUMI", "nGene"
  , "CDS_LENGTH", "GC_CONTENT", "PERCENT_MT")]
################################################################################

### PCA Correlation Matrix

## PCA of expression - raw counts - all genes

print("PCA of expression - raw counts - all genes")

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(exDF), scale = FALSE))
# Run PCA
pcAll <- prcomp_irlba(meanCenteredM, n = 200, center = FALSE);
topPCs <- pcAll$rotation[ ,1:5];
# Calculate variance explained by each PC
varExp <- (pcAll$sdev)^2 / sum(pcAll$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

## Correlation matrix of expression PCs and technical covariates

pairsDat <- data.frame(cvDF)
# Convert character columns to factor
pairsDat[sapply(pairsDat, is.character)] <- lapply(pairsDat[sapply(pairsDat
  , is.character)], as.factor)

cond <- labels2colors(cvDF$LIBRARY)  ## colors

png(paste0(outGraphs, "CovariatesMatrix_FtMm250_CountsAllGenes.png")
  , height = 12, width = 12, units = "in", res = 300)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = paste0(
    "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values"
    ,"\nRaw counts - All genes"))
dev.off()

## PCA of expression - raw counts - Seurat variable genes

print("PCA of expression - raw counts - Seurat variable genes")

# Subset expression dataframe to Seurat variable genes
varDF <- exDF[row.names(exDF) %in% centSO@var.genes, ]

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(varDF), scale = FALSE))
# Run PCA
pcVar <- prcomp_irlba(meanCenteredM, n = 200, center = FALSE);
topPCs <- pcVar$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pcVar$sdev)^2 / sum(pcVar$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

## Correlation matrix of expression PCs and technical covariates

pairsDat <- data.frame(cvDF)
# Convert character columns to factor
pairsDat[sapply(pairsDat, is.character)] <- lapply(pairsDat[sapply(pairsDat
  , is.character)], as.factor)

cond <- labels2colors(cvDF$LIBRARY)  ## colors

png(paste0(outGraphs, "CovariatesMatrix_FtMm250_CountsSeuratVarGenes.png")
  , height = 12, width = 12, units = "in", res = 300)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = paste0(
    "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values"
    ,"\nRaw counts - variable genes"))
dev.off()

## PCA of expression - Seurat normalized - Seurat variable genes

print("PCA of expression - Seurat normalized - Seurat variable genes")

# Subset expression dataframe to Seurat variable genes
nmVarDF <- centSO@scale.data[row.names(centSO@scale.data) %in% centSO@var.genes, ]
# Subset covariates to cells remaining after Seurat filters
nmCvDF <- cvDF[match(colnames(nmVarDF), row.names(cvDF)), ]

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(nmVarDF), scale = FALSE))
# Run PCA
pcNmVar <- prcomp_irlba(meanCenteredM, n = 200, center = FALSE);
topPCs <- pcNmVar$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pcNmVar$sdev)^2 / sum(pcNmVar$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

## Correlation matrix of expression PCs and technical covariates

pairsDat <- data.frame(nmCvDF)
# Convert character columns to factor
pairsDat[sapply(pairsDat, is.character)] <- lapply(pairsDat[sapply(pairsDat
  , is.character)], as.factor)

cond <- labels2colors(nmCvDF$LIBRARY)  ## colors

png(paste0(outGraphs, "CovariatesMatrix_FtMm250_SeuratNormVarGenes.png")
  , height = 12, width = 12, units = "in", res = 300)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = paste0(
    "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values"
    ,"\nSeurat normalized, centered, scaled - Variable genes"))
dev.off()

## Save PCA
save(pcAll, pcVar, pcNmVar, file = paste0(outData, "PCA.Rdata"))
################################################################################

## Bootstrap regression
# 
# covDF <- data.frame(metDatDF[ ,c("ExpCondition", "RIN.y", "X260.230"
#   , "X260.280", "ExtractionDate.y")], topPCdatSeq)
# covDF$ExpCondition <- as.numeric(covDF$ExpCondition)
# covDF$ExtractionDate.y <- as.numeric(covDF$ExtractionDate.y)
# 
# # Factors must be coded as numeric for bootstrap
# RegressCovariatesBootstrap <- function (exM, covDF) {
#   # 1) lme using covariates on unnormalized data
#   options(stringsAsFactors=FALSE)
#   boot <- FALSE
#   numboot <- 10
#   bs <- function(formula, data, indices) {
#     d <- data[indices,] # allows bootstrap function to select samples
#     fit <- lm(formula, data=d)
#     return(coef(fit))
#   }
#   exBootRegM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
#   rownames(exBootRegM) <- rownames(exM)
#   colnames(exBootRegM) <- colnames(exM)
#   # ncol(covDF)+1 when condition has 2 levels
#   coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
#   set.seed(8675309)
#   for (i in 1:nrow(exM)) {
#     if (i%%1000 == 0) {print(i)}
#     thisexp <- as.numeric(exM[i, ])
#     #       bs.results <- boot(data = data.frame(thisexp, covDF), statistic = bs,
#     #         R = numboot, formula = thisexp~. + age + sex + PMI)
#     bs.results <- boot(data = data.frame(thisexp, covDF), statistic = bs,
#       R = numboot, formula = thisexp~ExpCondition+RIN.y+X260.230+X260.280+ExtractionDate.y+Seq.PC1+Seq.PC2+Seq.PC3+Seq.PC4+Seq.PC5)
#     ## get the median - we can sometimes get NA values here... so let's exclude
#     ## these - old code #bs.stats <- apply(bs.results$t,2,median) 
#     bs.stats <- rep(NA, ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
#     for (n in 1:ncol(bs.results$t)) {
#       bs.stats[n] <- median(na.omit(bs.results$t[ ,n]))
#     }
#     coefmat[i,] <- bs.stats
#     exBootRegM[i,] <- thisexp - bs.stats[3]*covDF[ ,2] - bs.stats[4]*covDF[ ,3] - bs.stats[5]*covDF[ ,4] - bs.stats[6]*covDF[ ,5] - bs.stats[7]*covDF[ ,6] - bs.stats[8]*covDF[ ,7] - bs.stats[9]*covDF[ ,8] - bs.stats[10]*covDF[ ,9] - bs.stats[11]*covDF[ ,10]
#     # cat('Done for Gene',i,'\n')
#   }
# }
# 


