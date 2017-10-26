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
# require(boot)

### Load data and assign variables

## Load data

# DS-002
cs1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N701/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vs1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N702/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vh1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N703/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
ch1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N704/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
# DS-003
cs2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N706/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vs2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N707/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vh2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N708/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
ch2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N709/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)

# Gene lengths and GC content for Union Exon model
# load("../source/ENSEMBLhg19_UnionAnno.rda")

# Output
# Graphs
outGraphs <- "../analysis/graphs/Covariates_DS-002-003_"
dir.create(dirname(outGraphs), recursive = TRUE)
graphsTitle <- "Covariates_DS-002-003.R"
# Data
outData <- "../analysis/Covariates_DS-002-003_"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

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

### Compile sample data frames and covariates

### Sample data frames

## Add sample name to cell ID
names(vh1ExDF)[-1] <- paste0(names(vh1ExDF)[-1], "_VZH1")
names(ch1ExDF)[-1] <- paste0(names(ch1ExDF)[-1], "_CPH1")
names(vs1ExDF)[-1] <- paste0(names(vs1ExDF)[-1], "_VZS1")
names(cs1ExDF)[-1] <- paste0(names(cs1ExDF)[-1], "_CPS1")
names(vh2ExDF)[-1] <- paste0(names(vh2ExDF)[-1], "_VZH2")
names(ch2ExDF)[-1] <- paste0(names(ch2ExDF)[-1], "_CPH2")
names(vs2ExDF)[-1] <- paste0(names(vs2ExDF)[-1], "_VZS2")
names(cs2ExDF)[-1] <- paste0(names(cs2ExDF)[-1], "_CPS2")

## Combine samples into 1 dataframe
exLDF <- list(cs1ExDF, vs1ExDF, vh1ExDF, ch1ExDF, cs2ExDF, vs2ExDF, vh2ExDF, ch2ExDF)
# exLDF <- list(chExDF, vhExDF)
for (i in 1:length(exLDF)) {
  df <- exLDF[[i]]
  if (i == 1) {
    exDF <- df}
  else {
    exDF <- merge(exDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
  }
}
str(exDF)
exDF[is.na(exDF)] <- 0
row.names(exDF) <- exDF$GENE
exDF <- exDF[ ,-1]
print("Number of cells input:")
print(ncol(exDF))

### Covariates data frame

# Dissociation method
cvDF <- data.frame(DISSOCIATION = rep("hard", ncol(exDF)))
cvDF$DISSOCIATION <- as.character(cvDF$DISSOCIATION)
idx <- grep("S", names(exDF))
cvDF$DISSOCIATION[idx] <- "soft"

# Cell batch
cvDF$BATCH <- rep("1", ncol(exDF))
cvDF$BATCH <- as.character(cvDF$BATCH)
idx <- grep("2", names(exDF))
cvDF$BATCH[idx] <- "2"

# Brain region
cvDF$REGION <- rep("VZ", ncol(exDF))
cvDF$REGION <- as.character(cvDF$REGION)
idx <- grep("CP", names(exDF))
cvDF$REGION[idx] <- "CP"
################################################################################

### PCA Correlation Matrix

## PCA of expression

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(exDF), scale = FALSE))
# Run PCA
pCdat <- prcomp(meanCenteredM, center = FALSE);
topPCs <- pCdat$rotation[,1:10];
# Calculate variance explained by each PC
varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

## Correlation matrix of expression PCs and technical covariates

pairsDat <- data.frame(cvDF)
# Convert character columns to factor
pairsDat[sapply(pairsDat, is.character)] <- lapply(pairsDat[sapply(pairsDat
  , is.character)], as.factor)

cond <- labels2colors(cvDF$REGION)  ## colors

pdf(paste0(outGraphs, "CovariatesMatrix_FtMm250.pdf"), height = 12, width = 12)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()
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


