# Damon Polioudakis
# 2017-02-17
# Plot Seurat variable genes on UMI collapse / non collapsed vs mean expression
# graph to see importance of UMIs
################################################################################

require(ggplot2)

## Inputs

# DS-002
pntDir <- "../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/"
dirs <- c("N701", "N702", "N703", "N704")
# DGE
inD2DgeL <- paste0(pntDir, dirs, "/out_gene_exon_tagged_dge_FtMm250.txt")
d2DgeLDF <- lapply(inD2DgeL, function(inD2dge) read.table(inD2dge, header = TRUE))
# Counts (pre UMI collapse)
inD2CtL <- paste0(pntDir, dirs, "/out_gene_exon_tagged_counts_FtMm250.txt")
d2CtLDF <- lapply(inD2CtL, function(inD2Ct) read.table(inD2Ct, header = TRUE))

# DS-003
pntDir <- "../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/"
dirs <- c("N706", "N707", "N708", "N709")
# DGE
inD3DgeL <- paste0(pntDir, dirs, "/out_gene_exon_tagged_dge_FtMm250.txt")
d3DgeLDF <- lapply(inD3DgeL, function(inD3dge) read.table(inD3dge, header = TRUE))
# Counts (pre UMI collapse)
inD3CtL <- paste0(pntDir, dirs, "/out_gene_exon_tagged_counts_FtMm250.txt")
d3CtLDF <- lapply(inD3CtL, function(inD3Ct) read.table(inD3Ct, header = TRUE))

## Variables
outGraph <- "../analysis/graphs/UMI_Pre_Vs_Post_Collapse_DS002_003"
graphCodeTitle <- "UMI_Pre_Vs_Post_Collapse_DS002_003.R"

## ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 12))
################################################################################

## Combine samples into 1 data frame

# Function to merge list of data frames
Merge_Data_Frames <- function (dfL, mergeBy) {
  for (i in 1:length(dfL)) {
    df <- dfL[[i]]
    if (i == 1) {
      mergeDF <- df}
    else {
      mergeDF <- merge(mergeDF, df, by.x = mergeBy, by.y = mergeBy, all = TRUE)  
    }
  }
  return(mergeDF)
}

# Merge DGE data frames
dgeLDF <- append(d2DgeLDF, d3DgeLDF)
dgeDF <- Merge_Data_Frames(dgeLDF, "GENE")
str(dgeDF)
dgeDF[is.na(dgeDF)] <- 0
row.names(dgeDF) <- dgeDF$GENE
dgeDF <- dgeDF[ ,-1]
print("Number of cells input:")
print(ncol(dgeDF))

# Merge count data frames
ctLDF <- append(d2CtLDF, d3CtLDF)
ctDF <- Merge_Data_Frames(ctLDF, "GENE")
str(ctDF)
ctDF[is.na(ctDF)] <- 0
row.names(ctDF) <- ctDF$GENE
ctDF <- ctDF[ ,-1]
print("Number of cells input:")
print(ncol(ctDF))

# Calculate mean expression for each gene
mnDF <- data.frame(DGE = rowMeans(dgeDF), COUNTS = rowMeans(ctDF))

# Order by mean DGE
mnDF <- mnDF[order(mnDF$DGE), ]

# Log2 transform
mnDF <- log(mnDF + 0.01, 2)

# ggplot
ggplot(data = mnDF, aes(x = DGE, y = COUNTS)) +
  geom_point(alpha = 0.25) +
  scale_shape(solid = FALSE) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("log2(mean DGE)") +
  ylab("log2(mean counts)") +
  ggtitle(paste0(graphCodeTitle,
    "\nMean gene expression before and after UMI collapse"))
ggsave(paste0(outGraph, ".pdf"), width = 5, height = 5)

## Variance versus expression level

# Counts
# Variance
ctVrDF <- apply(ctDF, 1, var)
# Mean expression
ctMn <- rowMeans(ctDF)
# Variance / (mean expression)^2 and log2(mean expression)
ggDF <- data.frame(STANDARDIZED_VARIANCE = (ctVrDF/ ctMn^2), MEAN = log(ctMn, 2))

# ggplot
pdf(paste0(outGraph, "_Variance_Counts.pdf"), width = 5, height = 5)
ggplot(data = ggDF, aes(x = MEAN, y = STANDARDIZED_VARIANCE)) +
  geom_point(alpha = 0.25) +
  scale_shape(solid = FALSE) +
  xlab("log2(mean DGE)") +
  ylab("Standardized variance") +
  ggtitle(paste0(graphCodeTitle,
    "\nStandarized variance of counts versus mean expression counts"))
ggplot(data = ggDF, aes(x = MEAN, y = STANDARDIZED_VARIANCE)) +
  geom_point(alpha = 0.25) +
  scale_shape(solid = FALSE) +
  coord_cartesian(ylim = c(0, 50)) +
  xlab("log2(mean DGE)") +
  ylab("Standardized variance") +
  ggtitle(paste0(graphCodeTitle,
    "\nStandarized variance of counts versus mean expression counts"))
dev.off()

# DGE
# Variance
dgeVrDF <- apply(dgeDF, 1, var)
# Mean expression
dgeMn <- rowMeans(dgeDF)
# Variance / (mean expression)^2 and log2(mean expression)
ggDF <- data.frame(STANDARDIZED_VARIANCE = (dgeVrDF/ dgeMn^2), MEAN = log(dgeMn, 2))

# ggplot
pdf(paste0(outGraph, "_Variance_DGE.pdf"), width = 5, height = 5)
ggplot(data = ggDF, aes(x = MEAN, y = STANDARDIZED_VARIANCE)) +
  geom_point(alpha = 0.25) +
  scale_shape(solid=FALSE) +
  xlab("log2(mean DGE)") +
  ylab("Standardized variance") +
  ggtitle(paste0(graphCodeTitle,
    "\nStandarized variance of DGE versus mean expression DGE"))
ggplot(data = ggDF, aes(x = MEAN, y = STANDARDIZED_VARIANCE)) +
  geom_point(alpha = 0.25) +
  scale_shape(solid = FALSE) +
  coord_cartesian(ylim = c(0, 50)) +
  xlab("log2(mean DGE)") +
  ylab("Standardized variance") +
  ggtitle(paste0(graphCodeTitle,
    "\nStandarized variance of DGE versus mean expression DGE"))
dev.off()
################################################################################