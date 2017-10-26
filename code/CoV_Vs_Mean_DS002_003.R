# Damon Polioudakis
# 2017-02-21
# Plot mean expression vs CoV graph colored by detection rate

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

require(ggplot2)
require(Seurat)
require(VennDiagram)

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

# Seurat clustering (for variable genes)
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
serO <- fetb

# Cell cycle markers from Macosko 2015 Table S2
ccDF <- read.csv("../source/Macosko_2015_ST2_CellCycle.csv", header = TRUE
  , fill = TRUE)

## Variables
outGraphPfx <- "../analysis/graphs/CoV_Vs_Mean_DS002_003_"
graphCodeTitle <- "CoV_Vs_Mean_DS002_003.R"

## ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_set(theme_get() + theme(axis.line = element_line(colour = "black")
  , panel.border = element_blank()))
theme_update(plot.title = element_text(size = 12))
################################################################################

### Functions

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

# Calculate coefficient of variation
Coefficent_Of_Variation <- function(mean, sd){
  (sd/mean)*100
}
################################################################################

### Format data

## Combine samples into 1 data frame

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
################################################################################

### Detection rate, CoV, and mean expression

# Calculate mean expression for each gene
mnL <- rowMeans(dgeDF)

# Calculate CoV for each gene
# Standard deviation for each gene
sdL <- apply(dgeDF, 1, sd)
# CoV for each gene
Coefficent_Of_Variation(mnL[[1]], sdL[[1]])
cvL <- mapply(Coefficent_Of_Variation, mnL, sdL)

# Detection rate
dtL <- apply(dgeDF, 1, function(geneDge) {
  sum(geneDge >= 1) / ncol(dgeDF) * 100}
)

# Data frame of mean and CoV for ggplot
ggDF <- data.frame(MEAN = mnL, COV = cvL, DETECTION_RATE = dtL)

# ggplot mean, CoV, and detection
pdf(paste0(outGraphPfx, "Vs_DetectionRate.pdf"), width = 7, height = 5)
ggplot(data = ggDF, aes(x = MEAN, y = log(COV,2), col = DETECTION_RATE)) +
  geom_point(size = 0.75) +
  scale_color_gradient(name = "Detection rate (percent)", high = "red"
    , low = "blue") +
  scale_shape(solid = FALSE) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by detection rate"))
ggplot(data = ggDF, aes(x = MEAN, y = log(COV, 2), col = DETECTION_RATE)) +
  geom_point(size = 0.75) +
  scale_color_gradient(name = "Detection rate (percent)", high = "red"
    , low = "blue") +
  scale_shape(solid = FALSE) +
  coord_cartesian(xlim = c(0, 10)) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by detection rate"
    , "\nAxis limits adjusted"))
dev.off()

# Seurat variable genes
ggDF$SEURAT_VARIABLE_GENE <- "Not variable gene"
ggDF$SEURAT_VARIABLE_GENE[row.names(ggDF) %in% serO@var.genes] <- "Variable gene"

# ggplot mean, CoV, and Seurat variable genes
pdf(paste0(outGraphPfx, "Vs_VariableGenes.pdf"), width = 7, height = 5)
ggplot(data = ggDF, aes(x = MEAN, y = log(COV, 2), col = SEURAT_VARIABLE_GENE)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_discrete(name = "Variable genes") +
  scale_shape(solid = FALSE) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by Seurat variable genes"))
ggplot(data = ggDF, aes(x = MEAN, y = log(COV, 2), col = SEURAT_VARIABLE_GENE)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_discrete(name = "Variable genes") +
  scale_shape(solid = FALSE) +
  coord_cartesian(xlim = c(0, 10)) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by Seurat variable genes"
    , "\nAxis limits adjusted"))
dev.off()

# Macosko cell cycle genes
# Factors to character
ccDF <- data.frame(lapply(ccDF, as.character), stringsAsFactors = FALSE)
# Merge columns into one list
ccL <- do.call(c, ccDF)
# Clean up extra spaces in gene symbols
ccL <- gsub(" *", "", ccL)
# Add to ggplot data frame
ggDF$CELLCYCLE_GENE <- "Not cell cycle gene"
ggDF$CELLCYCLE_GENE[row.names(ggDF) %in% ccL] <- "Cell cycle gene"

# ggplot mean, CoV, and Macosko cell cycle genes
pdf(paste0(outGraphPfx, "Vs_CellCycleGenes.pdf"), width = 7, height = 5)
ggplot(data = ggDF, aes(x = MEAN, y = log(COV, 2), col = CELLCYCLE_GENE)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_discrete(name = "Cell cycle genes") +
  # scale_shape(solid = FALSE) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by Macosko cell cycle genes"))
ggplot(data = ggDF, aes(x = MEAN, y = log(COV, 2), col = CELLCYCLE_GENE)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_discrete(name = "Cell cycle genes") +
  # scale_shape(solid = FALSE) +
  coord_cartesian(xlim = c(0, 10)) +
  xlab("log2(mean DGE)") +
  ylab("log2 coefficient of variance") +
  ggtitle(paste0(graphCodeTitle
    , "\nLog2 mean DGE versus log2 coefficient of variance"
    , "\nColored by Macosko cell cycle genes"
    , "\nAxis limits adjusted"))
dev.off()


# Intersection variable genes and cell cycle genes: 174
length(intersect(row.names(ggDF)[ggDF$SEURAT_VARIABLE_GENE == "Variable gene"]
  , row.names(ggDF)[ggDF$CELLCYCLE_GENE == "Cell cycle gene"]))
# Variable genes: 4277
length(row.names(ggDF)[ggDF$SEURAT_VARIABLE_GENE == "Variable gene"])
# Cell cycle genes: 555
length(row.names(ggDF)[ggDF$CELLCYCLE_GENE == "Cell cycle gene"])
# Venn diagram
pdf(paste0(outGraphPfx, "Intersection_CellCycle_Variable_Genes.pdf"))
draw.pairwise.venn(area1 = length(row.names(ggDF)[ggDF$SEURAT_VARIABLE_GENE == "Variable gene"])
  , area2 = length(row.names(ggDF)[ggDF$CELLCYCLE_GENE == "Cell cycle gene"])
  , cross.area = length(intersect(row.names(ggDF)[ggDF$SEURAT_VARIABLE_GENE == "Variable gene"]
    , row.names(ggDF)[ggDF$CELLCYCLE_GENE == "Cell cycle gene"]))
  , category = c("Variable Genes", "Cell Cycle Genes")
  , cat.pos = c(0, 0)
  , fill = c("#a6cee3", "#fb9a99"))
dev.off()
################################################################################


