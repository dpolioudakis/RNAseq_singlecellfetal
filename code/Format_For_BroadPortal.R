# Damon Polioudakis
# 2017-02-21
# Test Seurat PCA method and clustering parameters

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(seurat)

## Input data

# Digital gene expression
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

# Clustering
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
seuratO <- fetb

## Output
################################################################################

### Compile gene expression data frames

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
# Change NAs to 0s (from merging samples and may not detect a gene in 1 sample)
exDF[is.na(exDF)] <- 0
print("Number of cells input:")
print(ncol(exDF))
# Write to table
write.table(exDF, file = "../analysis/Merged_DGE_DS002003.txt", quote = FALSE
  , sep = "\t", row.names = FALSE)
system(paste0("gzip ../analysis/Merged_DGE_DS002003.txt"))

### Metadata

## Add sample name to cell ID
metDF <- data.frame(NAME = c("TYPE", names(exDF[-1])))
metDF$REGION <- c("group"
  , rep("GZ", length(names(vh1ExDF)[-1]))
  , rep("CP", length(names(ch1ExDF)[-1]))
  , rep("GZ", length(names(vs1ExDF)[-1]))
  , rep("CP", length(names(cs1ExDF)[-1]))
  , rep("GZ", length(names(vh2ExDF)[-1]))
  , rep("CP", length(names(ch2ExDF)[-1]))
  , rep("GZ", length(names(vs2ExDF)[-1]))
  , rep("CP", length(names(cs2ExDF)[-1]))
)
metDF$LIBRARY <- c("group"
  , rep("1", length(names(vh1ExDF)[-1]))
  , rep("1", length(names(ch1ExDF)[-1]))
  , rep("1", length(names(vs1ExDF)[-1]))
  , rep("1", length(names(cs1ExDF)[-1]))
  , rep("2", length(names(vh2ExDF)[-1]))
  , rep("2", length(names(ch2ExDF)[-1]))
  , rep("2", length(names(vs2ExDF)[-1]))
  , rep("2", length(names(cs2ExDF)[-1]))
)
metDF$CELL_BATCH <- c("group"
  , rep("1", length(names(vh1ExDF)[-1]))
  , rep("2", length(names(ch1ExDF)[-1]))
  , rep("3", length(names(vs1ExDF)[-1]))
  , rep("4", length(names(cs1ExDF)[-1]))
  , rep("5", length(names(vh2ExDF)[-1]))
  , rep("6", length(names(ch2ExDF)[-1]))
  , rep("7", length(names(vs2ExDF)[-1]))
  , rep("8", length(names(cs2ExDF)[-1]))
)

dim(metDF)
head(metDF)
tail(metDF)
# Write to table
write.table(metDF, file = "../analysis/MetaData_DS002003.txt", quote = FALSE
  , sep = "\t", row.names = FALSE)

### Clusters

str(seuratO@tsne.rot)
clustersDF <- data.frame(
  NAME = c("TYPE", row.names(seuratO@tsne.rot))
  , X = c("numeric", seuratO@tsne.rot$tSNE_1)
  , Y = c("numeric", seuratO@tsne.rot$tSNE_2)
)
dim(clustersDF)
head(clustersDF)
tail(clustersDF)
write.table(clustersDF, file = "../analysis/Clusters_DS002003.txt", quote = FALSE
  , sep = "\t", row.names = FALSE)
################################################################################