
# Damon Polioudakis
# 2016-12-08
# Plot bulk RNAseq of VZ and CP from Luis and Jason's ATAC versus pooled
# scRNAseq VZ and CP
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(biomaRt)

### Load data and assign variables

## Load data

# DS-002 digital expression
cs1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N701/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vs1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N702/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vh1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N703/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
ch1ExDF <- read.table("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N704/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)

# DS-003
cs2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N706/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vs2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N707/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
vh2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N708/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)
ch2ExDF <- read.table("../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8/N709/out_gene_exon_tagged_dge_FtMm250.txt", header = TRUE)

# Bulk RNAseq
buExDF <- read.csv("../bulk_VZ_CP_from_ATAC/data/htseq/Exprs_HTSCexon.csv"
  , row.names = 1)
# Bulk metadata
buMetDF <- read.csv("../bulk_VZ_CP_from_ATAC/metadata/VZCP_sampleinfo.csv", header = TRUE)
# Picard Sequencing Statistics - bulk RNAseq
picStatsBuDF <- read.table("../../bulk_VZ_CP_from_ATAC/metadata/alignment_summary.txt", fill = TRUE
  , header = TRUE)


## Variables

graphCodeTitle <- "Compare_To_Bulk_DS002_003.R"
outGraph <- "../analysis/graphs/Compare_To_Bulk_DS002_003_"
outTable <- "../analysis/tables/Compare_To_Bulk_DS002_003_"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 12))
################################################################################

### Functions

# Function to convert list of Ensembl IDs to data frame of Ensembl IDs and Gene Symbols
ConvertEnsemblToGeneSym <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneDF <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol")
    , filters = "ensembl_gene_id"
    , values = geneDF
    , mart = ensembl
  )
  ensemblGeneSymDF
}
################################################################################

### Format

## Subset bulk to ensembl with hgnc_symbols

## Remove chromosome number from ensembl IDs
row.names(buExDF) <- gsub("\\..*", "", row.names(buExDF))

## Convert ensembl IDs to gene symbols
enHgDF <- ConvertEnsemblToGeneSym(row.names(buExDF))
# Remove genes with no hgnc symbol
enHgDF <- enHgDF[! enHgDF$hgnc_symbol == "", ]
# Remove duplicate hgnc symbols
enHgDF <- enHgDF[! duplicated(enHgDF$hgnc_symbol), ]
# Merge gene expression data frame and ensembl gene symbol data frame
buExDF <- merge(enHgDF, buExDF, by.x = "ensembl_gene_id", by.y = "row.names")


## DS-002 DS-003 combine VZ samples and CP samples into data frames

## Add sample name to cell ID
names(vh1ExDF)[-1] <- paste0(names(vh1ExDF)[-1], "_VZH")
names(ch1ExDF)[-1] <- paste0(names(ch1ExDF)[-1], "_CPH")
names(vs1ExDF)[-1] <- paste0(names(vs1ExDF)[-1], "_VZS")
names(cs1ExDF)[-1] <- paste0(names(cs1ExDF)[-1], "_CPS")
names(vh2ExDF)[-1] <- paste0(names(vh2ExDF)[-1], "_VZH")
names(ch2ExDF)[-1] <- paste0(names(ch2ExDF)[-1], "_CPH")
names(vs2ExDF)[-1] <- paste0(names(vs2ExDF)[-1], "_VZS")
names(cs2ExDF)[-1] <- paste0(names(cs2ExDF)[-1], "_CPS")

## Combine samples into 1 dataframe
# VZ
exLDF <- list(vs1ExDF, vh1ExDF, vs2ExDF, vh2ExDF)
# exLDF <- list(chExDF, vhExDF)
for (i in 1:length(exLDF)) {
  df <- exLDF[[i]]
  if (i == 1) {
    vzExDF <- df}
  else {
    vzExDF <- merge(vzExDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
  }
}
# CP
exLDF <- list(cs1ExDF, ch1ExDF, cs2ExDF, ch2ExDF)
# exLDF <- list(chExDF, vhExDF)
for (i in 1:length(exLDF)) {
  df <- exLDF[[i]]
  if (i == 1) {
    cpExDF <- df}
  else {
    cpExDF <- merge(cpExDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
  }
}

vzExDF[is.na(vzExDF)] <- 0
row.names(vzExDF) <- vzExDF$GENE
vzExDF <- vzExDF[ ,-1]
print("Number of cells input VZ:")
print(ncol(vzExDF))

cpExDF[is.na(cpExDF)] <- 0
row.names(cpExDF) <- cpExDF$GENE
cpExDF <- cpExDF[ ,-1]
print("Number of cells input CP:")
print(ncol(cpExDF))
################################################################################

### Process data - Mean counts, Filter, read depth normalize, and pool

## Change column names to brain region
colnames(buExDF)[-c(1,2)] <- as.character(buMetDF$ExpCondition)

## Split bulk into VZ and CP
vzBuExDF <- buExDF[ ,colnames(buExDF) == "VZ"]
row.names(vzBuExDF) <- buExDF$hgnc_symbol
cpBuExDF <- buExDF[ ,colnames(buExDF) == "CP"]
row.names(cpBuExDF) <- buExDF$hgnc_symbol

# ## Read depth normalize by number of mapped reads
# 
# # Read depth normalize bulk by number mapped reads
# buCpmDF <- buExDF / (picStatsBuDF$PF_READS_ALIGNED / 10^6)

# Pool scRNAseq
pVzExDF <- data.frame(GENE = row.names(vzExDF)
  , POOLED = apply(vzExDF, 1, sum))
head(pVzExDF)
tail(pVzExDF)
pCpExDF <- data.frame(GENE = row.names(cpExDF)
  , POOLED = apply(cpExDF, 1, sum))
head(pCpExDF)
tail(pCpExDF)


## Mean counts for bulk RNAseq

# # Mean counts
# mnBuEx <- apply(buCpmDF, 1, mean)
# head(mnBuEx)

#VZ
mnVzBuExDF <- data.frame(GENE = row.names(vzBuExDF)
  , MEAN_EXPRESSION = apply(vzBuExDF, 1, mean))
head(mnVzBuExDF)

#CP
mnCpBuExDF <- data.frame(GENE = row.names(cpBuExDF)
  , MEAN_EXPRESSION = apply(cpBuExDF, 1, mean))
head(mnCpBuExDF)

Plot_Bulk_SC_Correlation <- function (pooled, meanBulk, region, color) {
  # Merge pooled SC and bulk mean
  pbDF <- merge(pooled, meanBulk, by.x = "GENE", by.y = "GENE")
  # Log2(counts + 1)
  pbDF[ ,2:3] <- log(pbDF[ ,2:3] + 1, 2)
  # Order by pooled
  pbDF <- pbDF[order(pbDF$POOLED), ]
  pbDF$GENE <- factor(pbDF$GENE, levels = as.character(pbDF$GENE))
  # Spearman correlation of log2(counts + 1) pooled vs bulk mean
  sprCorMn <- round(cor(pbDF$POOLED, pbDF$MEAN_EXPRESSION
    , method = "spearman"), 2)
  # ggplot2
  p <- ggplot(pbDF, aes(x = POOLED, y = MEAN_EXPRESSION)) +
    geom_point(alpha = 0.4, shape = 1, size = 0.7, color = color) +
    theme_bw(base_size = 18) +
    ylab("Bulk RNA-seq: log2(Mean Counts + 1)") +
    xlab("scRNA-seq Pooled: log2(Counts + 1)") +
    ggtitle(paste0(graphCodeTitle
      , paste0("\nPooled scRNAseq vs Bulk RNAseq - ", region)
      , "\nSpearman correlation: ", sprCorMn))
  print(p)
}

pdf(paste0(outGraph, "Pooled_Vs_Bulk_log2p1.pdf"))
Plot_Bulk_SC_Correlation(pVzExDF, mnVzBuExDF, "VZ", "blue")
Plot_Bulk_SC_Correlation(pCpExDF, mnCpBuExDF, "CP", "red")
dev.off()


### Split into deciles and check spearman correlation for each decile

# Merge pooled SC and bulk mean
pbDF <- merge(pVzExDF, mnVzBuExDF, by.x = "GENE", by.y = "GENE")
# Log2(counts + 1)
pbDF[ ,2:3] <- log(pbDF[ ,2:3] + 1, 2)
pbDF <- pbDF[order(-pbDF$MEAN_EXPRESSION), ]

## Expression deciles
splits <- seq(0, range(pbDF$POOLED)[2], length.out = 11)

pbDF$split <- with(pbDF, cut(POOLED, 
  breaks = splits, 
  include.lowest = TRUE))
splitDF <- split(pbDF, pbDF$split)

# Spearman correlation for each range
df <- data.frame(sapply(splitDF, function(subgroup) cor(subgroup$POOLED
  , subgroup$MEAN_EXPRESSION, method = "spearman")))
colnames(df) <- "SPEARMAN"
df <- round(df, 2)
write.table(df, paste0(outTable, "Expression_Deciles.txt"), sep = "\t", quote = FALSE)


## Gene deciles ordered by bulk mean expression (i.e. Genes 1-1000, 2-2000, etc)

splitDF <- split(pbDF, rep(1:10, each = 1758))
df <- data.frame(sapply(splitDF, function(subgroup) cor(subgroup$POOLED
  , subgroup$MEAN_EXPRESSION, method = "spearman")))
colnames(df) <- "SPEARMAN"
df <- round(df, 2)
write.table(df, paste0(outTable, "Genes_Deciles_10.txt"), sep = "\t", quote = FALSE)

splitDF <- split(pbDF, rep(1:20, each = 879))
df <- data.frame(sapply(splitDF, function(subgroup) cor(subgroup$POOLED
  , subgroup$MEAN_EXPRESSION, method = "spearman")))
colnames(df) <- "SPEARMAN"
df <- round(df, 2)
write.table(df, paste0(outTable, "Genes_Deciles_20.txt"), sep = "\t", quote = FALSE)

splitDF <- split(pbDF, rep(1:5, each = 3516))
df <- data.frame(sapply(splitDF, function(subgroup) cor(subgroup$POOLED
  , subgroup$MEAN_EXPRESSION, method = "spearman")))
colnames(df) <- "SPEARMAN"
df <- round(df, 2)
write.table(df, paste0(outTable, "Genes_Deciles_5.txt"), sep = "\t", quote = FALSE)

