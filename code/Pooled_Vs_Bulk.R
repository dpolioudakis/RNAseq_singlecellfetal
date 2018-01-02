# Damon Polioudakis
# 2017-04-05
# Pooled versus bulk
# R 3.4.0
################################################################################

rm(list = ls())
sessionInfo()
set.seed(27)

# install.packages("devtools")
# library(devtools)
# install_github("satijalab/seurat", dependencies = TRUE, repos = 'http://cran.rstudio.com/')
require(Seurat)
require(biomaRt)
require(gProfileR)
require(reshape2)
require(eulerr)
require(gridExtra)

# load("../analysis/Pooled_Vs_Bulk_Workspace.RData")

# Seurat
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
dsExM <- as.matrix(centSO@raw.data)
centSO@raw.data <- NA
centSO@data <- NA
# load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# dsExM <- as.matrix(ssCentSO@raw.data)
# centSO <- ssCentSO

# Pollen
pnExDF <- read.csv("../pollen_2015/data/htseq/Exprs_HTSCexon.csv")
# # Picard Sequencing Statistics - bulk RNAseq
# picStatsBuDF <- read.table("../../kriegstein_2015/metadata/PicardToolsQC.csv")

# Fluidigm LT
flExDF <- read.csv("../C196-001_002_DP/data/htseq/merged/Exprs_HTSCexon.csv"
  , row.names = 1)
flMtDF <- read.table("../C196-001_002_DP/metadata/Compiled_Metadata_20160218.txt"
  , fill = TRUE, header = TRUE)

# Fluidigm HT
fhExDF <- read.csv(
  "../HT-003_DP/data/htseq/GRCh37.75_NoERCC/Exprs_HTSCexon_FtMm3e4.csv")

# Bulk
blExDF <- read.csv("../bulk_VZ_CP_from_ATAC/data/htseq/Exprs_HTSCexon.csv"
  , header = TRUE, row.names = 1)
# Metadata
blMtDF <- read.csv("../bulk_VZ_CP_from_ATAC/metadata/VZCP_sampleinfo.csv")

# biomart
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-01-05.csv", header = TRUE
  , fill = TRUE)

# # Markers from Seurat
# mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_Top20.txt"
#   , header = TRUE)

# Brain expressed genes from Vivek
brexpDF <- read.csv("../source/BrainExpressed_GTex.csv", header = FALSE)

## Variables
graphCodeTitle <- "Pooled_Vs_Bulk.R"
outGraph <- "../analysis/graphs/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outGraph), recursive = TRUE)
outTable <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outTable), recursive = TRUE)
outTableGpro <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler/Pooled_Vs_Bulk_gprofiler_"
dir.create("../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler", recursive = TRUE)
outTableGproIds <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler_goIDs/Pooled_Vs_Bulk_gprofiler_goIDs_"
dir.create("../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler_goIDs", recursive = TRUE)
outAnalysis <- "../analysis/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
dir.create(dirname(outAnalysis), recursive = TRUE)

## ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
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

Convert_Mixed_GeneSym_EnsID_To_EnsID <- function(ids){
  idx <- match(ids, bmDF$hgnc_symbol)
  ens <- bmDF$ensembl_gene_id[idx]
  ids[! is.na(ens)] <- as.character(ens[! is.na(ens)])
  return(ids)
}

# ensembl = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org")
# ensembl= useMart(biomart="ENSEMBL_MART_ENSEMBL"
#   , host="feb2014.archive.ensembl.org", path="/biomart/martservice"
#   , dataset="hsapiens_gene_ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
# listAttributes(ensembl)
# 
# listMarts(mart = NULL, host="feb2014.archive.ensembl.org", path="/biomart/martservice"
#   , port=80, includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE
#   , verbose = FALSE)
# 
# listAttributes(ensembl)[[1]][grep("gc", listAttributes(ensembl)[[1]])]

# Pools of different sizes correlations to bulk
# mnExDF: mean expression dataframe
# nSamp: Vector of sizes of pools
Pool_Size_Correlation <- function (mnExDF, nSamp) {
  sprCor <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    ssExDF <- mnExDF[ ,c(1:2, sample(3:ncol(mnExDF), i, replace = FALSE))]
    # Pool
    plExDF <- data.frame(ssExDF[1:2], POOLED = apply(ssExDF[3:ncol(ssExDF)], 1, sum))
    # # Read depth normalize pooled Drop-seq by number mapped to exons
    # rDep <- sum(plExDF$POOLED, na.rm = TRUE) / 10^6
    # Convert NAs to 0s
    plExDF[is.na(plExDF)] <- 0
    # # Read depth normalize
    # plExDF$POOLED <- plExDF$POOLED / rDep
    # Spearman correlation
    round(cor(plExDF$MEAN_BULK, plExDF$POOLED, method = "spearman"), 2)
  })
  sprDF <- data.frame(NUMBER_OF_CELLS = nSamp, SPEARMAN = sprCor)
  return(sprDF)
}

# Randomly sample and pool cells
# rNums: Columns to sample from
# nPool: Number of pools to split into
Pool_Cells <- function (exDF, nPool, rNums) {
  rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / nPool)))
  pExDF <- data.frame(lapply(rndmGroups
    , function (group) {apply(exDF[ ,group], 1, sum)}))
  print(head(pExDF, 10))
  return(pExDF)
}

# Subset to higher in bulk or higher in pooled
Subset_Two_SD <- function (dataset) {
  idx <- data.frame(
    HIGH_BULK = dataset < (mean(dataset) - 2*sd(dataset))
    , HIGH_POOL = dataset > (mean(dataset) + 2*sd(dataset))
  )
  return(idx)
}

# Plot gProfileR results
Graph_GO_Enrichment <- function (df, domain, depth, graphTitle) {
  ggDF <- subset(df, domain == domain & relative.depth == depth)
  ggDF <- ggDF[order(ggDF$p.value), ]
  # Top 30 by p-value
  ggDF <- ggDF[1:30, ]
  ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
  # Reverse order for ggplot so most significant at top of graph
  ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
  ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
  print(head(ggDF$term.name))
  # ggplot
  gg <- ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
    geom_bar(stat = "identity") +
    theme_set(theme_get() + theme(text = element_text(size = 11))) +
    coord_flip() +
    xlab(paste0(domain, " depth: ", depth)) +
    ylab("-log10(corrected p-value)") +
    ggtitle(graphTitle)
  # ggsave(paste0(outGraph, "AllGenes_ScatterSubsetHighBulk.pdf"))
  print(gg)
}

# CDS length from biomart, using longest CDS
Plot_CDS_Length <- function(geneSet, title) {
  gg <- ggplot(ggDF, aes_string(x = geneSet, y = "mnExDF.cds_length")) +
    geom_boxplot(fill = c("#e41a1c", "#377eb8", "#4daf4a")) +
    coord_cartesian(ylim = c(0, 5000)) +
    theme(text = element_text(size = 12)) +
    xlab("Gene subset") +
    ylab("CDS length") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"))
  return(gg)
}

# Percent GC from biomart
Plot_GC <- function(geneSet, title) {
  gg <- ggplot(ggDF, aes_string(x = geneSet, y = "mnExDF.percentage_gc_content")) +
    geom_boxplot(fill = c("#e41a1c", "#377eb8", "#4daf4a")) +
    coord_cartesian(ylim = c(0, 100)) +
    theme(text = element_text(size = 12)) +
    xlab("Gene subset") +
    ylab("Percent GC") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"))
  return(gg)
}

# Percentage of MT genes
Plot_MT_Percent <- function(geneSet, title) {
  # Percent
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "All other genes", ]
  genes <- df$HGNC_SYMBOL
  pct <- length(grep("^MT-", genes, value = TRUE)) / length(genes)
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Bulk high", ]
  genes <- df$HGNC_SYMBOL
  pct <- c(pct, length(grep("^MT-", genes, value = TRUE)) / length(genes))
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Pool high", ]
  genes <- df$HGNC_SYMBOL
  pct <- c(pct, length(grep("^MT-", genes, value = TRUE)) / length(genes))
  pct <- pct*100
  # Format for ggplot
  ggDF <- data.frame(SUBSET = c("All other genes", "Bulk high", "Pool high")
    , PERCENT_MT = pct)
  # ggplot
  gg <- ggplot(ggDF, aes(x = SUBSET, y = PERCENT_MT)) +
    geom_bar(stat = "identity", fill = c("#e41a1c", "#377eb8", "#4daf4a")) +
    coord_cartesian(ylim = c(0, 100)) +
    theme(text = element_text(size = 12)) +
    xlab("Gene subset") +
    ylab("Percent MT") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"))
  return(gg)
}
################################################################################

### Subset drop-seq raw counts matrix to cells that passed QC

dsExM <- dsExM[ ,colnames(dsExM) %in% row.names(centSO@meta.data)]
################################################################################

### Remove ERCCs and STAR stats from tail and move gene names to row.names

# Remove ERCCs and STAR stats from tail
# Bulk
tail(blExDF, 10)
blExDF <- head(blExDF, -97)
tail(blExDF, 5)
# Pollen
tail(pnExDF, 10)
pnExDF <- head(pnExDF, -97)
tail(pnExDF, 5)
# Fluidigm LT
tail(flExDF, 10)
flExDF <- head(flExDF, -97)
tail(flExDF, 5)
# Fluidigm HT
tail(fhExDF, 10)
fhExDF <- head(fhExDF, -5)
tail(fhExDF, 5)

# Move gene names
# Pollen
row.names(pnExDF) <- pnExDF$X
pnExDF <- pnExDF[ ,-1]
# Fluidigm HT
row.names(fhExDF) <- fhExDF$X
fhExDF <- fhExDF[ ,-1]
################################################################################

### Subset to GZ

# Bulk
# Check GZ sample IDs
blMtDF$RNAID[blMtDF$ExpCondition == "VZ"]
head(blExDF[ ,blMtDF$ExpCondition == "VZ"])
# Subset
gzBlExDF <- blExDF[ ,blMtDF$ExpCondition == "VZ"]

# Drop-seq
# # Check GZ sample IDs
# table(grepl("VZ", colnames(dsExM)))
# Subset
# gzDsExM <- dsExM[ ,grep("VZ", colnames(dsExM))]
ids <- row.names(centSO@meta.data)[centSO@meta.data$REGION == "GZ"]
gzDsExM <- dsExM[, colnames(dsExM) %in% ids]

# Fluidigm HT
gzFhExDF <- fhExDF[ ,grep(
  "COL01|COL02|COL03|COL04|COL05|COL06|COL07|COL08|COL09|COL10"
  , colnames(fhExDF))]
################################################################################

### Format

# Clean chr number of ens IDs for bulk, Pollen
# Bulk
row.names(blExDF) <- gsub("\\.[0-9]*", "", row.names(blExDF))
# Bulk GZ
row.names(gzBlExDF) <- gsub("\\.[0-9]*", "", row.names(gzBlExDF))
# Pollen
row.names(pnExDF) <- gsub("\\.[0-9]*", "", row.names(pnExDF))
# Fluidigm LT
row.names(flExDF) <- gsub("\\.[0-9]*", "", row.names(flExDF))
# Fluidigm HT
row.names(fhExDF) <- gsub("\\.[0-9]*", "", row.names(fhExDF))
# Fluidigm HT GZ
row.names(gzFhExDF) <- gsub("\\.[0-9]*", "", row.names(gzFhExDF))

# Convert gene symbols to Ensembl IDs for Drop-seq
# # All cells
# dim(dsExM)
# ensDF <- QueryBiomaRt(row.names(dsExM), "hgnc_symbol"
#   , c("hgnc_symbol", "ensembl_gene_id"))
# dsExDF <- merge(ensDF, dsExM, by.x = "hgnc_symbol", by.y = "row.names")
# dim(dsExDF)
# # GZ
# # Convert gene symbols to Ensembl IDs for Drop-seq
# dim(gzDsExM)
# ensDF <- QueryBiomaRt(row.names(gzDsExM), "hgnc_symbol"
#   , c("hgnc_symbol", "ensembl_gene_id"))
# gzDsExDF <- merge(ensDF, gzDsExM, by.x = "hgnc_symbol", by.y = "row.names")
# dim(gzDsExDF)
dsExDF <- data.frame(dsExM)
dsExDF$ensembl_gene_id <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(dsExDF))
gzDsExDF <- data.frame(gzDsExM)
gzDsExDF$ensembl_gene_id <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(gzDsExDF))
################################################################################

### Pools of different sizes and read depth normalize

# Read depth normalize bulk by number mapped to exons
# All cells
rDep <- (apply(blExDF, 2, sum) / 10^6)
buCpmDF <- blExDF / rDep
# Bulk mean CPM
mnBuExDF <- data.frame(MEAN_BULK = apply(buCpmDF, 1, mean))
# GZ
rDep <- (apply(gzBlExDF, 2, sum) / 10^6)
gzBuCpmDF <- gzBlExDF / rDep
# Bulk mean CPM
mnGzBuExDF <- data.frame(MEAN_BULK = apply(gzBuCpmDF, 1, mean))

## Pools of different sizes correlated to bulk

# Drop-seq
# All cells
# Merge bulk and pooled means, include genes not in both
mnBuDsExDF <- merge(mnBuExDF, dsExDF, by.x = "row.names"
  , by.y = "ensembl_gene_id", all = TRUE)
names(mnBuDsExDF)[1] <- "ensembl_gene_id"
mnBuDsExDF <- mnBuDsExDF[ ,! names(mnBuDsExDF) == "hgnc_symbol"]
# Convert NAs to 0s
mnBuDsExDF[is.na(mnBuDsExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
  , seq(1000, (ncol(dsExM)-3), 500))
sprDsDF <- Pool_Size_Correlation(mnBuDsExDF, nSamp)
# GZ
# Merge bulk and pooled means, include genes not in both
mnGzBuDsExDF <- merge(mnGzBuExDF, gzDsExDF, by.x = "row.names"
  , by.y = "ensembl_gene_id", all = TRUE)
names(mnGzBuDsExDF)[1] <- "ensembl_gene_id"
mnGzBuDsExDF <- mnGzBuDsExDF[ ,! names(mnGzBuDsExDF) == "hgnc_symbol"]
# Convert NAs to 0s
mnGzBuDsExDF[is.na(mnGzBuDsExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
  , seq(1000, (ncol(gzDsExM)-3), 500))
sprGzDsDF <- Pool_Size_Correlation(mnGzBuDsExDF, nSamp)

# Pollen
# Merge bulk and pooled means, include genes not in both
mnBuPnExDF <- merge(mnGzBuExDF, pnExDF, by = "row.names", all = TRUE)
names(mnBuPnExDF)[1] <- "ensembl_gene_id"
# Convert NAs to 0s
mnBuPnExDF[is.na(mnBuPnExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(pnExDF)-2), 100))
sprPnDF <- Pool_Size_Correlation(mnBuPnExDF, nSamp)

# Fluidigm HT
# All cells
# Merge bulk and pooled means, include genes not in both
mnBuFhExDF <- merge(mnBuExDF, fhExDF, by = "row.names", all = TRUE)
names(mnBuFhExDF)[1] <- "ensembl_gene_id"
# Convert NAs to 0s
mnBuFhExDF[is.na(mnBuFhExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(fhExDF)-2), 100))
sprFhDF <- Pool_Size_Correlation(mnBuFhExDF, nSamp)
# GZ
# Merge bulk and pooled means, include genes not in both
mnGzBuFhExDF <- merge(mnGzBuExDF, gzFhExDF, by = "row.names", all = TRUE)
names(mnGzBuFhExDF)[1] <- "ensembl_gene_id"
# Convert NAs to 0s
mnGzBuFhExDF[is.na(mnGzBuFhExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(gzFhExDF)-2), 100))
sprGzFhDF <- Pool_Size_Correlation(mnGzBuFhExDF, nSamp)

# Fluidigm LT
# All cells
# Merge bulk and pooled means, include genes not in both
mnBuFlExDF <- merge(mnBuExDF, flExDF, by = "row.names", all = TRUE)
names(mnBuFlExDF)[1] <- "ensembl_gene_id"
# Convert NAs to 0s
mnBuFlExDF[is.na(mnBuFlExDF)] <- 0
# Pool sizes
nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(flExDF)-2), 50))
sprFlDF <- Pool_Size_Correlation(mnBuFlExDF, nSamp)

## Plot as grid
# All cells
sprDsDF$DATASET <- "DROP-SEQ"
sprFhDF$DATASET <- "FLUIDIGM_HT"
sprFlDF$DATASET <- "FLUIDIGM_LT"
ggDF <- rbind(sprDsDF, sprFhDF, sprFlDF)
ggplot(ggDF, aes(x = NUMBER_OF_CELLS, y = SPEARMAN)) +
  facet_wrap( ~ DATASET, ncol = 2, scales = "free_x") +
  geom_line() +
  xlab("Number of cells in pool") +
  ylab("Spearman") +
  ylim(0, 1) +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCorrelation of different size pools to bulk"
    , "\nGZ and CP fetal human brain"
    , "\n"
  ))
ggsave(paste0(outGraph, "Pool_Size_Spearman.pdf"))
# GZ
sprGzDsDF$DATASET <- "DROP-SEQ"
sprPnDF$DATASET <- "POLLEN_2015"
sprGzFhDF$DATASET <- "FLUIDIGM_HT"
ggDF <- rbind(sprGzDsDF, sprPnDF, sprGzFhDF)
ggplot(ggDF, aes(x = NUMBER_OF_CELLS, y = SPEARMAN)) +
  facet_wrap( ~ DATASET, ncol = 2, scales = "free_x") +
  geom_line() +
  xlab("Number of cells in pool") +
  ylab("Spearman") +
  ylim(0, 1) +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCorrelation of different size pools to bulk"
    , "\nGZ fetal human brain"
    , "\n"
  ))
ggsave(paste0(outGraph, "Pool_Size_Spearman_GZ.pdf"))

# Number of genes detected as pool size increases
Number_Genes_Detected <- function (exDF, nSamp, datasetName) {
  # Genes with >= 1 count (raw counts)
  print("Genes with >= 1 count (raw counts)")
  tfExDF <- exDF >= 1
  nGene1 <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    ssExDF <- tfExDF[ ,c(sample(1:ncol(tfExDF), i, replace = FALSE))]
    # rowSums throws error if only 1 column (nSamp = 1)
    if (i == 1) {
      print("nSamp = 1")
      sum(ssExDF >= 1)
    }
    else {
      # Sum genes with >= 1 count
      rSums <- rowSums(ssExDF)
      sum(rSums >= 1)  
    }
  })
  # Genes with >= 2 count (raw counts)
  print("Genes with >= 2 count (raw counts)")
  tfExDF <- exDF >= 2
  nGene2 <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    ssExDF <- tfExDF[ ,c(sample(1:ncol(tfExDF), i, replace = FALSE))]
    # rowSums throws error if only 1 column (nSamp = 1)
    if (i == 1) {
      print("nSamp = 1")
      sum(ssExDF >= 1)
    }
    else {
      # Sum genes with >= 1 count
      rSums <- rowSums(ssExDF)
      sum(rSums >= 1)  
    }
  })
  # Genes with >= 3 count (raw counts)
  print("Genes with >= 3 count (raw counts)")
  tfexDF <- exDF >= 3
  nGene3 <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    ssExDF <- tfexDF[ ,c(sample(1:ncol(tfexDF), i, replace = FALSE))]
    # rowSums throws error if only 1 column (nSamp = 1)
    if (i == 1) {
      print("nSamp = 1")
      sum(ssExDF >= 1)
    }
    else {
      # Sum genes with >= 3 count
      rSums <- rowSums(ssExDF)
      sum(rSums >= 1)  
    }
  })
  # Genes with >= 10 count (raw counts)
  print("Genes with >= 10 count (raw counts)")
  tfexDF <- exDF >= 10
  nGene10 <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    ssExDF <- tfexDF[ ,c(sample(1:ncol(tfexDF), i, replace = FALSE))]
    # rowSums throws error if only 1 column (nSamp = 1)
    if (i == 1) {
      print("nSamp = 1")
      sum(ssExDF >= 1)
    }
    else {
      # Sum genes with >= 10 count
      rSums <- rowSums(ssExDF)
      sum(rSums >= 1)  
    }
  })
  nGeneDF <- data.frame(
    NUMBER_OF_CELLS = nSamp
    , DATASET = datasetName
    , NUMBER_GENES_DETECTED_1 = nGene1
    , NUMBER_GENES_DETECTED_2 = nGene2
    , NUMBER_GENES_DETECTED_3 = nGene3
    , NUMBER_GENES_DETECTED_10 = nGene10)
  return(nGeneDF)
}

### Number of genes detected as pool size increases
ngdDF <- rbind(
  # Drop-seq
  Number_Genes_Detected(dsExM
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
      , seq(1000, (ncol(dsExM)), 500))
    , "DROPSEQ")
  # Fluidigm HT
  , Number_Genes_Detected(fhExDF
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(fhExDF)), 100))
    , "FLUIDIGM_HT")
  # Fluidigm LT
  , Number_Genes_Detected(flExDF
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(flExDF)), 50))
    , "FLUIDIGM_LT")
  # Drop-seq GZ
  , Number_Genes_Detected(gzDsExM
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
      , seq(1000, (ncol(gzDsExM)), 500))
    , "DROPSEQ_GZ")
  # Fluidigm HT GZ
  , Number_Genes_Detected(gzFhExDF
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(gzFhExDF)), 100))
    , "FLUIDIGM_HT_GZ")
  , Number_Genes_Detected(pnExDF
    , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(pnExDF)), 100))
    , "POLLEN")
)

ggDF <- melt(ngdDF, measure.vars = c("NUMBER_GENES_DETECTED_1"
  , "NUMBER_GENES_DETECTED_2", "NUMBER_GENES_DETECTED_3"
  , "NUMBER_GENES_DETECTED_10"))

## Plot as grid
ggplot(ggDF, aes(x = NUMBER_OF_CELLS, y = value, col = variable)) +
  facet_wrap( ~ DATASET, ncol = 3, scales = "free_x") +
  geom_line() +
  xlab("Number of cells in pool") +
  ylab("Number of genes detected") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nNumber genes detected as pool size increases"
    , "\nGZ and CP fetal human brain"
    , "\n"
  ))
ggsave(paste0(outGraph, "Pool_Size_Genes_Detected.pdf"), height = 7, width = 11)
################################################################################

### Pool scRNA-seq and read depth normalize

# Drop-seq randomly split into pooled groups
# All cells - 310 pooled groups of 100 cells each
ncol(dsExM)
nPool <- 310
rNums <- sample(1:31000, 31000, replace = FALSE)
pDsExDF <- Pool_Cells(dsExM, nPool, rNums)
# GZ - 150 pooled groups of 100 cells each
ncol(gzDsExM)
nPool <- 150
rNums <- sample(1:15000, 15000, replace = FALSE)
pGzDsExDF <- Pool_Cells(gzDsExM, nPool, rNums)

# # For testing
# ncol(dsExM)
# nPool <- 20
# rNums <- sample(1:2000, 2000, replace = FALSE)
# pDsExDF <- Pool_Cells(dsExM, nPool, rNums)
# # GZ - 34 pooled groups of 100 cells each
# ncol(gzDsExM)
# nPool <- 10
# rNums <- sample(1:1000, 1000, replace = FALSE)
# pGzDsExDF <- Pool_Cells(gzDsExM, nPool, rNums)

# Pollen - randomly split into 10 pooled groups of 39 cells each
ncol(pnExDF)
nPool <- 10
rNums <- sample(1:390, 390, replace = FALSE)
pPnExDF <- Pool_Cells(pnExDF, nPool, rNums)

# Fluidigm HT - randomly split into pooled groups
# All cells - 12 pooled groups of 48 cells each
ncol(fhExDF)
nPool <- 12
rNums <- sample(1:576, 576, replace = FALSE)
pFhExDF <- Pool_Cells(fhExDF, nPool, rNums)
# GZ - 6 pooled groups of 48 cells each
ncol(gzFhExDF)
nPool <- 6
rNums <- sample(1:288, 288, replace = FALSE)
pGzFhExDF <- Pool_Cells(gzFhExDF, nPool, rNums)

# Fluidigm LT - randomly split into pooled groups
# All cells - 6 pooled groups of 30 cells each
ncol(flExDF)
nPool <- 6
rNums <- sample(1:180, 180, replace = FALSE)
pFlExDF <- Pool_Cells(flExDF, nPool, rNums)

## Read depth normalize

# Read depth normalize pooled Drop-seq by number mapped to exons
# All cells
rDep <- (apply(pDsExDF, 2, sum) / 10^6)
pDsCpmDF <- pDsExDF / rDep
# GZ
rDep <- (apply(pGzDsExDF, 2, sum) / 10^6)
pGzDsCpmDF <- pGzDsExDF / rDep

# Read depth normalize pooled Pollen by number mapped to exons
rDep <- (apply(pPnExDF, 2, sum) / 10^6)
pPnCpmDF <- pPnExDF / rDep

# Read depth normalize pooled Fluidigm HT by number mapped to exons
# All cells
rDep <- (apply(pFhExDF, 2, sum) / 10^6)
pFhCpmDF <- pFhExDF / rDep
# GZ
rDep <- (apply(pGzFhExDF, 2, sum) / 10^6)
pGzFhCpmDF <- pGzFhExDF / rDep

# Read depth normalize pooled Fluidigm LT by number mapped to exons
# All cells
rDep <- (apply(pFlExDF, 2, sum) / 10^6)
pFlCpmDF <- pFlExDF / rDep
################################################################################

### Mean counts and merge datasets to align genes

## Mean counts

# # Mean counts for bulk RNAseq
# mnBuExDF <- data.frame(MEAN_BULK = apply(buCpmDF, 1, mean))

# Mean counts for pooled Drop-seq groups
# All cells
mnPoDsExDF <- data.frame(MEAN_POOLED_DROPSEQ = apply(pDsCpmDF, 1, mean))
# GZ
mnPoGzDsExDF <- data.frame(MEAN_POOLED_DROPSEQ_GZ = apply(pGzDsCpmDF, 1, mean))

# Mean counts for pooled Pollen groups
mnPoPnExDF <- data.frame(MEAN_POOLED_POLLEN = apply(pPnCpmDF, 1, mean))

# Mean counts for pooled Fluidigm HT groups
# All cells
mnPoFhExDF <- data.frame(MEAN_POOLED_FLUIDIGM_HT = apply(pFhCpmDF, 1, mean))
# GZ
mnPoGzFhExDF <- data.frame(MEAN_POOLED_FLUIDIGM_HT_GZ = apply(pGzFhCpmDF, 1, mean))

# Mean counts for pooled Fluidigm LT groups
# All cells
mnPoFlExDF <- data.frame(MEAN_POOLED_FLUIDIGM_LT = apply(pFlCpmDF, 1, mean))

# Convert gene symbols to Ensembl IDs for Drop-seq

# # All cells
# dim(mnPoDsExDF)
# ensDF <- QueryBiomaRt(row.names(mnPoDsExDF), "hgnc_symbol"
#   , c("hgnc_symbol", "ensembl_gene_id"))
# mnPoDsExDF <- merge(ensDF, mnPoDsExDF, by.x = "hgnc_symbol", by.y = "row.names")
# dim(mnPoDsExDF)
# # GZ
# dim(mnPoGzDsExDF)
# ensDF <- QueryBiomaRt(row.names(mnPoGzDsExDF), "hgnc_symbol"
#   , c("hgnc_symbol", "ensembl_gene_id"))
# mnPoGzDsExDF <- merge(ensDF, mnPoGzDsExDF, by.x = "hgnc_symbol", by.y = "row.names")
# dim(mnPoGzDsExDF)
row.names(mnPoDsExDF) <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(mnPoDsExDF))
mnPoDsExDF$ensembl_gene_id <- row.names(mnPoDsExDF)
row.names(mnPoGzDsExDF) <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(mnPoGzDsExDF))
mnPoGzDsExDF$ensembl_gene_id <- row.names(mnPoGzDsExDF)

## Combine

# Merge bulk and pooled means, include genes not in both
# Bulk GZ
names(mnGzBuExDF) <- "MEAN_BULK_GZ"
mnExDF <- merge(mnBuExDF, mnGzBuExDF, by.x = "row.names"
  , by.y = "row.names", all = TRUE)
# DS all cells
mnExDF <- merge(mnExDF, mnPoDsExDF, by.x = "Row.names"
  , by.y = "ensembl_gene_id", all = TRUE)
# DS GZ
mnExDF <- merge(mnExDF, mnPoGzDsExDF, by.x = "Row.names"
  , by.y = "ensembl_gene_id", all = TRUE)
# Pollen
mnExDF <- merge(mnExDF, mnPoPnExDF, by.x = "Row.names", by.y = "row.names"
  , all = TRUE)
# Fluidigm HT all cells
mnExDF <- merge(mnExDF, mnPoFhExDF, by.x = "Row.names", by.y = "row.names"
  , all = TRUE)
# Fluidigm HT GZ
mnExDF <- merge(mnExDF, mnPoGzFhExDF, by.x = "Row.names", by.y = "row.names"
  , all = TRUE)
# Fluidigm LT all cells
mnExDF <- merge(mnExDF, mnPoFlExDF, by.x = "Row.names", by.y = "row.names"
  , all = TRUE)
# Convert NAs to 0s
mnExDF[is.na(mnExDF)] <- 0
head(mnExDF)

# Convert ensembl ID to gene symbol and subset to ensembl IDs found in biomaRt
df <- QueryBiomaRt(mnExDF$Row.names, "ensembl_gene_id"
  , c("ensembl_gene_id", "hgnc_symbol"))
mnExDF <- merge(df, mnExDF, by.x = "ensembl_gene_id", by.y = "Row.names")
# Add gene information
df <- QueryBiomaRt(mnExDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "gene_biotype", "percentage_gc_content"))
mnExDF <- merge(mnExDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(mnExDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "cds_length"))
# Take longest CDS
df <- aggregate(cds_length ~ ensembl_gene_id, data = df, max)
mnExDF <- merge(mnExDF, df, by = "ensembl_gene_id", all = TRUE)
str(mnExDF)
head(mnExDF)
# # Convert ensembl ID to gene symbol
# df <- EnsemblToGeneSymbol(mnExDF$Row.names)
# mnExDF$hgnc_symbol <- df$hgnc_symbol[match(mnExDF$Row.names, df$ensembl_gene_id)]
# head(mnExDF)
# Convert NAs to 0s
mnExDF[is.na(mnExDF)] <- 0
head(mnExDF)

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))

# Format for ggplot
df1 <- melt(mnExDF[c("MEAN_BULK", "MEAN_BULK_GZ", "MEAN_POOLED_DROPSEQ"
  , "MEAN_POOLED_FLUIDIGM_HT", "MEAN_POOLED_FLUIDIGM_LT")]
  , measure.vars = c("MEAN_POOLED_DROPSEQ"
  , "MEAN_POOLED_FLUIDIGM_HT", "MEAN_POOLED_FLUIDIGM_LT"))
df2 <- melt(mnExDF[c("MEAN_BULK", "MEAN_BULK_GZ", "MEAN_POOLED_DROPSEQ_GZ"
  , "MEAN_POOLED_POLLEN", "MEAN_POOLED_FLUIDIGM_HT_GZ")]
  , measure.vars = c("MEAN_POOLED_DROPSEQ_GZ"
  , "MEAN_POOLED_POLLEN", "MEAN_POOLED_FLUIDIGM_HT_GZ"))
df2$MEAN_BULK <- df2$MEAN_BULK_GZ 
ggDF <- rbind(df1, df2)
# Spearman
sprL <- lapply(split(ggDF, ggDF$variable), function(df) {
  cor(df$MEAN_BULK, df$value, method = "spearman")})
sprL <- lapply(sprL, round, 2)
ggDF$value <- log(ggDF$value + 1, 2)
ggDF$MEAN_BULK <- log(ggDF$MEAN_BULK + 1, 2)
labelsDatasetSpr <- c(paste0(unique(ggDF$variable), "\nSpearman: ", sprL))
names(labelsDatasetSpr) <- unique(ggDF$variable)
# ggplot
ggplot(ggDF, aes(x = value, y = MEAN_BULK)) +
  facet_wrap("variable", ncol = 2
    , labeller = labeller(variable = labelsDatasetSpr)) +
  geom_point(alpha = 0.15, shape = 1, size = 0.25) +
  stat_smooth() +
  ylab("Bulk: log2(mean CPM + 1)") +
  xlab("Pooled: log2(mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPooled scRNA-seq vs Bulk RNA-seq"
    , "\n"
    , "\nMean of CPM across samples and pools"))
ggsave(paste0(outGraph, "Scatter.png"), dpi = 600, width = 6.5, height = 12)

## Drop-seq union genes versus intersection genes
# Intersection (only genes in both)
mnBothDF <- merge(mnBuExDF, mnPoDsExDF, by.x = "row.names"
  , by.y = "ensembl_gene_id")
# Convert NAs to 0s
mnBothDF[is.na(mnBothDF)] <- 0
# Spearman correlation
# Union
sprUn <- round(cor(mnExDF$MEAN_BULK, mnExDF$MEAN_POOLED_DROPSEQ
  , method = "spearman"), 2)
sprUn
# Intersection
sprIn <- round(cor(mnBothDF$MEAN_BULK, mnBothDF$MEAN_POOLED
  , method = "spearman"), 2)
sprIn
# Format for ggplot2
df1 <- data.frame(Pooled = log(mnExDF$MEAN_BULK + 1, 2)
  , Bulk = log(mnExDF$MEAN_POOLED_DROPSEQ + 1, 2))
df1$DATASET <- "Union genes"
df2 <- data.frame(Pooled = log(mnBothDF$MEAN_BULK + 1, 2)
  , Bulk = log(mnBothDF$MEAN_POOLED + 1, 2))
df2$DATASET <- "Intersection genes"
ggDF <- rbind(df1, df2)
# Add spearman to labels
sprL <- c(sprUn, sprIn)
labelsDatasetSpr <- c(paste0(unique(ggDF$DATASET), "\nSpearman: ", sprL))
names(labelsDatasetSpr) <- unique(ggDF$DATASET)
# ggplot
ggplot(ggDF, aes(x = Pooled, y = Bulk, col = DATASET)) +
  facet_wrap("DATASET", ncol = 2
    , labeller = labeller(DATASET = labelsDatasetSpr)) +
  geom_point(alpha = 0.5, shape = 1, size = 0.25) +
  stat_smooth(col = "black") +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  theme(legend.position = "none") +
  ylab("Bulk: log2(Mean CPM + 1)") +
  xlab("Pooled: log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPooled Drop-seq vs Bulk RNAseq"
    , "\n"
    , "\nIntersection (genes only detected in both) or union (genes in either)"
    , "\nHuman Fetal Brain CP + GZ"
    , "\nMean of CPM across samples and pools"))
ggsave(paste0(outGraph, "ScatterIntersectionUnion.png"), width = 8, height = 5
  , dpi = 600)
################################################################################

### Subset to higher in pooled or higher in bulk

# log2(mean CPM + 1)
l2MnCpmDF <- data.frame(ENSEMBL_ID = mnExDF$ensembl_gene_id
  , HGNC_SYMBOL = mnExDF$hgnc_symbol
  , MEAN_POOLED_DROPSEQ = log(mnExDF$MEAN_POOLED_DROPSEQ + 1, 2)
  , MEAN_POOLED_DROPSEQ_GZ = log(mnExDF$MEAN_POOLED_DROPSEQ_GZ + 1, 2)
  , MEAN_POOLED_POLLEN = log(mnExDF$MEAN_POOLED_POLLEN + 1, 2)
  , MEAN_POOLED_FLUIDIGM_HT = log(mnExDF$MEAN_POOLED_FLUIDIGM_HT + 1, 2)
  , MEAN_POOLED_FLUIDIGM_HT_GZ = log(mnExDF$MEAN_POOLED_FLUIDIGM_HT_GZ + 1, 2)
  , MEAN_POOLED_FLUIDIGM_LT = log(mnExDF$MEAN_POOLED_FLUIDIGM_LT + 1, 2)
  , MEAN_BULK = log(mnExDF$MEAN_BULK + 1, 2)
  , MEAN_BULK_GZ = log(mnExDF$MEAN_BULK_GZ + 1, 2))

# Ratio pooled vs bulk
l2MnCpmDF$RATIO_DROPSEQ <- l2MnCpmDF$MEAN_POOLED_DROPSEQ - l2MnCpmDF$MEAN_BULK
l2MnCpmDF$RATIO_DROPSEQ_GZ <- l2MnCpmDF$MEAN_POOLED_DROPSEQ_GZ - l2MnCpmDF$MEAN_BULK_GZ
l2MnCpmDF$RATIO_POLLEN <- l2MnCpmDF$MEAN_POOLED_POLLEN - l2MnCpmDF$MEAN_BULK_GZ
l2MnCpmDF$RATIO_FLUIDIGM_HT <- l2MnCpmDF$MEAN_POOLED_FLUIDIGM_HT - l2MnCpmDF$MEAN_BULK
l2MnCpmDF$RATIO_FLUIDIGM_HT_GZ <- l2MnCpmDF$MEAN_POOLED_FLUIDIGM_HT_GZ - l2MnCpmDF$MEAN_BULK_GZ
l2MnCpmDF$RATIO_FLUIDIGM_LT <- l2MnCpmDF$MEAN_POOLED_FLUIDIGM_LT - l2MnCpmDF$MEAN_BULK

# # Subset to genes not detected in pooled and top 10% quantile in bulk
# # ssP0b20df <- subset(l2MnCpmDF, BULK >= quantile(BULK, 0.9) & POOLED == 0)
# ssPlBhdf <- subset(l2MnCpmDF, BULK >= quantile(BULK, 0.6) & POOLED <= quantile(POOLED, 0.4))
# ssPlBhdf <- subset(l2MnCpmDF, RATIO >= quantile(RATIO, 0.1, na.rm = TRUE))
# # ssPlBhdf <- subset(l2MnCpmDF, RATIO < (mean(RATIO) - 2*sd(RATIO)) & BULK >= quantile(BULK, 0.75))

# Subset 
ssl2MnCpmDF <- l2MnCpmDF
ssl2MnCpmDF$SUBSET_DROPSEQ <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_DROPSEQ)
ssl2MnCpmDF$SUBSET_DROPSEQ[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_DROPSEQ[idxDF$HIGH_BULK] <- "Bulk high"

ssl2MnCpmDF$SUBSET_DROPSEQ_0 <- "All other genes"
ssl2MnCpmDF$SUBSET_DROPSEQ_0[ssl2MnCpmDF$MEAN_BULK == 0 & ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high"] <- "Pool high"
ssl2MnCpmDF$SUBSET_DROPSEQ_0[ssl2MnCpmDF$MEAN_POOLED_DROPSEQ == 0 & ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high"] <- "Bulk high"

ssl2MnCpmDF$SUBSET_FLUIDIGM_HT <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_FLUIDIGM_HT)
ssl2MnCpmDF$SUBSET_FLUIDIGM_HT[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_FLUIDIGM_HT[idxDF$HIGH_BULK] <- "Bulk high"

ssl2MnCpmDF$SUBSET_FLUIDIGM_LT <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_FLUIDIGM_LT)
ssl2MnCpmDF$SUBSET_FLUIDIGM_LT[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_FLUIDIGM_LT[idxDF$HIGH_BULK] <- "Bulk high"

ssl2MnCpmDF$SUBSET_DROPSEQ_GZ <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_DROPSEQ_GZ)
ssl2MnCpmDF$SUBSET_DROPSEQ_GZ[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_DROPSEQ_GZ[idxDF$HIGH_BULK] <- "Bulk high"

ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_FLUIDIGM_HT_GZ)
ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ[idxDF$HIGH_BULK] <- "Bulk high"

ssl2MnCpmDF$SUBSET_POLLEN <- "All other genes"
idxDF <- Subset_Two_SD(l2MnCpmDF$RATIO_POLLEN)
ssl2MnCpmDF$SUBSET_POLLEN[idxDF$HIGH_POOL] <- "Pool high"
ssl2MnCpmDF$SUBSET_POLLEN[idxDF$HIGH_BULK] <- "Bulk high"

## Table of number of genes in each subset
df <- cbind(data.frame(table(ssl2MnCpmDF$SUBSET_DROPSEQ))
  , data.frame(table(ssl2MnCpmDF$SUBSET_DROPSEQ_0))[ ,2]
  , data.frame(table(ssl2MnCpmDF$SUBSET_FLUIDIGM_HT))[ ,2]
  , data.frame(table(ssl2MnCpmDF$SUBSET_FLUIDIGM_LT))[ ,2]
  , data.frame(table(ssl2MnCpmDF$SUBSET_DROPSEQ_GZ))[ ,2]
  , data.frame(table(ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ))[ ,2]
  , data.frame(table(ssl2MnCpmDF$SUBSET_POLLEN))[ ,2])
colnames(df) <- c("GENE_SET", "DROPSEQ", "DROPSEQ_0", "FLUIDIGM_HT"
  , "FLUIDIGM_LT", "DROP_SEQ_GZ", "FLUIDIGM_HT_GZ", "POLLEN")
write.csv(df, file = paste0(outTable, "_NumberGenesSubset.csv"), quote = FALSE)

## Plot subset scatter plots
# Format for ggplot
df1 <- ssl2MnCpmDF[c("MEAN_BULK", "MEAN_POOLED_DROPSEQ", "SUBSET_DROPSEQ")]
names(df1)[2:3] <- c("MEAN_POOLED", "SUBSET")
df1$DATASET <- "DROPSEQ"
df7 <- ssl2MnCpmDF[c("MEAN_BULK", "MEAN_POOLED_DROPSEQ", "SUBSET_DROPSEQ_0")]
names(df7)[2:3] <- c("MEAN_POOLED", "SUBSET")
df7$DATASET <- "DROPSEQ_0"
df2 <- ssl2MnCpmDF[c("MEAN_BULK", "MEAN_POOLED_FLUIDIGM_HT", "SUBSET_FLUIDIGM_HT")]
names(df2)[2:3] <- c("MEAN_POOLED", "SUBSET")
df2$DATASET <- "FLUIDIGM_HT"
df3 <- ssl2MnCpmDF[c("MEAN_BULK", "MEAN_POOLED_FLUIDIGM_LT", "SUBSET_FLUIDIGM_LT")]
names(df3)[2:3] <- c("MEAN_POOLED", "SUBSET")
df3$DATASET <- "FLUIDIGM_LT"
df4 <- ssl2MnCpmDF[c("MEAN_BULK_GZ", "MEAN_POOLED_DROPSEQ_GZ", "SUBSET_DROPSEQ_GZ")]
names(df4)[1:3] <- c("MEAN_BULK", "MEAN_POOLED", "SUBSET")
df4$DATASET <- "DROPSEQ_GZ"
df5 <- ssl2MnCpmDF[c("MEAN_BULK_GZ", "MEAN_POOLED_FLUIDIGM_HT_GZ", "SUBSET_FLUIDIGM_HT_GZ")]
names(df5)[1:3] <- c("MEAN_BULK", "MEAN_POOLED", "SUBSET")
df5$DATASET <- "FLUIDIGM_HT_GZ"
df6 <- ssl2MnCpmDF[c("MEAN_BULK_GZ", "MEAN_POOLED_POLLEN", "SUBSET_POLLEN")]
names(df6)[1:3] <- c("MEAN_BULK", "MEAN_POOLED", "SUBSET")
df6$DATASET <- "POLLEN"
# Rbind
ggDF <- rbind(df1, df7, df2, df3, df4, df5, df6)
# ggplot
ggplot(ggDF, aes(x = MEAN_POOLED, y = MEAN_BULK)) +
  facet_wrap("DATASET", ncol = 2) +
  geom_point(alpha = 0.15, shape = 1, size = 0.25, aes(col = SUBSET)) +
  scale_colour_manual(values = c("#bdbdbd","#1f78b4", "#e31a1c")) +
  ylab("Bulk: log2(Mean CPM + 1)") +
  xlab("Pooled: log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSubsetting pooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain GZ and CP"
    , "\nMean of CPM across samples and pools"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"))
ggsave(paste0(outGraph, "ScatterSubset.png"), dpi = 600, width = 8, height = 12)

## Plot subset density plot

ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  facet_grid(variable ~ DATASET) +
  geom_density(aes(col = SUBSET)) +
  scale_colour_manual(values = c("grey","#1f78b4", "#e31a1c")) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Density") +
  xlab("log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSubsetting pooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain GZ and CP"
    , "\nMean of CPM across samples and pools"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"))
ggsave(paste0(outGraph, "DensitySubset.pdf"), width = 14, height = 8)

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))
################################################################################

### Gene Biotype frequency

ssl2MnCpmDF <- data.frame(ssl2MnCpmDF, gene_biotype = mnExDF$gene_biotype)
ggDF <- ssl2MnCpmDF
# Remove TR_ and IG_ biotypes
ggDF <- ggDF[-grep("TR_|IG_", ggDF$gene_biotype), ]
unique(ggDF$gene_biotype)
ggDF$gene_biotype <- factor(ggDF$gene_biotype, levels = unique(ggDF$gene_biotype))
# Frequency
ggDF <- data.frame(
  data.frame(table(ggDF$gene_biotype))
  , DROPSEQ_HIGH = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_DROPSEQ == "Pool high"]))$Freq
  # , DROPSEQ_ALL_OTHER = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_DROPSEQ == "All other genes"]))$Freq
  , DROPSEQ_LOW = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_DROPSEQ == "Bulk high"]))$Freq
  , FLUIDIGM_HT_HIGH = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_FLUIDIGM_HT == "Pool high"]))$Freq
  , FLUIDIGM_HT_LOW = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_FLUIDIGM_HT == "Bulk high"]))$Freq
  , FLUIDIGM_LT_HIGH = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_FLUIDIGM_LT == "Pool high"]))$Freq
  , FLUIDIGM_LT_LOW = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_FLUIDIGM_LT == "Bulk high"]))$Freq
  , POLLEN_HIGH = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_POLLEN == "Pool high"]))$Freq
  , POLLEN_LOW = data.frame(table(ggDF$gene_biotype[ggDF$SUBSET_POLLEN == "Bulk high"]))$Freq
)
# Calculate percent of total biotypes
ggDF[2:10] <- apply(ggDF[2:10], 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 100
})
names(ggDF)[1:2] <- c("BIOTYPE", "ALL_GENES")
ggDF <- melt(ggDF)
ggDF$variable <- factor(ggDF$variable, levels = c("DROPSEQ_HIGH"
  , "DROPSEQ_LOW", "FLUIDIGM_HT_HIGH", "FLUIDIGM_HT_LOW", "FLUIDIGM_LT_HIGH"
  , "FLUIDIGM_LT_LOW" ,"POLLEN_HIGH", "POLLEN_LOW", "ALL_GENES"))
colnames(ggDF)[2] <- "DATASET"
# ggplot
ggplot(ggDF, aes(x = BIOTYPE, y = value)) +
  geom_bar(aes(fill = DATASET), position = "dodge", stat = "identity") +
  scale_fill_brewer(type = "qual", palette = 3) +
  coord_flip() +
  ylab("Gene biotypes") +
  xlab("Pooled: log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nFrequency of gene biotypes"
    , "\npooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"
    , "\nPercent of biotype relative to all genes in dataset"))
ggsave(paste0(outGraph, "Biotypes.pdf"), width = 8, height = 8)
################################################################################

### Gene intersections high in bulk or high in pooled

Venn_Diagram_3 <- function (v1, v2, v3, labels, title) {
  print(length(intersect(v1, intersect(v2, v3))))
  venn <- euler(c(
    "A" = length(v1)
    , "B" = length(v2)
    , "C" = length(v3)
    , "A&B" = length(intersect(v1, v2))
    , "A&C" = length(intersect(v1, v3))
    , "B&C" = length(intersect(v2, v3))
    , "A&B&C" = length(intersect(v1, intersect(v2, v3)))))
  plot(venn, counts = TRUE, labels = labels
    , fill = c("#e41a1c", "#377eb8", "#4daf4a")
    , main = paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"
      , "\nHuman fetal brain"))
}
Venn_Diagram_4 <- function (v1, v2, v3, v4, labels, title) {
  print(length(intersect(v4, intersect(v1, intersect(v2, v3)))))
  venn <- euler(c(
    "A" = length(v1)
    , "B" = length(v2)
    , "C" = length(v3)
    , "D" = length(v4)
    , "A&B" = length(intersect(v1, v2))
    , "A&C" = length(intersect(v1, v3))
    , "A&D" = length(intersect(v1, v4))
    , "B&C" = length(intersect(v2, v3))
    , "B&D" = length(intersect(v2, v4))
    , "C&D" = length(intersect(v3, v4))
    , "A&B&C" = length(intersect(v1, intersect(v2, v3)))
    , "A&B&D" = length(intersect(v1, intersect(v2, v4)))
    , "B&C&D" = length(intersect(v2, intersect(v3, v4)))
    , "A&B&C&D" = length(intersect(v1, intersect(v2, intersect(v3, v4))))
    ))
  plot(venn, counts = TRUE, labels = labels
    , fill = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
    , main = paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"
      , "\nHuman fetal brain"))
}

pdf(paste0(outGraph, "VennDiagrams.pdf"))
# GZ
# Bulk high
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high", ]$ENSEMBL_ID
  , c("Drop-seq GZ", "Fluidigm HT GZ", "Pollen")
  , "Genes high in bulk versus pooled"
)
# Pool high
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
  , c("Drop-seq GZ", "Fluidigm HT GZ", "Pollen")
  , "Genes high in pool versus bulk"
)
# GZ CP
# Bulk high
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high", ]$ENSEMBL_ID
  , c("Drop-seq", "Fluidigm HT", "Fluidigm LT")
  , "Genes high in bulk versus pooled"
)
# Pool high
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high", ]$ENSEMBL_ID
  , c("Drop-seq", "Fluidigm HT", "Fluidigm LT")
  , "Genes high in pool versus bulk"
)
# Drop-seq GZ CP + Pollen
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high", ]$ENSEMBL_ID
  , c("Drop-seq", "Drop-seq GZ", "Pollen")
  , "Genes high in Bulk versus pool"
)
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high", ]$ENSEMBL_ID
  , c("Drop-seq", "Drop-seq GZ", "Pollen")
  , "Genes high in Bulk versus pool"
)
# Fluidigm HT GZ CP + Pollen
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
  , c("Fluidigm HT", "Fluidigm HT GZ", "Pollen")
  , "Genes high in pool versus bulk"
)
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
  , c("Fluidigm HT", "Fluidigm HT GZ", "Pollen")
  , "Genes high in pool versus bulk"
)
# Drop-seq + Fluidigm LT GZ CP + Pollen
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
  , c("Drop-seq", "Fluidigm LT", "Pollen")
  , "Genes high in pool versus bulk"
)
Venn_Diagram_3(
  ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high", ]$ENSEMBL_ID
  , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
  , c("Drop-seq", "Fluidigm LT", "Pollen")
  , "Genes high in pool versus bulk"
)
# # Drop-seq GZ, Fluidigm HT GZ, Fluidigm LT, Pollen
#   # Bulk high
# Venn_Diagram_4(
#   ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high", ]$ENSEMBL_ID
#   , c("Drop-seq GZ", "Fluidigm HT GZ", "Fluidigm LT", "Pollen")
#   , "Genes high in bulk versus pooled"
# )
#   # Pool high
# Venn_Diagram_4(
#   ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high", ]$ENSEMBL_ID
#   , ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high", ]$ENSEMBL_ID
#   , c("Drop-seq GZ", "Fluidigm HT GZ", "Fluidigm LT", "Pollen")
#   , "Genes high in pool versus bulk"
# )
dev.off()

intersect(ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high", ]$ENSEMBL_ID
, ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high", ]$ENSEMBL_ID)
################################################################################

### GC, length, MT genes, brain expressed, and biotype of gene intersections
### high in bulk or pooled

# Intersect Dropseq GZ, Fluidigm HT GZ, Pollen
ssl2MnCpmDF$SINGLE_CELL_BIAS_GZ <- "All other genes"
ssl2MnCpmDF$SINGLE_CELL_BIAS_GZ[
  ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high" &
    ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high"] <- "Bulk high"
ssl2MnCpmDF$SINGLE_CELL_BIAS_GZ[
  ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high" &
    ssl2MnCpmDF$SUBSET_POLLEN == "Pool high"] <- "Pool high"

# Intersect Dropseq, Fluidigm HT, Fluidigm LT
ssl2MnCpmDF$SINGLE_CELL_BIAS <- "All other genes"
ssl2MnCpmDF$SINGLE_CELL_BIAS[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Bulk high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high"] <- "Bulk high"
ssl2MnCpmDF$SINGLE_CELL_BIAS[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high" &
    ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high"] <- "Pool high"

# # Dropseq GZ, Fluidigm HT GZ, Fluidigm LT, Exclude Pollen
# ssl2MnCpmDF$LAB_BIAS_GESCHWIND <- "All other genes"
# ssl2MnCpmDF$LAB_BIAS_GESCHWIND[
#   ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high" &
#     ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high" &
#     ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high" &
#     ! ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high"] <- "Bulk high"
# ssl2MnCpmDF$LAB_BIAS_GESCHWIND[
#   ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high" &
#     ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high" &
#     ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high" &
#     ! ssl2MnCpmDF$SUBSET_POLLEN == "Pool high"] <- "Pool high"
# 
# # Include Pollen, Exclude Dropseq GZ, Fluidigm HT GZ, Fluidigm LT
# ssl2MnCpmDF$LAB_BIAS_POLLEN <- "All other genes"
# ssl2MnCpmDF$LAB_BIAS_POLLEN[
#   ! ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high" &
#     ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high"] <- "Bulk high"
# ssl2MnCpmDF$LAB_BIAS_POLLEN[
#   ! ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high" &
#     ssl2MnCpmDF$SUBSET_POLLEN == "Pool high"] <- "Pool high"
# 
# # Include Drop-seq, Exclude Fluidigm HT GZ, Fluidigm LT
# ssl2MnCpmDF$DROPSEQ_BIAS <- "All other genes"
# ssl2MnCpmDF$DROPSEQ_BIAS[
#   ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Bulk high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high"] <- "Bulk high"
# ssl2MnCpmDF$DROPSEQ_BIAS[
#   ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high" &
#     ! ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high"] <- "Pool high"

# Format for ggplot
ggDF <- data.frame(ssl2MnCpmDF, mnExDF$percentage_gc_content, mnExDF$cds_length)
ggDF <- ggDF[! ggDF$mnExDF.cds_length == 0, ]

## Length
ggL <- list(
  Plot_CDS_Length("SINGLE_CELL_BIAS"
  , "Single-cell bias\nIntersect Drop-seq, Fluidigm HT, Fluidigm LT")
  , Plot_CDS_Length("SINGLE_CELL_BIAS_GZ"
    , "Single-cell bias\nIntersect Drop-seq GZ, Fluidigm HT GZ, Pollen")
  , Plot_CDS_Length("SUBSET_DROPSEQ", "Drop-seq")
  , Plot_CDS_Length("SUBSET_DROPSEQ_0", "Drop-seq not detected")
  , Plot_CDS_Length("SUBSET_DROPSEQ_GZ", "Drop-seq GZ")
  , Plot_CDS_Length("SUBSET_FLUIDIGM_HT", "Fluidigm HT")
  , Plot_CDS_Length("SUBSET_FLUIDIGM_HT_GZ", "Fluidigm HT GZ")
  , Plot_CDS_Length("SUBSET_FLUIDIGM_LT", "Fluidigm LT")
  , Plot_CDS_Length("SUBSET_POLLEN", "Pollen")
  )
pdf(paste0(outGraph, "Length.pdf"), width = 8, height = 12)
do.call("grid.arrange", c(ggL, ncol = 2))
dev.off()

# CDS length from biomart, using longest CDS
Plot_CDS_Length_Histogram <- function(exM, title, geneIdType){
  cGene <- rowSums(exM)
  ggDF <- data.frame(GENE = names(cGene), COUNTS = cGene)
  ggDF <- merge(ggDF, mnExDF[c(geneIdType, "cds_length")]
    , by.x = "GENE", by.y = geneIdType)
  ggDF <- data.frame(CDS_LENGTH = rep(ggDF$cds_length, ggDF$COUNTS))
  gg <- ggplot(ggDF, aes(x = CDS_LENGTH)) +
    geom_histogram(binwidth = 500) +
    # coord_cartesian(xlim = c(0, 10000)) +
    scale_x_continuous(limits = c(0, 10000)) +
    theme(text = element_text(size = 12)) +
    xlab("CDS length") +
    ylab("Counts") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\nHistogram of raw counts versus CDS length"
      , "\n", title
      , "\n"))
  return(gg)  
}

ggL <- list(
  Plot_CDS_Length_Histogram(blExDF, "Bulk RNA-seq", "ensembl_gene_id")
  # , Plot_CDS_Length_Histogram(dsExM, "Drop-seq", "hgnc_symbol")
  # , Plot_CDS_Length_Histogram(gzDsExM, "Drop-seq GZ", "hgnc_symbol")
  , Plot_CDS_Length_Histogram(dsExM, "Drop-seq", "ensembl_gene_id")
  , Plot_CDS_Length_Histogram(gzDsExM, "Drop-seq GZ", "ensembl_gene_id")
  , Plot_CDS_Length_Histogram(fhExDF, "Fluidigm HT", "ensembl_gene_id")
  , Plot_CDS_Length_Histogram(gzFhExDF, "Fluidigm HT GZ", "ensembl_gene_id")
  , Plot_CDS_Length_Histogram(flExDF, "Fluidigm LT", "ensembl_gene_id")
  , Plot_CDS_Length_Histogram(pnExDF, "Pollen", "ensembl_gene_id")
)
# pdf(paste0(outGraph, "LengthCounts_Histogram.pdf"), width = 8, height = 12)
# ggL[[1]]
# dev.off()

png(paste0(outGraph, "LengthCounts_Histogram.png"), width = 8, height = 20
  , units = "in", res = 300)
do.call("grid.arrange", c(ggL, ncol = 2))
dev.off()


# CDS length from biomart, using longest CDS
Plot_CDS_Length <- function(geneSet, title) {
  gg <- ggplot(ggDF, aes_string(x = geneSet, y = "mnExDF.cds_length")) +
    geom_boxplot(fill = c("#e41a1c", "#377eb8", "#4daf4a")) +
    coord_cartesian(ylim = c(0, 5000)) +
    theme(text = element_text(size = 12)) +
    xlab("Gene subset") +
    ylab("CDS length") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"))
  return(gg)
}

## GC
ggL <- list(
  Plot_GC("SINGLE_CELL_BIAS"
    , "Single-cell bias\nIntersect Drop-seq, Fluidigm HT, Fluidigm LT")
  , Plot_GC("SINGLE_CELL_BIAS_GZ"
    , "Single-cell bias GZ\nIntersect Drop-seq GZ, Fluidigm HT GZ, Pollen")
  , Plot_GC("SUBSET_DROPSEQ", "Drop-seq")
  , Plot_GC("SUBSET_DROPSEQ_0", "Drop-seq not detected")
  , Plot_GC("SUBSET_DROPSEQ_GZ", "Drop-seq GZ")
  , Plot_GC("SUBSET_FLUIDIGM_HT", "Fluidigm HT")
  , Plot_GC("SUBSET_FLUIDIGM_HT_GZ", "Fluidigm HT GZ")
  , Plot_GC("SUBSET_FLUIDIGM_LT", "Fluidigm LT")
  , Plot_GC("SUBSET_POLLEN", "Pollen")
)
pdf(paste0(outGraph, "GC.pdf"), width = 8, height = 12)
do.call("grid.arrange", c(ggL, ncol = 2))
dev.off()

## Biotype
# Format for ggplot2 and calculate percent of each biotype in each subset
ggDF <- data.frame(ssl2MnCpmDF, mnExDF$gene_biotype)
ggDF <- do.call("rbind", lapply(c("SINGLE_CELL_BIAS", "SINGLE_CELL_BIAS_GZ"
  , "SUBSET_DROPSEQ", "SUBSET_DROPSEQ_0", "SUBSET_DROPSEQ_GZ"
  , "SUBSET_FLUIDIGM_HT", "SUBSET_FLUIDIGM_HT_GZ", "SUBSET_FLUIDIGM_LT"
  , "SUBSET_POLLEN"), function (dataset){
  ldf <- split(ggDF, ggDF[ ,c(dataset)])
  ggDF <- do.call("rbind", (lapply(names(ldf), function (name) {
    print(name)
    df <- ldf[[name]]
    df <- data.frame(table(df$mnExDF.gene_biotype))
    df$PERCENT <- (df$Freq / sum(df$Freq)) * 100
    df$SUBSET <- name
    return(df)
  })))
  ggDF$DATASET <- dataset
  return(ggDF)
}))
# ggplot
ggplot(ggDF, aes(x = Var1, y = PERCENT)) +
  facet_wrap("DATASET") +
  geom_bar(aes(fill = SUBSET), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  coord_flip() +
  theme(text = element_text(size = 12)) +
  xlab("Gene biotypes") +
  ylab("Percent of genes") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nFrequency of gene biotypes"
    , "\npooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"
    , "\nPercent of biotype relative to all genes in subset"))
ggsave(paste0(outGraph, "Biotypes.pdf"), width = 9.5, height = 14)

## Percent MT

ggL <- list(
  Plot_MT_Percent("SINGLE_CELL_BIAS"
    , "Single-cell bias\nIntersect Drop-seq, Fluidigm HT, Fluidigm LT")
  , Plot_MT_Percent("SINGLE_CELL_BIAS_GZ"
    , "Single-cell bias GZ\nIntersect Drop-seq GZ, Fluidigm HT GZ, Pollen")
  , Plot_MT_Percent("SUBSET_DROPSEQ", "Drop-seq")
  , Plot_MT_Percent("SUBSET_DROPSEQ_0", "Drop-seq not detected")
  , Plot_MT_Percent("SUBSET_DROPSEQ_GZ", "Drop-seq GZ")
  , Plot_MT_Percent("SUBSET_FLUIDIGM_HT", "Fluidigm HT")
  , Plot_MT_Percent("SUBSET_FLUIDIGM_HT_GZ", "Fluidigm HT GZ")
  , Plot_MT_Percent("SUBSET_FLUIDIGM_LT", "Fluidigm LT")
  , Plot_MT_Percent("SUBSET_POLLEN", "Pollen")
)
pdf(paste0(outGraph, "MT_Percent.pdf"), width = 9, height = 12)
do.call("grid.arrange", c(ggL, ncol = 2))
dev.off()

## Overlap with brain expressed genes

# Percentage of brain expressed genes
Plot_Brain_Expressed_Percent <- function(geneSet, title) {
  # Percent
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "All other genes", ]
  genes <- df$ENSEMBL_ID
  pct <- sum(brexpDF[ ,1] %in% genes) / length(genes)
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Bulk high", ]
  genes <- df$ENSEMBL_ID
  pct <- c(pct, sum(brexpDF[ ,1] %in% genes) / length(genes))
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Pool high", ]
  genes <- df$ENSEMBL_ID
  pct <- c(pct, sum(brexpDF[ ,1] %in% genes) / length(genes))
  pct <- pct*100
  # Format for ggplot
  ggDF <- data.frame(SUBSET = c("All other genes", "Bulk high", "Pool high")
    , PERCENT_MT = pct)
  # ggplot
  gg <- ggplot(ggDF, aes(x = SUBSET, y = PERCENT_MT)) +
    geom_bar(stat = "identity", fill = c("#e41a1c", "#377eb8", "#4daf4a")) +
    coord_cartesian(ylim = c(0, 100)) +
    theme(text = element_text(size = 12)) +
    xlab("Gene subset") +
    ylab("Percent brain expressed") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", title
      , "\n"))
  return(gg)
}

ggL <- list(
  Plot_Brain_Expressed_Percent("SINGLE_CELL_BIAS"
    , "Single-cell bias\nIntersect Drop-seq, Fluidigm HT, Fluidigm LT")
  , Plot_Brain_Expressed_Percent("SINGLE_CELL_BIAS_GZ"
    , "Single-cell bias GZ\nIntersect Drop-seq GZ, Fluidigm HT GZ, Pollen")
  , Plot_Brain_Expressed_Percent("SUBSET_DROPSEQ", "Drop-seq")
  , Plot_Brain_Expressed_Percent("SUBSET_DROPSEQ_0", "Drop-seq not detected")
  , Plot_Brain_Expressed_Percent("SUBSET_DROPSEQ_GZ", "Drop-seq GZ")
  , Plot_Brain_Expressed_Percent("SUBSET_FLUIDIGM_HT", "Fluidigm HT")
  , Plot_Brain_Expressed_Percent("SUBSET_FLUIDIGM_HT_GZ", "Fluidigm HT GZ")
  , Plot_Brain_Expressed_Percent("SUBSET_FLUIDIGM_LT", "Fluidigm LT")
  , Plot_Brain_Expressed_Percent("SUBSET_POLLEN", "Pollen")
)
pdf(paste0(outGraph, "BrainExpressed_Percent.pdf"), width = 9, height = 12)
do.call("grid.arrange", c(ggL, ncol = 2))
dev.off()

## HGNC in Drop-seq 0 high in bulk

df <- data.frame(table(ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_0 == "Bulk high", 2] == 0)
  , table(ssl2MnCpmDF[ssl2MnCpmDF$SUBSET_DROPSEQ_0 == "Pool high", 2] == 0))
df <- df[c(2,4)]
# All pool high have hgnc symbols so filling in 0 for no hgnc symbol in table
df[2,2] <- 0
colnames(df) <- c("BULK_HIGH", "POOL_HIGH")
row.names(df) <- c("HGNC_SYMBOL", "NO_HGNC_SYMBOL")
write.csv(df, paste0(outTable, "Dropseq0_Number_HGNCsymbol.csv"))
################################################################################

### Write table of subset gene lists

# List of dataframes, each dataframe is a subset
cols <- grep("SUBSET", names(ssl2MnCpmDF), value = TRUE)
names(cols) <- cols
l <- lapply(cols, function(dfcol) {
  l <- list(
    as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Bulk high", 1])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Bulk high", 2])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Pool high", 1])
    , as.character(ssl2MnCpmDF[ssl2MnCpmDF[ ,dfcol] == "Pool high", 2])
    )
  # Data frame of lists of different lengths
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  mat <- sapply(l, "[", i = seq.max)
  print(str(mat))
  colnames(mat) <- c("BULK_HIGH_ENSEMBL", "BULK_HIGH_HGNC", "POOL_HIGH_ENSEMBL", "POOL_HIGH_HGNC") 
  # Convert dataframe to lists
  l <- as.list(as.data.frame(mat, stringsAsFactors = FALSE))
  return(l)
})
# Format column names
ll <- lapply(names(l), function(datasetName) {
  dataset <- l[[datasetName]]
  colNames <- paste0(datasetName, "_", names(dataset))
  names(dataset) <- colNames
  print(str(dataset))
  return(dataset)
})
# Append NAs to dataframes to make same length
n.obs <- sapply(ll, function(l) sapply(l, length))
seq.max <- seq_len(max(n.obs[1, ]))
lm <- lapply(ll, function(l) sapply(l, "[", i = seq.max))
m <- do.call("cbind", lm)
m[m == "0"] <- NA

# Write out
write.csv(m, file = paste0(outTable, "Subset_Gene_Lists.csv"))
################################################################################

### GO and pathway analysis

## gprofiler
gproLDF <- list(
  # High Drop-seq, Fluidigm HT, Fluidigm LT
  Single_Cell_Bias_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SINGLE_CELL_BIAS == "Pool high"]))
  # Low 
  , Single_Cell_Bias_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SINGLE_CELL_BIAS == "Bulk high"]))
  # High Drop-seq GZ, Fluidigm HT GZ, Pollen
  , Single_Cell_Bias_GZ_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SINGLE_CELL_BIAS_GZ == "Pool high"]))
  # Low 
  , Single_Cell_Bias_GZ_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SINGLE_CELL_BIAS_GZ == "Bulk high"]))
  # High Drop-seq
  , Dropseq_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high"]))
  # Low Drop-seq
  , Dropseq_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high"]))
  # High Drop-seq 0
  , Dropseq_0_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ_0 == "Pool high"]))
  # Low Drop-seq 0
  , Dropseq_0_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ_0 == "Bulk high"]))
  # High Fluidigm HT
  , Fluidigm_HT_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Pool high"]))
  # Low Fluidigm HT
  , Fluidigm_HT_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT == "Bulk high"]))
  # High Fluidigm LT
  , Fluidigm_LT_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Pool high"]))
  # Low Fluidigm LT
  , Fluidigm_HT_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_LT == "Bulk high"]))
  # High Drop-seq GZ
  , Dropseq_GZ_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Pool high"]))
  # Low Drop-seq GZ
  , Dropseq_GZ_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_DROPSEQ_GZ == "Bulk high"]))
  # High Fluidigm HT GZ
  , Fluidigm_HT_GZ_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Pool high"]))
  # Low Fluidigm HT GZ
  , Fluidigm_HT_GZ_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_FLUIDIGM_HT_GZ == "Bulk high"]))
  # High Pollen
  , Pollen_High_Pool = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_POLLEN == "Pool high"]))
  # Low Pollen
  , Pollen_HT_High_Bulk = gprofiler(query = as.character(ssl2MnCpmDF$HGNC_SYMBOL[ssl2MnCpmDF$SUBSET_POLLEN == "Bulk high"]))
)

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))

# Write out for REVIGO
# Drop-seq, Fluidigm HT, Fluidigm LT
# High bulk
lapply(names(gproLDF), function(dataset) {
  # Full gprofiler table
  goproDF <- gproLDF[[dataset]]
  write.csv(goproDF, paste0(outTableGpro, "GO_", dataset, ".csv")
    , quote = FALSE, row.names = FALSE)
  # Term IDs for inputing into REVIGO
  go <- gproLDF[[dataset]]$term.id[
    gproLDF[[dataset]]$domain == "BP" |
      gproLDF[[dataset]]$domain == "MF" |
      gproLDF[[dataset]]$domain == "CC"]
  write.table(go, paste0(outTableGproIds, "GO_TermID_", dataset, ".txt")
    , quote = FALSE, row.names = FALSE)
})

# ## Intersect GO terms
# 
# # Drop-seq GZ, Fluidigm HT GZ, Pollen
# # High pool
# go <- intersect(
#   intersect(gpPhiBloDsGzDF$term.id[
#       gpPhiBloDsGzDF$domain == "BP" |
#       gpPhiBloDsGzDF$domain == "MF" |
#       gpPhiBloDsGzDF$domain == "CC"]
#     , gpPhiBloFhGzDF$term.id)
#   , gpPhiBloPnDF$term.id)
# write.table(go, paste0(outTable, "GO_DsGz_FhGz_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(gpPloBhiDsGzDF$term.id[
#       gpPloBhiDsGzDF$domain == "BP" |
#       gpPloBhiDsGzDF$domain == "MF" |
#       gpPloBhiDsGzDF$domain == "CC"]
#     , gpPloBhiFhGzDF$term.id)
#   , gpPloBhiPnDF$term.id)
# write.table(go, paste0(outTable, "GO_DsGz_FhGz_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# 
# # Drop-seq, Fluidigm HT, Fluidigm LT
# # High pool
# go <- intersect(
#   intersect(gpPhiBloDsDF$term.id[
#       gpPhiBloDsDF$domain == "BP" |
#       gpPhiBloDsDF$domain == "MF" |
#       gpPhiBloDsDF$domain == "CC"]
#     , gpPhiBloFhDF$term.id)
#   , gpPhiBloFlDF$term.id)
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(gpPloBhiDsDF$term.id[
#       gpPloBhiDsDF$domain == "BP" |
#       gpPloBhiDsDF$domain == "MF" |
#       gpPloBhiDsDF$domain == "CC"]
#     , gpPloBhiFhDF$term.id)
#   , gpPloBhiFlDF$term.id)
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# 
# # Drop-seq, Fluidigm HT, Fluidigm LT, Pollen
# # High pool
# go <- intersect(
#   intersect(
#     intersect(gpPhiBloDsDF$term.id[
#         gpPhiBloDsDF$domain == "BP" |
#         gpPhiBloDsDF$domain == "MF" |
#         gpPhiBloDsDF$domain == "CC"]
#       , gpPhiBloFhDF$term.id)
#     , gpPhiBloFlDF$term.id)
#   , gpPhiBloPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- intersect(
#   intersect(
#     intersect(gpPloBhiDsDF$term.id[
#         gpPloBhiDsDF$domain == "BP" |
#         gpPloBhiDsDF$domain == "MF" |
#         gpPloBhiDsDF$domain == "CC"]
#       , gpPloBhiFhDF$term.id)
#     , gpPloBhiFlDF$term.id)
#   , gpPloBhiPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_Fh_Fl_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# 
# # Drop-seq, not in Fluidigm HT, Fluidigm LT
# # High pool
# go <- setdiff(
#   setdiff(gpPhiBloDsDF$term.id[
#       gpPhiBloDsDF$domain == "BP" |
#       gpPhiBloDsDF$domain == "MF" |
#       gpPhiBloDsDF$domain == "CC"]
#     , gpPhiBloFhDF$term.id)
#   , gpPhiBloPnDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_NotInFhFl_Pollen_PoolHigh.txt")
#   , quote = FALSE, row.names = FALSE)
# # High bulk
# go <- setdiff(
#   setdiff(gpPloBhiDsDF$term.id[
#       gpPloBhiDsDF$domain == "BP" |
#       gpPloBhiDsDF$domain == "MF" |
#       gpPloBhiDsDF$domain == "CC"]
#     , gpPloBhiFhDF$term.id)
#   , gpPloBhiFlDF$term.id
# )
# write.table(go, paste0(outTable, "GO_Ds_NotInFhFl_Pollen_BulkHigh.txt")
#   , quote = FALSE, row.names = FALSE)

## Gene length for each GO ID

dir.create("../analysis/graphs/Pooled_Vs_Bulk_GO_CDS_Length", recursive = TRUE)
# Loop through gprofiler results for each subset
lapply(names(gproLDF), function(name) {
  print(name)
  gproDF <- gproLDF[[name]]
  if (nrow(gproDF) > 0) {
    # Only BP GO terms
    gproDF <- gproDF[gproDF$domain == "BP", ]
    # Genes
    l <- strsplit(gproDF$intersection, ",")
    names(l) <- gproDF$term.name
    # Data frame of genes and term names
    ggDF <- do.call("rbind", lapply(names(l), function(name) {
      data.frame(hgnc_symbol = l[[name]], term.name = name)  
    }))
    # Add CDS length
    ggDF <- merge(ggDF, mnExDF[c("hgnc_symbol", "cds_length")], by = "hgnc_symbol")
    # Take longest CDS (some hgnc_symbol occur twice in from matching multiple
    # ensembl IDs)
    ggDF <- aggregate(cds_length ~ hgnc_symbol + term.name, data = ggDF, max)
    # Add all genes
    df <- mnExDF[mnExDF$cds_length > 0, c("hgnc_symbol", "cds_length")]
    df[ ,3] <- df$cds_length
    df[ ,2] <- "All genes"
    names(df)[2:3] <- c("term.name", "cds_length")
    ggDF <- rbind(ggDF, df)
    # Term to color all genes box only
    ggDF$AREA_COLOR <- "Genes GO intersection"
    ggDF$AREA_COLOR[ggDF$term.name == "All genes"] <- "All genes"
    # ggplot
    ggplot(ggDF, aes(y = cds_length, x = term.name, fill = AREA_COLOR)) +
      geom_boxplot() +
      scale_fill_brewer(type = "qual", palette = 6) +
      coord_flip(ylim = c(0, 6000)) +
      theme(text = element_text(size = 12))
    ggsave(paste0(
      "../analysis/graphs/Pooled_Vs_Bulk_GO_CDS_Length/Pooled_Vs_Bulk_GO_CDS_Length_"
      , name, ".pdf"), height = (length(unique(ggDF$term.name))*.1 + 4)
      , width = 12)
  }
  else {print(paste0("No significant GO terms for ", name))}
})

# # "BP"  "CC"  "MF"  "cor" "hp"  "hpa" "keg" "mi"  "omi" "rea" "tf" 
# 
# pdf(paste0(outGraph, "gProfileR.pdf"), width = 9, height = 7.5)
# # Title page
# plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
# text(5, 8, "GO term enrichment using gProfileR of genes")
# text(5, 7, "higher in pooled versus bulk or higher in bulk versus pooled")
# text(5, 6, "Multiple hypothesis correction using g:SCS")
# text(5, 5, "g:SCS stringency is in between FDR and Bonferroni and was designed for GO term enrichment analysis")
# text(5, 4, "All human genes used as background")
# # Drop-seq
# Graph_GO_Enrichment(gpPloBhiDsDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Drop-seq"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloDsDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Drop-seq"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Fluidigm HT
# Graph_GO_Enrichment(gpPloBhiFhDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm HT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFhDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidigm HT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Fluidigm LT
# Graph_GO_Enrichment(gpPloBhiFlDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm LT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFlDF, "BP", 3
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidgm LT"
#     , "\nHuman Fetal Brain CP GZ"
#     , "\ngProfiler")
# )
# # Drop-seq GZ
# Graph_GO_Enrichment(gpPloBhiDsGzDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Drop-seq"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloDsGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Drop-seq"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# # Fluidigm HT GZ
# Graph_GO_Enrichment(gpPloBhiFhGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Fluidigm HT"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloFhGzDF, "BP", 5
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Fluidigm HT"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# # Pollen
# Graph_GO_Enrichment(gpPloBhiPnDF, "BP", 4
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nLow in Pollen"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# Graph_GO_Enrichment(gpPhiBloPnDF, "BP", 2
#   , paste0(graphCodeTitle
#     , "\n"
#     , "\nGO Enrichment: Pooled scRNAseq vs Bulk RNAseq"
#     , "\n"
#     , "\nHigh in Pollen"
#     , "\nHuman Fetal Brain GZ"
#     , "\ngProfiler")
# )
# dev.off()

# # GO Terms molecular function
# ggDF <- subset(gpPlBhDF, domain == "MF" & relative.depth == 5)
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
# 
# # Kegg
# ggDF <- subset(gpPlBhDF, domain == "keg")
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
# 
# # Omi
# ggDF <- subset(gpPlBhDF, domain == "omi")
# ggDF <- ggDF[order(ggDF$p.value), ]
# # Top 30 by p-value
# ggDF <- ggDF[1:30, ]
# ggDF$NEG_LOG10_PVALUE <- -log(ggDF$p.value, 10)
# # Reverse order for ggplot so most significant at top of graph
# ggDF <- ggDF[order(ggDF$NEG_LOG10_PVALUE), ]
# ggDF$term.name <- factor(ggDF$term.name, levels = ggDF$term.name)
# # ggplot
# ggplot(ggDF, aes(y = NEG_LOG10_PVALUE, x = term.name)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
################################################################################

### Expression of genes biased towards pool or bulk by cluster

Heatmap_Expression_Clusters_Format <- function (geneListDF, expLimitHigh
  , expLimitLow){
  # Subset to genes, merge to keep duplicate genes
  # e.g. if plotting cluster markers, may be the same marker of multiple clusters
  ggDF <- merge(geneListDF, centSO@scale.data, by.x = 1, by.y = "row.names"
    , all.x = TRUE)
  row.names(ggDF) <- paste0(length(ggDF[ ,1]):1, "_", ggDF[ ,1])
  ggDF <- ggDF[ ,-c(1)]
  # Order by clustering
  idx <- match(colnames(ggDF), names(sort(centSO@ident)))
  ggDF <- ggDF[ ,idx]
  # Format for ggplot2
  ggDF <- as.matrix(ggDF)
  ggDF <- melt(ggDF)
  # Add clusters
  idx <- match(ggDF$Var2, names(centSO@ident))
  ggDF$CLUSTERS <- centSO@ident[idx]
  # Set sample order by clustering
  ggDF$Var2 <- factor(ggDF$Var2, levels = names(sort(centSO@ident)))
  # Set expression limits
  ggDF$value[ggDF$value > expLimitHigh] <- expLimitHigh
  ggDF$value[ggDF$value < expLimitLow] <- expLimitLow
  # Return data frame formatted for ggplot
  return(ggDF)
}

Heatmap_Expression_Clusters <- function (ggDF, graphTitle) {
  # ggplot
  ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    facet_grid(~CLUSTERS, space = "free", scales = "free") +
    # scale_fill_gradient2(high = "#d7191c", low = "#2c7bb6")
    scale_y_discrete(labels = gsub(".*_", "", ggDF$Var1)) + 
    scale_fill_distiller(name = "Normalized\nexpression", type = "div"
      , palette = 5, direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(text = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    ylab("Genes") +
    xlab("Cells") +
    ggtitle(paste0(graphCodeTitle
      , "\n"
      , "\n", graphTitle
      , "\nMean centered and scaled normalized expression"
      , "\n"))
}

geneListDF <- data.frame(GENE = ssl2MnCpmDF$HGNC_SYMBOL[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Bulk high"])
geneListDF <- data.frame(GENE = geneListDF[! geneListDF$GENE == 0, ])
ggDF <- Heatmap_Expression_Clusters_Format(geneListDF, 5, -1)
gg <- Heatmap_Expression_Clusters(
  ggDF, "Genes higher in bulk versus pooled Drop-seq")
ggsave(paste0(outGraph, "Heatmap_Dropseq_BulkHigh.png")
  , width = 10, height = 10, dpi = 600)

geneListDF <- data.frame(GENE = ssl2MnCpmDF$HGNC_SYMBOL[
  ssl2MnCpmDF$SUBSET_DROPSEQ == "Pool high"])
geneListDF <- data.frame(GENE = geneListDF[! geneListDF$GENE == 0, ])
ggDF <- Heatmap_Expression_Clusters_Format(geneListDF, 5, -1)
gg <- Heatmap_Expression_Clusters(
  ggDF, "Genes higher in pooled Drop-seq versus bulk")
ggsave(paste0(outGraph, "Heatmap_Dropseq_PoolHigh.png")
  , width = 10, height = 10, dpi = 600)
################################################################################

### Known cell type markers from Luis

## Intersection with pooled vs bulk subsets
Intersect_KnownMarks_Subsets <- function(sub, bias) {
  genes <- ssl2MnCpmDF[ssl2MnCpmDF[ ,c(sub)] == bias, ]$HGNC_SYMBOL
  kmDF$Gene.Symbol %in% genes
}
l <- list(Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Pool high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Pool high")
)
nr <- max(sapply(l, length))
nc <- length(ldf)
m <- matrix(NA, nr, nc)
for(i in 1:length(l)){
  m[ ,i] <- l[[i]]
}
colnames(m) <- c(
  "SUBSET_DROPSEQ_BULK_HIGH"
  , "SUBSET_DROPSEQ_POOL_HIGH"
  , "SUBSET_DROPSEQ_0_BULK_HIGH"
  , "SUBSET_DROPSEQ_0_POOL_HIGH"
  , "SINGLE_CELL_BIAS_GZ_BULK_HIGH"
  , "SINGLE_CELL_BIAS_GZ_POOL_HIGH"
  , "SINGLE_CELL_BIAS_BULK_HIGH"
  , "SINGLE_CELL_BIAS_POOL_HIGH"
)
rnames <- paste0(kmDF$Group, " ", kmDF$Gene.Symbol)
row.names(m) <- rnames
ggDF <- melt(m)
ggDF$Group <- gsub(" .*", "", ggDF$Var1)
ggDF$Var1 <- gsub(".* ", "", ggDF$Var1)
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_wrap("Group", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  ylab("Variant location") +
  xlab("Cell ID") +
  ggtitle(paste0(graphCodeTitle
    , "\nVariants called for each Cell: Fluidigm HT, Human VZ and CP"
    , "\n"))
ggsave(paste0(outGraph, "KnownMarkers_In_Subsets_Heatmap.pdf"), height = 24, width = 16)


## Intersection with pooled vs bulk subsets
Intersect_KnownMarks_Subsets <- function(sub, bias) {
  genes <- ssl2MnCpmDF[ssl2MnCpmDF[ ,c(sub)] == bias, ]$HGNC_SYMBOL
  as.list(kmDF[kmDF$Gene.Symbol %in% genes, 1:3])
}
ll <- list(Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ", "Pool high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Bulk high")
  , Intersect_KnownMarks_Subsets("SUBSET_DROPSEQ_0", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS_GZ", "Pool high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Bulk high")
  , Intersect_KnownMarks_Subsets("SINGLE_CELL_BIAS", "Pool high")
  )
l <- unlist(ll, recursive = FALSE)
nr <- max(sapply(l, length))
nc <- 3*length(ldf)
m <- matrix(NA, nr, nc)
for(i in 1:length(l)){
  m[ ,i] <- as.character(l[[i]][1:nr])  
}
m[is.na(m)] <- ""
colnames(m) <- c(
  rep("SUBSET_DROPSEQ_BULK_HIGH", 3)
  , rep("SUBSET_DROPSEQ_POOL_HIGH", 3)
  , rep("SUBSET_DROPSEQ_0_BULK_HIGH", 3)
  , rep("SUBSET_DROPSEQ_0_POOL_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_GZ_BULK_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_GZ_POOL_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_BULK_HIGH", 3)
  , rep("SINGLE_CELL_BIAS_POOL_HIGH", 3)
)
write.csv(m, paste0(outTable, "KnownMarkers_In_Subsets.csv"), quote = FALSE
  , row.names = FALSE)

## Expression pooled vs bulk

# Z-score
sdDs <- sd(l2MnCpmDF$RATIO_DROPSEQ)
mnDs <- mean(l2MnCpmDF$RATIO_DROPSEQ)
# Subset to marker genes of interest for Luis' excel file
df <- merge(l2MnCpmDF, kmDF, by.x = "HGNC_SYMBOL", by.y = "Gene.Symbol")
# df <- l2MnCpmDF[l2MnCpmDF$HGNC_SYMBOL %in% kmDF$Gene.Symbol, ]
df$Z_SCORE <- (df$RATIO_DROPSEQ - mnDs) / sdDs
split(df, df$Grouping)

# Markers individually bar graph
ggplot(df, aes(x = HGNC_SYMBOL, y = Z_SCORE)) +
  facet_wrap("Grouping", scales = "free_x") +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(text = element_text(size = 12)) +
  ylab("Z-score (<- Bulk, Pooled ->)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nEnrichment or depletion of cell type markers in pooled versus bulk"
    , "\n"
    , "\nDrop-seq GZ CP"
    , "\nRatio of log2(CPM + 1) pooled versus bulk"))
ggsave(paste0(outGraph, "KnownMarkers.pdf"), width = 12, height = 12)

# Markers combined as boxplot
ggplot(df, aes(x = Grouping, y = Z_SCORE)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Z-score (<- Bulk, Pooled ->)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nEnrichment or depletion of cell type markers in pooled versus bulk"
    , "\n"
    , "\nDrop-seq GZ CP"
    , "\nRatio of log2(CPM + 1) pooled versus bulk"))
ggsave(paste0(outGraph, "KnownMarkersCombined.pdf"), width = 8, height = 6)
################################################################################

# ### Cluster markers from Seurat
# 
# # Z-score
# sdDs <- sd(l2MnCpmDF$RATIO_DROPSEQ)
# mnDs <- mean(l2MnCpmDF$RATIO_DROPSEQ)
# # Subset to marker genes of interest for Luis' excel file
# df <- merge(l2MnCpmDF, mkDF, by.x = "HGNC_SYMBOL", by.y = "gene")
# # df <- l2MnCpmDF[l2MnCpmDF$HGNC_SYMBOL %in% kmDF$Gene.Symbol, ]
# df$Z_SCORE <- (df$RATIO_DROPSEQ - mnDs) / sdDs
# df$cluster <- as.factor(df$cluster)
# 
# # Markers individually bar graph
# ggplot(df, aes(x = HGNC_SYMBOL, y = Z_SCORE)) +
#   facet_wrap("Grouping", scales = "free_x") +
#   geom_bar(stat = "identity") +
#   coord_cartesian(ylim = c(-5, 5)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   theme(text = element_text(size = 12)) +
#   ylab("Z-score (<- Bulk, Pooled ->)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nEnrichment or depletion of cluster markers in pooled versus bulk"
#     , "\n"
#     , "\nDrop-seq GZ CP"
#     , "\nRatio of log2(CPM + 1) pooled versus bulk"))
# ggsave(paste0(outGraph, "SeuratMarkers.pdf"), width = 12, height = 12)
# 
# # Markers combined as boxplot
# ggplot(df, aes(x = cluster, y = Z_SCORE)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ylab("Z-score (<- Bulk, Pooled ->)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nEnrichment or depletion of cluster markers in pooled versus bulk"
#     , "\n"
#     , "\nDrop-seq GZ CP"
#     , "\nRatio of log2(CPM + 1) pooled versus bulk"))
# ggsave(paste0(outGraph, "SeuratMarkersCombined.pdf"), width = 8, height = 6)
# ################################################################################