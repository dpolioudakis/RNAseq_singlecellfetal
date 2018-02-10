# Damon Polioudakis
# 2017-04-05
# Pooled versus bulk
# R 3.4.0
################################################################################

rm(list = ls())
sessionInfo()
set.seed(27)

require(Seurat)
require(biomaRt)
require(gProfileR)
require(reshape2)
require(eulerr)
require(gridExtra)
require(cqn)
require(cowplot)
source("Function_Library.R")

# load("../analysis/analyzed_data/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_Workspace.RData")

# Seurat
load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")
ds_ex_DF <- as.data.frame(as.matrix(centSO@raw.data))
centSO@raw.data <- NA
centSO@data <- NA
# load("../analysis/analyzed_data/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_TEST_seuratO.Robj")
# ds_ex_DF <- as.data.frame(as.matrix(ssCentSO@raw.data))
# centSO <- ssCentSO

# Pollen
pollen_ex_DF <- read.csv("../pollen_2015/data/htseq/Exprs_HTSCexon.csv")
# # Picard Sequencing Statistics - bulk RNAseq
# picStatsBuDF <- read.table("../../kriegstein_2015/metadata/PicardToolsQC.csv")

# Fluidigm LT
fldm_LT_ex_DF <- read.csv("../C196-001_002_DP/data/htseq/merged/Exprs_HTSCexon.csv"
  , row.names = 1)
fldm_LT_mt_DF <- read.table("../C196-001_002_DP/metadata/Compiled_Metadata_20160218.txt"
  , fill = TRUE, header = TRUE)

# Fluidigm HT
fldm_HT_ex_DF <- read.csv(
  "../HT-003_DP/data/htseq/GRCh37.75_NoERCC/Exprs_HTSCexon_FtMm3e4.csv")

# Bulk
bulk_ex_DF <- read.csv("../bulk_VZ_CP_from_ATAC/data/htseq/Exprs_HTSCexon.csv"
  , header = TRUE, row.names = 1)
# Metadata
bulk_mt_DF <- read.csv("../bulk_VZ_CP_from_ATAC/metadata/VZCP_sampleinfo.csv")

# biomart
bmDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , header = TRUE)

# Gene lengths and GC content for Union Exon model to normalize bulk by length
load("../source/ENSEMBLhg19_UnionAnno.rda")

# Known cell type markers from Luis
kmDF <- read.csv("../source/MarkersforSingleCell_2017-10-11_Markers.csv"
  , header = TRUE, fill = TRUE)

# # Markers from Seurat
# mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_Top20.txt"
#   , header = TRUE)

# Brain expressed genes from Vivek
brain_exp_DF <- read.csv("../source/BrainExpressed_GTex.csv", header = FALSE)

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
outAnalysis <- "../analysis/analyzed_data/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_"
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

# Subset to higher in bulk or higher in pooled
Subset_Two_SD <- function(dataset) {
  idx <- data.frame(
    High_Bulk = dataset < (mean(dataset) - 2*sd(dataset))
    , High_Pool = dataset > (mean(dataset) + 2*sd(dataset))
  )
  return(idx)
}

Venn_Diagram_3 <- function(v1, v2, v3, labels, title) {
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

Venn_Diagram_3_Wrapper <- function(
  datasets
  , bias_direction
  , title = paste0(paste(datasets, collapse = ", "), ": ", bias_direction)){

  ensembl_gene_id_L <- lapply(datasets, function(dataset){
    print(dataset)
    ensembl_gene_id <- row.names(bias_flag_DF)[
      bias_flag_DF[[dataset]] == bias_direction]
    return(ensembl_gene_id)
  })
  print(head(ensembl_gene_id_L[[1]]))
  Venn_Diagram_3(
    v1 = ensembl_gene_id_L[[1]]
    , v2 = ensembl_gene_id_L[[2]]
    , v3 = ensembl_gene_id_L[[3]]
    , labels = datasets
    , title = title
  )
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

### Cleaning and formatting

## Subset drop-seq raw counts matrix to cells that passed QC

ds_ex_DF <- ds_ex_DF[ ,colnames(ds_ex_DF) %in% row.names(centSO@meta.data)]

## Remove ERCCs and STAR stats from tail and move gene names to row.names

# Remove ERCCs and STAR stats from tail
# Bulk
tail(bulk_ex_DF, 10)
bulk_ex_DF <- head(bulk_ex_DF, -97)
tail(bulk_ex_DF, 5)
# Pollen
tail(pollen_ex_DF, 10)
pollen_ex_DF <- head(pollen_ex_DF, -97)
tail(pollen_ex_DF, 5)
# Fluidigm LT
tail(fldm_LT_ex_DF, 10)
fldm_LT_ex_DF <- head(fldm_LT_ex_DF, -97)
tail(fldm_LT_ex_DF, 5)
# Fluidigm HT
tail(fldm_HT_ex_DF, 10)
fldm_HT_ex_DF <- head(fldm_HT_ex_DF, -5)
tail(fldm_HT_ex_DF, 5)

# Move gene names
# Pollen
row.names(pollen_ex_DF) <- pollen_ex_DF$X
pollen_ex_DF <- pollen_ex_DF[ ,-1]
# Fluidigm HT
row.names(fldm_HT_ex_DF) <- fldm_HT_ex_DF$X
fldm_HT_ex_DF <- fldm_HT_ex_DF[ ,-1]

## Clean chr number of ens IDs for bulk, Pollen
# Bulk
row.names(bulk_ex_DF) <- gsub("\\.[0-9]*", "", row.names(bulk_ex_DF))
# Pollen
row.names(pollen_ex_DF) <- gsub("\\.[0-9]*", "", row.names(pollen_ex_DF))
# Fluidigm LT
row.names(fldm_LT_ex_DF) <- gsub("\\.[0-9]*", "", row.names(fldm_LT_ex_DF))
# Fluidigm HT
row.names(fldm_HT_ex_DF) <- gsub("\\.[0-9]*", "", row.names(fldm_HT_ex_DF))

## Convert gene symbols to Ensembl IDs for Drop-seq
row.names(ds_ex_DF) <- Convert_Mixed_GeneSym_EnsID_To_EnsID(row.names(ds_ex_DF))

## Subset to GZ
# Bulk
gz_bulk_ex_DF <- bulk_ex_DF[ ,bulk_mt_DF$ExpCondition == "VZ"]
# Drop-seq
ids <- row.names(centSO@meta.data)[centSO@meta.data$REGION == "GZ"]
gz_df_ex_DF <- ds_ex_DF[, colnames(ds_ex_DF) %in% ids]
# Fluidigm HT
gz_fldm_HT_ex_DF <- fldm_HT_ex_DF[ ,grep(
  "COL01|COL02|COL03|COL04|COL05|COL06|COL07|COL08|COL09|COL10"
  , colnames(fldm_HT_ex_DF))]

# ## CQN Normalization for length
#
# # CQN fits this model:
# # log2(RPM) = s(x) + s(log2(length))
# # Since x is required, input log2(length) for x
#
# # Use CQN to normalize for GC content and gene length
# RunCQN <- function (exprDatDF, lengthGcDF) {
#   # Remove genes not in GC and gene length annotation file
#   keepV <- intersect(rownames(lengthGcDF), rownames(exprDatDF))
#   geneAnno <- lengthGcDF[match(keepV, rownames(lengthGcDF)), ]
#   # Set genes with length 0 to length 1 - why? - from Vivek's code
#   geneAnno[geneAnno[,1] == 0] <- 1
#   exprDatDF <- exprDatDF[match(keepV, rownames(exprDatDF)), ]
#
#   # Run CQN with specified depths and no quantile normalization
#   cqnDat <- cqn(exprDatDF, lengths = 1
#                 , x = log(as.numeric(geneAnno[ ,1]), 2)
#                 , lengthMethod = c("fixed")
#                 , sqn = FALSE)
#   # Get the log2(Normalized FPKM) values
#   cqnDat <- cqnDat$y + cqnDat$offset
#   cqnDat
# }
#
# bulk_ex_DF <- RunCQN(bulk_ex_DF, ENSEMBLhg19.70UnionAnno)
################################################################################

### Pools of different sizes and read depth normalize

## Pools of different sizes correlated to bulk

# Pools of different sizes correlations to bulk
# mean_bulk_ex_DF: mean expression dataframe
# sc_ex_DF: single-cell expression data frame
# nSamp: Vector of sizes of pools
Pool_Size_Correlation <- function(
  mean_bulk_ex_DF, sc_ex_DF, nSamp, Dataset_label = NA) {

  print(Dataset_label)

  # Do not exceed number of cells in dataset with sampling
  nSamp <- nSamp[nSamp < ncol(sc_ex_DF)]

  # Format
  # Merge bulk and pooled means, use union of genes
  combo_mean_ex_DF <- merge(mean_bulk_ex_DF, sc_ex_DF
    , by = "row.names", all = TRUE)
  names(combo_mean_ex_DF)[1] <- "ensembl_gene_id"
  combo_mean_ex_DF <- combo_mean_ex_DF[
    , ! names(combo_mean_ex_DF) == "hgnc_symbol"]
  # Convert NAs to 0s
  combo_mean_ex_DF[is.na(combo_mean_ex_DF)] <- 0

  # Sample, pool, and correlation
  sprman_cor <- sapply(nSamp, function(i) {
    print(i)
    # Sample
    sampled_DF <- combo_mean_ex_DF[
      , c(1:2, sample(3:ncol(combo_mean_ex_DF), i, replace = FALSE))]
    # Pool
    sampled_DF <- data.frame(sampled_DF[1:2]
      , Pooled = apply(sampled_DF[3:ncol(sampled_DF)], 1, sum)
    )
    # # Read depth normalize pooled by number mapped to exons
    # total_counts <- sum(plExDF$POOLED, na.rm = TRUE) / 10^6
    # Convert NAs to 0s
    sampled_DF[is.na(sampled_DF)] <- 0
    # # Read depth normalize
    # plExDF$POOLED <- plExDF$POOLED / total_counts

    # Spearman correlation
    round(cor(sampled_DF$Mean_Bulk, sampled_DF$Pooled, method = "spearman"), 2)
  })

  # Format ouput correlation data frame
  sprman_cor_DF <- data.frame(Number_of_cells = nSamp, Spearman = sprman_cor)
  sprman_cor_DF$Dataset_label <- Dataset_label
  row.names(sprman_cor_DF) <- paste0(Dataset_label, row.names(sprman_cor_DF))

  return(sprman_cor_DF)
}


Read_Depth_Normalize_To_CPM <- function(exM) {
  # Read depth normalize bulk by number mapped to exons
  # All cells
  total_counts <- (apply(exM, 2, sum) / 10^6)
  exM <- exM / total_counts
  return(exM)
}

Pool_Size_Correlation_Run <- function(){

  # Read depth normalize bulk by number mapped to exons
  # CP GZ
  bulk_ex_DF <- Read_Depth_Normalize_To_CPM(bulk_ex_DF)
  # GZ
  gz_bulk_ex_DF <- Read_Depth_Normalize_To_CPM(gz_bulk_ex_DF)

  # Bulk mean CPM
  # CP GZ
  mean_ex_DF <- data.frame(
    Mean_Bulk = apply(bulk_ex_DF, 1, mean)
    , Mean_Bulk_GZ = apply(gz_bulk_ex_DF, 1, mean)
  )

  # Pool sizes
  nSamp <- c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
    , seq(1000, 40000, 500))

  pool_corr_DF <- rbind(
    # Drop-seq
    Pool_Size_Correlation(
      mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk")]
      , sc_ex_DF = ds_ex_DF
      , nSamp = nSamp
      , Dataset_label = "Dropseq"
    )
    # Drop-seq GZ
    , Pool_Size_Correlation(
        mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk_GZ")]
        , sc_ex_DF = gz_df_ex_DF
        , nSamp = nSamp
        , Dataset_label = "Dropseq GZ"
    )
    # Pollen
    , Pool_Size_Correlation(
        mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk_GZ")]
        , sc_ex_DF = pollen_ex_DF
        , nSamp = nSamp
        , Dataset_label = "Pollen"
    )
    # Fluidigm HT
    , Pool_Size_Correlation(
        mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk")]
        , sc_ex_DF = fldm_HT_ex_DF
        , nSamp = nSamp
        , Dataset_label = "Fluidigm HT"
    )
    # Fluidigm HT GZ
    , Pool_Size_Correlation(
        mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk_GZ")]
        , sc_ex_DF = gz_fldm_HT_ex_DF
        , nSamp = nSamp
        , Dataset_label = "Fluidigm HT GZ"
    )
    # Fluidigm LT
    , Pool_Size_Correlation(
        mean_bulk_ex_DF = mean_ex_DF[c("Mean_Bulk")]
        , sc_ex_DF = fldm_LT_ex_DF
        , nSamp = nSamp
        , Dataset_label = "Fluidigm LT"
    )
  )
  return(pool_corr_DF)
}

# Number of genes detected as pool size increases
Number_Genes_Detected <- function(exDF, nSamp, datasetName) {
  # Genes with >= 1 count (raw counts)
  nGeneCutoff <- c(0,1,2,10)
  names(nGeneCutoff) <- paste0("Detection_threshold_", nGeneCutoff)
  nGeneL <- lapply(nGeneCutoff, function(nGeneCutoff){
    print(paste0("Genes with >= ", nGeneCutoff, " count (raw counts)"))
    tfExDF <- exDF >= nGeneCutoff
    nGene <- sapply(nSamp, function(i) {
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
    return(nGene)
  })
  # Format into data drame
  nGeneDF <- as.data.frame(do.call("cbind", nGeneL))
  nGeneDF$Number_of_cells <- nSamp
  nGeneDF$Dataset <- datasetName

  return(nGeneDF)
}

### Number of genes detected as pool size increases
Number_Genes_Detected_Run <- function(){
  ngd_DF <- rbind(
    # Drop-seq
    Number_Genes_Detected(ds_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
        , seq(1000, (ncol(ds_ex_DF)), 500))
      , "DROPSEQ")
    # Fluidigm HT
    , Number_Genes_Detected(fldm_HT_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(fldm_HT_ex_DF)), 100))
      , "FLUIDIGM_HT")
    # Fluidigm LT
    , Number_Genes_Detected(fldm_LT_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(fldm_LT_ex_DF)), 50))
      , "FLUIDIGM_LT")
    # Drop-seq GZ
    , Number_Genes_Detected(gz_df_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100)
        , seq(1000, (ncol(gz_df_ex_DF)), 500))
      , "DROPSEQ_GZ")
    # Fluidigm HT GZ
    , Number_Genes_Detected(gz_fldm_HT_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(gz_fldm_HT_ex_DF)), 100))
      , "FLUIDIGM_HT_GZ")
    , Number_Genes_Detected(pollen_ex_DF
      , c(seq(1, 9, 1), seq(10, 90, 10), seq(100, (ncol(pollen_ex_DF)), 100))
      , "POLLEN")
  )
  return(ngd_DF)
}

pool_corr_DF <- Pool_Size_Correlation_Run()
number_gene_detected_DF <- Number_Genes_Detected_Run()
################################################################################

### Pool scRNA-seq and read depth normalize

# Randomly sample and pool cells
# rNums: Columns to sample from
# nPool: Number of pools to split into
Pool_Cells <- function(exDF, nPool, rNums) {
  rNums <- rNums[rNums < ncol(exDF)]
  rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / nPool)))
  pExDF <- data.frame(lapply(rndmGroups, function (group){
    apply(exDF[ ,group], 1, sum)
  }))
  print(head(pExDF, 10))
  return(pExDF)
}

Pool_To_CPM <- function(exM, nPools, rNums){
  exM <- Pool_Cells(exM, nPools, rNums)
  total_counts <- (apply(exM, 2, sum) / 10^6)
  cpm_DF <- exM / total_counts
  return(cpm_DF)
}

Pool_To_CPM_Run <- function(){
  pools_cpm_DFL <- list(
    # Dropseq
    # 310 pooled groups of 100 cells each
    Dropseq = Pool_To_CPM(
      exM = ds_ex_DF
      , nPools = 310
      , rNums = sample(1:31000, 31000, replace = FALSE)
    )
    # Dropseq GZ
    # 150 pooled groups of 100 cells each
    , Dropseq_GZ = Pool_To_CPM(
      exM = gz_df_ex_DF
      , nPools = 150
      , rNums = sample(1:15000, 15000, replace = FALSE)
    )
    # Fluidigm HT
    # 12 pooled groups of 48 cells each
    , Fluidigm_HT = Pool_To_CPM(
      exM = fldm_HT_ex_DF
      , nPools = 12
      , rNums <- sample(1:48, 48, replace = FALSE)
    )
    # Fluidigm HT GZ
    # 6 pooled groups of 48 cells each
    , Fluidigm_HT_GZ = Pool_To_CPM(
      exM = gz_fldm_HT_ex_DF
      , nPools = 6
      , rNums <- sample(1:288, 288, replace = FALSE)
    )
    # Pollen
    # 10 pooled groups of 39 cells each
    , Pollen = Pool_To_CPM(
      exM = pollen_ex_DF
      , nPools = 10
      , rNums = sample(1:390, 390, replace = FALSE)
    )
    # Fluidigm LT
    # 6 pooled groups of 30 cells each
    , Fluidigm_LT = Pool_To_CPM(
      exM = fldm_LT_ex_DF
      , nPools = 6
      , rNums = sample(1:180, 180, replace = FALSE)
    )
  )
  return(pools_cpm_DFL)
}

pools_cpm_DFL <- Pool_To_CPM_Run()
################################################################################

### Mean CPM for each pooled dataset and bulk

Mean_CPM <- function(){
  # Mean counts
  mean_cpm_DFL <- lapply(pools_cpm_DFL, function(pools_cpm_DF) {
    mean_cpm <- apply(pools_cpm_DF, 1, mean)
    return(data.frame(mean_cpm))
  })

  # Read depth normalize bulk by number mapped to exons
  # CP GZ
  bulk_ex_DF <- Read_Depth_Normalize_To_CPM(bulk_ex_DF)
  # GZ
  gz_bulk_ex_DF <- Read_Depth_Normalize_To_CPM(gz_bulk_ex_DF)
  # Bulk mean CPM
  mean_cpm_DFL$Bulk = data.frame(mean_cpm = apply(bulk_ex_DF, 1, mean))
  mean_cpm_DFL$Bulk_GZ = data.frame(mean_cpm = apply(gz_bulk_ex_DF, 1, mean))

  # Merge keeping all genes (union)
  mean_cpm_DF <- mean_cpm_DFL[[1]]
  names(mean_cpm_DF) <- names(mean_cpm_DFL)[1]
  mean_cpm_DF$ensembl_gene_id <- row.names(mean_cpm_DF)
  for (name in names(mean_cpm_DFL[-1])) {
    print(name)
    df1 <- mean_cpm_DFL[[name]]
    names(df1) <- name
    df1$ensembl_gene_id <- row.names(df1)
    mean_cpm_DF <- merge(mean_cpm_DF, df1
      , by = "ensembl_gene_id", all = TRUE
    )
  }
  # Convert NAs to 0s
  mean_cpm_DF[is.na(mean_cpm_DF)] <- 0

  return(mean_cpm_DF)
}

Subset_To_Biomart_EnsIDs <- function(){
  mean_cpm_DF <- mean_cpm_DF[
    mean_cpm_DF$ensembl_gene_id %in% bmDF$ensembl_gene_id, ]
  return(mean_cpm_DF)
}

Add_Gene_Info <- function(){
  # Subset to ensembl IDs in biomart and add gene information
  idx <- match(mean_cpm_DF$ensembl_gene_id, bmDF$ensembl_gene_id)
  mean_cpm_DF$ensembl_gene_id <- bmDF$ensembl_gene_id[idx]
  mean_cpm_DF$gene_biotype <- bmDF$gene_biotype[idx]
  mean_cpm_DF$percentage_gc_content <- bmDF$percentage_gc_content[idx]
  mean_cpm_DF$cds_length <- bmDF$cds_length[idx]

  return(mean_cpm_DF)
}

mean_cpm_DF <- Mean_CPM()
mean_cpm_DF <- Subset_To_Biomart_EnsIDs()
mean_cpm_DF <- Add_Gene_Info()

save(list = ls(), file = paste0(outAnalysis, "Workspace_TEST.RData"))
################################################################################

### Subset to higher in pooled or higher in bulk

# Log2(mean CPM + 1)
Log2_Mean_CPM <- function(){
  log2_mean_cpm_DF <- apply(mean_cpm_DF[ ,2:9], 2, function(y){
    log2_mean_cpm <- log(y + 1, 2)
    return(log2_mean_cpm)
  })
  log2_mean_cpm_DF <- as.data.frame(log2_mean_cpm_DF)
  row.names(log2_mean_cpm_DF) <- mean_cpm_DF$ensembl_gene_id
  return(log2_mean_cpm_DF)
}

Ratio_Log2_Mean_CPM <- function(){
  # Ratio of single-cell vs bulk
  df1 <- apply(log2_mean_cpm_DF[ ,c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")], 2
    , function(y){
      ratio <- (y - log2_mean_cpm_DF[ ,"Bulk"])
      return(ratio)
  })
  df2 <- apply(log2_mean_cpm_DF[ ,c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen")], 2
    , function(y){
      ratio <- y - log2_mean_cpm_DF[ ,"Bulk_GZ"]
      return(ratio)
  })
  ratio_log2_mean_cpm_DF <- cbind(df1, df2)
  ratio_log2_mean_cpm_DF <- as.data.frame(ratio_log2_mean_cpm_DF)
  row.names(ratio_log2_mean_cpm_DF) <- row.names(log2_mean_cpm_DF)
  return(ratio_log2_mean_cpm_DF)
}

# Flag as high in single-cell or bulk
Bias_Flag_DF <- function(){

  bias_flag_DF <- apply(ratio_log2_mean_cpm_DF, 2, function(y){
    bias_DF <- Subset_Two_SD(y)
    bias <- rep("None", nrow(bias_DF))
    bias[bias_DF$High_Bulk == TRUE] <- "Bulk high"
    bias[bias_DF$High_Pool == TRUE] <- "Pool high"
    return(bias)
  })
  bias_flag_DF <- as.data.frame(bias_flag_DF)
  row.names(bias_flag_DF) <- row.names(ratio_log2_mean_cpm_DF)

  # Flag 0 detection for drop-seq
  bias_flag_DF$Dropseq_0 <- "None"
  bias_flag_DF$Dropseq_0[
    log2_mean_cpm_DF$Bulk == 0 & bias_flag_DF$Dropseq == "Pool high"
    ] <- "Pool high"
  bias_flag_DF$Dropseq_0[
    log2_mean_cpm_DF$Dropseq == 0 & bias_flag_DF$Dropseq == "Bulk high"
    ] <- "Bulk high"
  bias_flag_DF$Dropseq_0 <- as.factor(bias_flag_DF$Dropseq_0)

  # Intersect Dropseq GZ, Fluidigm HT GZ, Pollen
  bias_flag_DF$GZ_Datasets <- "None"
  bias_flag_DF$GZ_Datasets[
    bias_flag_DF$Dropseq == "Pool high" &
    bias_flag_DF$Fluidigm_HT == "Pool high" &
    bias_flag_DF$Fluidigm_LT == "Pool high"
    ] <- "Pool high"
  bias_flag_DF$GZ_Datasets[
    bias_flag_DF$Dropseq == "Bulk high" &
    bias_flag_DF$Fluidigm_HT == "Bulk high" &
    bias_flag_DF$Fluidigm_LT == "Bulk high"
    ] <- "Bulk high"
  bias_flag_DF$GZ_Datasets <- as.factor(bias_flag_DF$GZ_Datasets)

  # Intersect Dropseq GZ, Fluidigm HT GZ, Pollen
  bias_flag_DF$GZCP_Datasets <- "None"
  bias_flag_DF$GZCP_Datasets[
    bias_flag_DF$Dropseq_GZ == "Pool high" &
    bias_flag_DF$Fluidigm_HT_GZ == "Pool high" &
    bias_flag_DF$Pollen == "Pool high"
    ] <- "Pool high"
  bias_flag_DF$GZCP_Datasets[
    bias_flag_DF$Dropseq_GZ == "Bulk high" &
    bias_flag_DF$Fluidigm_HT_GZ == "Bulk high" &
    bias_flag_DF$Pollen == "Bulk high"
    ] <- "Bulk high"
  bias_flag_DF$GZCP_Datasets <- as.factor(bias_flag_DF$GZCP_Datasets)

  return(bias_flag_DF)
}

log2_mean_cpm_DF <- Log2_Mean_CPM()
ratio_log2_mean_cpm_DF <- Ratio_Log2_Mean_CPM()
bias_flag_DF <- Bias_Flag_DF()

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))
################################################################################

### Plots

## Correlation to bulk as pool size increases
ggplot(pool_corr_DF, aes(
  x = Number_of_cells, y = Spearman, color = Dataset_label)) +
  facet_wrap(~Dataset_label, ncol = 2, scales = "free") +
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

## Number of genes detected as pool size increases
ggDF <- melt(number_gene_detected_DF
  , measure.vars = c("Detection_threshold_0"
  , "Detection_threshold_1", "Detection_threshold_2"
  , "Detection_threshold_10")
)
# Plot as grid
ggplot(ggDF, aes(x = Number_of_cells, y = value, col = variable)) +
  facet_wrap(~DATASET, ncol = 3, scales = "free_x") +
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

## Scatter plots of pooled vs bulk with fit and spearman
# Format for ggplot
df1 <- melt(mnExDF[
      c("Bulk", "Bulk_GZ", "Dropseq", "Fluidigm_HT", "Fluidigm_LT")
      ]
    , measure.vars = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
  )
df2 <- melt(mnExDF[
    c("Bulk", "Bulk_GZ", "Dropseq_GZ" "Pollen", "Fluidigm_HT_GZ")
  ]
  , measure.vars = c("Dropseq_GZ", "Pollen", "Fluidigm_HT_GZ")
)
df2$Bulk <- df2$Bulk_GZ
ggDF <- rbind(df1, df2)
# Spearman
sprL <- lapply(split(ggDF, ggDF$variable), function(df) {
  cor(df$Bulk, df$value, method = "spearman")})
sprL <- lapply(sprL, round, 2)
ggDF$value <- log(ggDF$value + 1, 2)
ggDF$Bulk <- log(ggDF$Bulk + 1, 2)
labelsDatasetSpr <- c(paste0(unique(ggDF$variable), "\nSpearman: ", sprL))
names(labelsDatasetSpr) <- unique(ggDF$variable)
# ggplot
ggplot(ggDF, aes(x = value, y = Bulk, color = variable)) +
  facet_wrap("variable", ncol = 2
    , labeller = labeller(variable = labelsDatasetSpr)) +
  geom_point(alpha = 0.15, shape = 1, size = 0.25) +
  stat_smooth(col = "black") +
  ylab("Bulk: log2(mean CPM + 1)") +
  xlab("Pooled: log2(mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nPooled scRNA-seq vs Bulk RNA-seq"
    , "\n"
    , "\nMean of CPM across samples and pools"))
ggsave(paste0(outGraph, "Scatter.png"), dpi = 600, width = 6.5, height = 12)

## Table of number of genes biased toward bulk or single-cell
dfl1 <- lapply(names(bias_flag_DF), function(name){
  y <- bias_flag_DF[name]
  df1 <- data.frame(table(y))
  row.names(df1) <- df1$y
  df1 <- df1[ -1]
  names(df1) <- name
  return(df1)
})
df1 <- do.call("cbind", dfl1)
write.csv(df, file = paste0(outTable, "_NumberGenesSubset.csv"), quote = FALSE)

## Plot subset scatter plots and density plots
ldf1 <- lapply(c("Dropseq", "Fluidigm_HT", "Fluidigm_LT"), function(name){
  data.frame(
    Mean_Pooled = log2_mean_cpm_DF[[name]]
    , Mean_Bulk = log2_mean_cpm_DF$Bulk
    , Subset = bias_flag_DF[[name]]
    , Dataset = name
  )
})
ldf2 <- lapply(c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen"), function(name){
  data.frame(
    Mean_Pooled = log2_mean_cpm_DF[[name]]
    , Mean_Bulk = log2_mean_cpm_DF$Bulk_GZ
    , Subset = bias_flag_DF[[name]]
    , Dataset = name
  )
})
ggDF <- rbind(do.call("rbind", ldf1)
  , do.call("rbind", ldf2)
)
# scatter
ggplot(ggDF, aes(x = Mean_Pooled, y = Mean_Bulk)) +
  facet_wrap("Dataset", ncol = 2) +
  geom_point(alpha = 0.15, shape = 16, size = 0.25, aes(color = Subset)) +
  scale_colour_manual(values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ylab("Bulk: log2(Mean CPM + 1)") +
  xlab("Pooled: log2(Mean CPM + 1)") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nSubsetting pooled scRNAseq vs bulk RNAseq"
    , "\n"
    , "\nHuman fetal brain GZ and CP"
    , "\nMean of CPM across samples and pools"
    , "\nSubset 2 SD from mean ratio of pooled vs bulk"))
ggsave(paste0(outGraph, "ScatterSubset.png"), dpi = 150, width = 6, height = 9)
# density plot
ggDF <- melt(ggDF)
ggplot(ggDF, aes(x = value)) +
  facet_grid(variable ~ Dataset) +
  geom_density(aes(col = Subset)) +
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

## Venn diagrams of gene intersections high in bulk or high in pooled
pdf(paste0(outGraph, "VennDiagrams.pdf"))
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
    , bias_direction = "Bulk high"
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq", "Fluidigm_HT", "Fluidigm_LT")
    , bias_direction = "Pool high"
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen")
    , bias_direction = "Bulk high"
  )
  Venn_Diagram_3_Wrapper(
    datasets = c("Dropseq_GZ", "Fluidigm_HT_GZ", "Pollen")
    , bias_direction = "Pool high"
  )
dev.off()

# ## Drop-seq union genes versus intersection genes
# # Intersection (only genes in both)
# mnBothDF <- merge(mean_bulk_ex_DF, mnPods_ex_DF, by.x = "row.names"
#   , by.y = "ensembl_gene_id")
# # Convert NAs to 0s
# mnBothDF[is.na(mnBothDF)] <- 0
# # Spearman correlation
# # Union
# sprUn <- round(cor(mnExDF$Mean_Bulk, mnExDF$MEAN_POOLED_DROPSEQ
#   , method = "spearman"), 2)
# sprUn
# # Intersection
# sprIn <- round(cor(mnBothDF$Mean_Bulk, mnBothDF$MEAN_POOLED
#   , method = "spearman"), 2)
# sprIn
# # Format for ggplot2
# df1 <- data.frame(Pooled = log(mnExDF$Mean_Bulk + 1, 2)
#   , Bulk = log(mnExDF$MEAN_POOLED_DROPSEQ + 1, 2))
# df1$DATASET <- "Union genes"
# df2 <- data.frame(Pooled = log(mnBothDF$Mean_Bulk + 1, 2)
#   , Bulk = log(mnBothDF$MEAN_POOLED + 1, 2))
# df2$DATASET <- "Intersection genes"
# ggDF <- rbind(df1, df2)
# # Add spearman to labels
# sprL <- c(sprUn, sprIn)
# labelsDatasetSpr <- c(paste0(unique(ggDF$DATASET), "\nSpearman: ", sprL))
# names(labelsDatasetSpr) <- unique(ggDF$DATASET)
# # ggplot
# ggplot(ggDF, aes(x = Pooled, y = Bulk, col = DATASET)) +
#   facet_wrap("DATASET", ncol = 2
#     , labeller = labeller(DATASET = labelsDatasetSpr)) +
#   geom_point(alpha = 0.5, shape = 1, size = 0.25) +
#   stat_smooth(col = "black") +
#   scale_color_manual(values = c("#e41a1c", "#377eb8")) +
#   theme(legend.position = "none") +
#   ylab("Bulk: log2(Mean CPM + 1)") +
#   xlab("Pooled: log2(Mean CPM + 1)") +
#   ggtitle(paste0(graphCodeTitle
#     , "\n"
#     , "\nPooled Drop-seq vs Bulk RNAseq"
#     , "\n"
#     , "\nIntersection (genes only detected in both) or union (genes in either)"
#     , "\nHuman Fetal Brain CP + GZ"
#     , "\nMean of CPM across samples and pools"))
# ggsave(paste0(outGraph, "ScatterIntersectionUnion.png"), width = 8, height = 5
#   , dpi = 600)
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

### GC, length, MT genes, brain expressed, and biotype of gene intersections
### high in bulk or pooled

## CDS length of biased gene sets
# Format for ggplot
ggDF <- bias_flag_DF
ggDF$ensembl_gene_id <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "ensembl_gene_id")
idx <- match(ggDF$ensembl_gene_id, mean_cpm_DF$ensembl_gene_id)
ggDF$cds_length <- mean_cpm_DF$cds_length[idx]
ggDF <- ggDF[! is.na(ggDF$cds_length), ]
ggDF$value[ggDF$value == "None"] <- "All other genes"
ggDF$value <- as.factor(ggDF$value)
# Plot
ggplot(ggDF, aes(x = variable, y = cds_length, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 5000)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("CDS length") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "Length_Boxplot.pdf"), width = 8, height = 5)
# Subset for paper
ggDF <- ggDF[ggDF$variable %in% c("GZ_Datasets", "GZCP_Datasets"), ]
ggplot(ggDF, aes(x = variable, y = cds_length, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 5000)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("CDS length") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "Length_Boxplot_paper.pdf"), width = 5, height = 5)

# CDS length for each dataset
Plot_CDS_Length_Histogram <- function(exM, title, geneIdType){
  cGene <- rowSums(exM)
  ggDF <- data.frame(GENE = names(cGene), COUNTS = cGene)
  ggDF <- merge(ggDF, mean_cpm_DF[c("ensembl_gene_id", "cds_length")]
    , by.x = "GENE", by.y = "ensembl_gene_id")
  ggDF <- data.frame(CDS_LENGTH = rep(ggDF$cds_length, ggDF$COUNTS))
  gg <- ggplot(ggDF, aes(x = CDS_LENGTH)) +
    geom_histogram(binwidth = 500) +
    # coord_cartesian(xlim = c(0, 10000)) +
    # scale_x_continuous(limits = c(0, 10000)) +
    xlab("CDS length") +
    ylab("Counts") +
    ggtitle(title)
  return(gg)
}
ggL <- list(
  Plot_CDS_Length_Histogram(bulk_ex_DF, "Bulk RNA-seq")
  , Plot_CDS_Length_Histogram(ds_ex_DF, "Drop-seq")
  , Plot_CDS_Length_Histogram(gz_df_ex_DF, "Drop-seq GZ")
  , Plot_CDS_Length_Histogram(fldm_HT_ex_DF, "Fluidigm HT")
  , Plot_CDS_Length_Histogram(gz_fldm_HT_ex_DF, "Fluidigm HT GZ")
  , Plot_CDS_Length_Histogram(fldm_LT_ex_DF, "Fluidigm LT")
  , Plot_CDS_Length_Histogram(pollen_ex_DF, "Pollen")
)
Plot_Grid(ggL, ncol = 3, rel_height = 0.1, title = paste0(graphCodeTitle
  , "\n"
  , "\nHistogram of raw counts versus CDS length"
  , "\n")
)
ggsave(paste0(outGraph, "LengthCounts_Histogram.png"), width = 8, height = 20)

## GC content for each dataset
ggDF <- bias_flag_DF
ggDF$percentage_gc_content <- row.names(ggDF)
ggDF <- melt(ggDF, id.vars = "percentage_gc_content")
idx <- match(ggDF$percentage_gc_content, mean_cpm_DF$percentage_gc_content)
ggDF$cds_length <- mean_cpm_DF$cds_length[idx]
ggDF <- ggDF[! is.na(ggDF$cds_length), ]
ggDF$value[ggDF$value == "None"] <- "All other genes"
ggDF$value <- as.factor(ggDF$value)
# Plot
ggplot(ggDF, aes(x = variable, y = cds_length, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 5000)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("Percent GC") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "GC_boxplot.pdf"), width = 8, height = 5)
# Subset for paper
ggDF <- ggDF[ggDF$variable %in% c("GZ_Datasets", "GZCP_Datasets"), ]
ggplot(ggDF, aes(x = variable, y = cds_length, fill = value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Gene subset"
    , values = c("#e31a1c", "#bdbdbd", "#1f78b4")) +
  coord_cartesian(ylim = c(0, 100)) +
  # theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("Percent GC") +
  ggtitle(paste0(graphCodeTitle
    , "\n"
    , "\nCDS length of genes biased towards single-cell or bulk RNA-seq"
    , "\n")
  )
ggsave(paste0(outGraph, "GC_boxplot_paper.pdf"), width = 5, height = 5)

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
  pct <- sum(brain_exp_DF[ ,1] %in% genes) / length(genes)
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Bulk high", ]
  genes <- df$ENSEMBL_ID
  pct <- c(pct, sum(brain_exp_DF[ ,1] %in% genes) / length(genes))
  df <- ssl2MnCpmDF[ssl2MnCpmDF[[geneSet]] == "Pool high", ]
  genes <- df$ENSEMBL_ID
  pct <- c(pct, sum(brain_exp_DF[ ,1] %in% genes) / length(genes))
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

print("### GO and pathway analysis")

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

dir.create("../analysis/graphs/Pooled_Vs_Bulk/DS2-11/Pooled_Vs_Bulk_GO_CDS_Length", recursive = TRUE)
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
nc <- length(l)
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
nr <- max(sapply(ll, length))
nc <- 3*length(ll)
m <- matrix(NA, nr, nc)
for(i in 1:length(ll)){
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
# Label non-marker genes
df$Grouping <- as.character(df$Grouping)
df$Grouping[df$Grouping == ""] <- "Non-marker gene"
df$Grouping <- factor(df$Grouping
  , levels = c("Non-marker gene", as.character(unique(kmDF$Grouping))))

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
ggsave(paste0(outGraph, "KnownMarkersCombined.pdf"), width = 10, height = 7)
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
