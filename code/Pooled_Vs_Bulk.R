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
# nowakowski_2017_DF <- read.table("../nowakowski_2017/geneMatrix.tsv"
#   , header = TRUE)

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
dir.create(dirname(outTableGpro), recursive = TRUE)
outTableGproIds <- "../analysis/tables/Pooled_Vs_Bulk/DS2-11/gprofiler_goIDs/Pooled_Vs_Bulk_gprofiler_goIDs_"
dir.create(dirname(outTableGproIds), recursive = TRUE)
outTable_GOelite_background <- "/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/analysis/tables/Pooled_Vs_Bulk/DS2-11/GOElite/background"
outTable_GOelite_genes <- "/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/analysis/tables/Pooled_Vs_Bulk/DS2-11/GOElite/genes"
outTable_GOelite_output <- "/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/analysis/tables/Pooled_Vs_Bulk/DS2-11/GOElite/output"
dir.create(outTable_GOelite_background, recursive = TRUE)
dir.create(outTable_GOelite_genes, recursive = TRUE)
dir.create(outTable_GOelite_output, recursive = TRUE)
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

# Subset to higher in bulk or higher in pooled
Subset_Two_SD <- function(dataset) {
  idx <- data.frame(
    High_Bulk = dataset < (mean(dataset) - 2*sd(dataset))
    , High_Pool = dataset > (mean(dataset) + 2*sd(dataset))
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

  # # Merge bulk and pooled means, use union of genes
  # combo_mean_ex_DF <- merge(mean_bulk_ex_DF, sc_ex_DF
  #   , by = "row.names", all = TRUE)
  # Merge bulk and pooled means, use intersection of genes
  combo_mean_ex_DF <- merge(mean_bulk_ex_DF, sc_ex_DF
    , by = "row.names")
  # Format
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
  mean_cpm_DF$hgnc_symbol <- bmDF$hgnc_symbol[idx]

  return(mean_cpm_DF)
}

Mean_CPM_Run <- function(){
  mean_cpm_DF <- Mean_CPM()
  mean_cpm_DF <- Subset_To_Biomart_EnsIDs()
  mean_cpm_DF <- Add_Gene_Info()
  return(mean_cpm_DF)
}
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
################################################################################

### GO and pathway analysis

## gprofiler
Gprofiler_Run <- function(){
  print("Gprofiler_Run")
  gpro_results_DFLL <- lapply(bias_flag_DF, function(dataset){
    idx <- match(row.names(bias_flag_DF), mean_cpm_DF$ensembl_gene_id)
    dataset <- data.frame(Flag = dataset
      , ensembl_gene_id = row.names(bias_flag_DF)
      , hgnc_symbol = mean_cpm_DF$hgnc_symbol[idx]
    )
    # Filter genes with no hgnc_symbol
    dataset <- dataset[! dataset$hgnc_symbol == "", ]
    # Background gene list
    background_gene_list <- as.character(dataset$hgnc_symbol)
    # Bulk high
    genes_bulk_high <- as.character(dataset$hgnc_symbol[
      dataset$Flag == "Bulk high"])
    results1 <- gprofiler(
      query = genes_bulk_high
      , custom_bg = background_gene_list
    )
    # Pool high
    genes_pool_high <- as.character(dataset$hgnc_symbol[
      dataset$Flag == "Pool high"])
    results2 <- gprofiler(
      query = genes_pool_high
      , custom_bg = background_gene_list
    )
    # Output results
    gpro_results_L <- list("Bulk high" = results1, "Pool high" = results2)
    return(gpro_results_L)
  })
  return(gpro_results_DFLL)
}

## Run GO elite as nohupped shell script:
Run_GOelite <- function(){
  print("Run_GOelite")

  # Intersected genes from GZCP datasets
  genes <- row.names(bias_flag_DF)[bias_flag_DF["GZCP_Datasets"] == "Bulk high"]
  background_genes <- row.names(bias_flag_DF)
  genes_path <- paste0(outTable_GOelite_genes, "/GZCP_Datasets/genes.txt")
  dir.create(dirname(genes_path), recursive = TRUE)
  background_path <- paste0(outTable_GOelite_background
    , "/GZCP_Datasets/genes.txt")
  dir.create(dirname(background_path), recursive = TRUE)
  goelite_output_path <- paste0(outTable_GOelite_output, "/GZCP_Datasets")
  dir.create(goelite_output_path, recursive = TRUE)
  write.table(genes, file = genes_path, quote = FALSE, row.names = FALSE
    , col.names = FALSE)
  write.table(background_genes, file = background_path, quote = FALSE
    , row.names = FALSE, col.names = FALSE)
  nperm <- as.integer(10) # or 50000
  system(paste0(
    "python ~/bin/GO-Elite_v.1.2.5-Py/GO_Elite.py --species Hs --mod Ensembl --permutations ", nperm, "  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",dirname(genes_path)," --denom ", dirname(background_path), " --output ", goelite_output_path," &", sep=""))

  # Intersected genes from GZ datasets
  genes <- row.names(bias_flag_DF)[bias_flag_DF["GZ_Datasets"] == "Bulk high"]
  background_genes <- row.names(bias_flag_DF)
  genes_path <- paste0(outTable_GOelite_genes, "/GZ_Datasets/genes.txt")
  dir.create(dirname(genes_path), recursive = TRUE)
  background_path <- paste0(outTable_GOelite_background
    , "/GZ_Datasets/genes.txt")
  dir.create(dirname(background_path), recursive = TRUE)
  goelite_output_path <- paste0(outTable_GOelite_output, "/GZ_Datasets")
  dir.create(goelite_output_path, recursive = TRUE)
  write.table(genes, file = genes_path, quote = FALSE, row.names = FALSE
    , col.names = FALSE)
  write.table(background_genes, file = background_path, quote = FALSE
    , row.names = FALSE, col.names = FALSE)
  nperm <- as.integer(10) # or 50000
  system(paste0(
    "python ~/bin/GO-Elite_v.1.2.5-Py/GO_Elite.py --species Hs --mod Ensembl --permutations ", nperm, "  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",dirname(genes_path)," --denom ", dirname(background_path), " --output ", goelite_output_path," &", sep=""))
}
################################################################################

### Process data

pool_corr_DF <- Pool_Size_Correlation_Run()
number_gene_detected_DF <- Number_Genes_Detected_Run()
pools_cpm_DFL <- Pool_To_CPM_Run()
mean_cpm_DF <- Mean_CPM_Run()
log2_mean_cpm_DF <- Log2_Mean_CPM()
ratio_log2_mean_cpm_DF <- Ratio_Log2_Mean_CPM()
bias_flag_DF <- Bias_Flag_DF()
gpro_results_DFLL <- Gprofiler_Run()

save(pool_corr_DF
  , number_gene_detected_DF
  , pools_cpm_DFL
  , mean_cpm_DF
  , log2_mean_cpm_DF
  , ratio_log2_mean_cpm_DF
  , bias_flag_DF
  , gpro_results_DFLL
  , file = paste0(outAnalysis, "Processed_Data.RData")
)

Run_GOelite()
################################################################################
