# Damon Polioudakis
# 2016-12-06
# Plot read depth versus number of genes detected for DS-002 and HT-003
# downsampled read depth
################################################################################

require(ggplot2)

# DS-002
dirs <- list.files("../DS-002_DP/data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1")
# dir <- c("N702")

# Fluidigm HT-003
ht08DF <- read.table("../HT-003_DP/data/htseq/GRCh37.75_NoERCC/downsampleReads/Exprs_HTSCexon_downSample08p_FtMm3e4_ens.csv"
  , header= TRUE)
colnames(ht08DF) <- gsub("_downSample08p", "", colnames(ht08DF))
metDF <- read.table("../HT-003_DP/metadata/GRCh37.75_NoERCC/RNAstar_Stats.txt"
  , header = TRUE, fill = TRUE, sep = "\t")

dir.create("../analysis/graphs", recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 10))
################################################################################

## Combine all DS-002 samples
ngrLDF <- lapply(dirs, function(dir) {
  # Output from Reads_Per_Cell_Histogram.sh scripts
  # Using output for map quality 10 reads (I think about equivalent to uniquely mapped reads)
  # Select 800 most numerous cell barcodes (corresponding to 800 expected cells per sample)
  rdDF <- read.table(paste0("../DS-002_DP/data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/"
    , dir
    , "/out_cell_readcounts_readQ10.txt.gz")
    , header = FALSE, stringsAsFactors = F, nrows = 800)
  colnames(rdDF) <- c("READ_DEPTH", "CELL_BARCODE")
  print(head(rdDF))
  print(dim(rdDF))
  print(mean(rdDF$READ_DEPTH))
  print(sum(rdDF$READ_DEPTH))
  
  # Load genes detected statistics from digital gene expression script
  dgeDF <- read.table(paste0("../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/"
    , dir
    , "/out_gene_exon_tagged_dge.summary_FtMm250.txt")
    , header = TRUE, stringsAsFactors=F)
  print(head(dgeDF))
  
  df <- merge(rdDF, dgeDF, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
  print(head(df))
  return(df)
})
ngrDF <- do.call(rbind, ngrLDF)
ngrDF <- ngrDF[order(ngrDF$READ_DEPTH), ]
ngrDF <- ngrDF[ ,-4]
colnames(ngrDF) <- c("CELL_ID", "READ_DEPTH", "NUM_GENES")
ngrDF <- data.frame(ngrDF$CELL_ID, ngrDF$NUM_GENES, ngrDF$READ_DEPTH)
colnames(ngrDF) <- c("CELL_ID", "NUM_GENES", "READ_DEPTH")
ngrDF$EXPERIMENT <- "Drop-seq DS-002"


## HT-003
Calc_Genes_Detected_Fluidigm <- function (exDF) {
  nGenesHtDF <- data.frame(colSums(exDF[ ,-c(1:2)] >= 1))
  colnames(nGenesHtDF) <- c("NUM_GENES")
  return(nGenesHtDF)
}
# Genes detected
ngrHtDF <- Calc_Genes_Detected_Fluidigm(ht08DF)
# Add Uniquely.mapped.reads.number for each capture well
ngrHtDF <- merge(ngrHtDF, metDF[ ,c(1,5)], by.x = "row.names", by.y = "SampleID")
colnames(ngrHtDF) <- c("CELL_ID", "NUM_GENES", "READ_DEPTH")
ngrHtDF$READ_DEPTH <- ngrHtDF$READ_DEPTH * 0.08
ngrHtDF$EXPERIMENT <- "Fluidigm HT-003 downsampled aligned reads to 8%"


## Combine DS-002 and HT-003 number genes read depth data frames for ggplot2
ggDF <- rbind(ngrDF, ngrHtDF)


# ggplot with adjusted x-axis limits
ggplot(ggDF, aes(x = READ_DEPTH, y = NUM_GENES, color = EXPERIMENT)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  xlab("Aligned reads") +
  ylab("Number of genes detected") +
  coord_cartesian(xlim = c(0, 10^5), ylim = c(0, 5000)) +
  # theme(plot.margin = unit(c(0.2,1,0.2,0.2), "in")) +
  ggtitle(paste0("Read_Depth_Vs_Genes_Detected_Combined_Samples.R"
    , "\nX-limit adjusted"
    , "\n"))
ggsave("../analysis/graphs/Read_Depth_Vs_Genes_Detected_Combined_Samples.pdf"
  , width = 14, height = 6)