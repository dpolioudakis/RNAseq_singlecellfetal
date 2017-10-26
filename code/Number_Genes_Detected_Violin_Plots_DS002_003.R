# Damon Polioudakis
# 2016-11-20
# Plot number of genes detected from Drop-seq digital gene expression
# after down-sampling number of reads or aligning to exons + introns

# On hoffman2:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(reshape2)

## Load data and assign variables

## Inputs

# DS-003 digital gene expression parent directory
inPntDir <- "../DS-002_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1"
# DS-002 digital gene expression sample directories
inDirsL <- list.files(inPntDir)
# DS-002 digital gene expression matrices file names
inFilesL <- list.files(paste0(inPntDir, "/", inDirsL[1]))
inFilesL <- inFilesL[grep("out_gene_exon_tagged_dge_FtMm250.txt", inFilesL)]
# Make dataframe of sample folder names and digital gene expression file names
inPathsDF <- data.frame(lapply(inDirsL
  , function(inDir) paste0(inPntDir, "/", inDir, "/", inFilesL))
  , stringsAsFactors = FALSE)
# Remove N705
inPaths02DF <- inPathsDF[ ,-5]
colnames(inPaths02DF) <- c("1", "2", "3", "4")

# DS-003 digital gene expression parent directory
inPntDir <- "../DS-003_DP/data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsVAP048L8"
# DS-002 digital gene expression sample directories
inDirsL <- list.files(inPntDir)
# DS-002 digital gene expression matrices file names
inFilesL <- list.files(paste0(inPntDir, "/", inDirsL[1]))
inFilesL <- inFilesL[grep("out_gene_exon_tagged_dge_FtMm250.txt", inFilesL)]
# Make dataframe of sample folder names and digital gene expression file names
inPaths03DF <- data.frame(lapply(inDirsL
  , function(inDir) paste0(inPntDir, "/", inDir, "/", inFilesL))
  , stringsAsFactors = FALSE)
colnames(inPaths03DF) <- c("1", "2", "3", "4")

# Combine DS-002 and DS-003
inPathsDF <- cbind(inPaths02DF, inPaths03DF)

# Make output graphs directory
dir.create("../analysis/graphs", recursive = TRUE)

## Variables
graphCodeTitle <- "Number_Genes_Detected_Violin_Plots_DS002_003.R"
outGraphPfx <- "../analysis/graphs/Number_Genes_Detected_Violin_Plots_DS002_003"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 12))
################################################################################

## Load digital gene expression matrices and calculate number of gene detected
##  per cell

# Function to load digital expression matrix and sum genes for each cell with
# >= 1 count
# Input is path to digital expession matrix
# Outputs dataframe with cell barcodes as row names and number genes detected
# column
Calc_Genes_Detected <- function (inDeDF) {
  if (grepl("txt", inDeDF)) {
    deDF <- read.table(inDeDF, header = TRUE, stringsAsFactors = FALSE
      , row.names = 1)
  }
  if (grepl("csv", inDeDF)) {
    deDF <- read.csv(inDeDF, header = TRUE, stringsAsFactors = FALSE
      , row.names = 1)
  }
  # Remove gene name column
  deDF <- deDF[ ,-1]
  print(dim(deDF))
  return(data.frame(colSums(deDF >= 1)))
}

# Loop down rows of data frame with paths to digital expression matrices
# Rows are different filters or refFlats
# Columns are different samples (N701, N702, etc)
# Outputs list of data frames for each filter type with genes detected per cell
# for all samples. Each sample has same number of cells across filters because
# same Hs Mm alignment was used to remove cells
nGenesLDF <- apply(inPathsDF, 1, function(inDeL) {
  print(basename(inDeL[1]))
  # Loop through each column (sample N701, N702, etc) of row
  deLDF <- lapply(inDeL, function(inDe) Calc_Genes_Detected(inDe))
  str(deLDF)
  # Rbind samples (each row is genes detected per cell for a cell for a sample)
  deDF <- do.call("rbind", deLDF)
  # Name column with filter type
  colnames(deDF) <- basename(inDeL[1])
  return(deDF)
})

## Format for ggplot2 and calculate mean genes detected

## DS-002
# Combine list data frames for each filter type into 1 data frame
nGenesDF <- do.call("cbind", nGenesLDF)
colnames(nGenesDF) <- c("DS-002", "DS-003")

# Format for ggplot2
ggDF <- melt(nGenesDF, variable.name = "SAMPLE", value.name = "NUM_GENES")
# 
# ## HT-003
# Calc_Genes_Detected_Fluidigm <- function (exDF, sampleID, ggDF) {
#   exDF <- exDF[ ,-c(1:2)]
#   nGenesHtDF <- data.frame(colSums(exDF >= 1))
#   nGenesHtDF$SAMPLE <- sampleID
#   colnames(nGenesHtDF) <- c("NUM_GENES", "SAMPLE")
#   ggDF <- rbind(ggDF, nGenesHtDF)
# }
# # HT-003
# ggDF <- Calc_Genes_Detected_Fluidigm(inHtDF, "Fluidigm HT-003", ggDF)
# # HT-003 downsampled read depth to 8%
# ggDF <- Calc_Genes_Detected_Fluidigm(inHt08DF, "Fluidigm HT-003 downsample 8%", ggDF)
# # HT-003 downsampled read depth to 49%
# ggDF <- Calc_Genes_Detected_Fluidigm(inHt49DF, "Fluidigm HT-003 downsample 49%", ggDF)

# Mean genes detected for each filter
mnDF <- aggregate(NUM_GENES ~ SAMPLE, ggDF, mean)
colnames(mnDF) <- c("SAMPLE", "MEAN_NUMBER_GENES_DETECTED")
mnDF$MEAN_NUMBER_GENES_DETECTED <- round(mnDF$MEAN_NUMBER_GENES_DETECTED, 2)


## ggplot2

# Violin plots of number of genes detected for each sample after UMI collapse
ggplot(ggDF, aes(x = SAMPLE, y = NUM_GENES, fill = SAMPLE)) +
  geom_violin() +
  # geom_text(data = mnDF
  #   , aes(label = paste("Mean:\n", MEAN_NUMBER_GENES_DETECTED), y = 3000)
  #   , col = "black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.2,0.6,0.2,0.2), "in")) +
  xlab("Sample") +
  ylab("Number of genes detected") +
  ggtitle(paste0(graphCodeTitle
    , "\nViolin plot of number of genes with >=1 count"
    , "\nFiltered cells with >250 counts mapping to mouse exons for each sample"
    , "\n"))
ggsave(paste0(outGraphPfx, ".pdf"), width = 6, height = 8)
################################################################################