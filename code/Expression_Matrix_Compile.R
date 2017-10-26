# Damon Polioudakis
# 2016-10-12
# Clustering of Drop-seq cells by digital gene expression

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())
sessionInfo()

### Digital gene expression - merged lanes

# Input directory paths
inSeqDir <- c("DS-002-011", "DS-003-004", "DS-005-006-007-008", "DS-008-011", "DS-009")
inDirs <- list.dirs(
  paste0("../", inSeqDir, "/data/digital_gene_expression/GRCh38_Gencode25"))
inDirs <- inDirs[c(grep("merged", inDirs, perl = TRUE)
  )]
inDirs <- inDirs[grep("N7[0-2][0-9]$", inDirs, perl = TRUE)]
# Don't use DS-002 N705
# inDirs <- inDirs[grep("SxaQSEQsXbp083L1/N705", inDirs, perl = TRUE, invert = TRUE)]

# Read in each lane as list of data frames
exLDF <- lapply(inDirs, function(inDirs) {
  if (grepl("DS-008-011|DS-009", inDirs)){
    read.table(paste0(inDirs, "/out_gene_exon_tagged_dge.txt.gz"), header = TRUE)
  }
  else {
    read.table(paste0(inDirs, "/out_gene_exon_tagged_dge_FtMm250.txt")) 
  }
})

## Combine samples into 1 data frame
for (i in 1:length(exLDF)) {
  df <- exLDF[[i]]
  if (i == 1) {
    exDF <- df}
  else {
    exDF <- merge(exDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
  }
}
str(exDF)
# 40k cells fails with 48G of RAM, works with 64G
exDF[is.na(exDF)] <- 0
row.names(exDF) <- exDF$GENE
exDF <- exDF[ ,-1]
print("Number of cells input:")
print(ncol(exDF))

## Metadata
# Nextera index
nexIdx <- c()
for (i in 1:length(inDirs)) {
  nexIdx <- c(nexIdx, rep(basename(inDirs)[i], ncol(exLDF[[i]])-1))
}
# Sequencing run
seq_run <- c()
for (i in 1:length(inDirs)) {
  seq_run <- c(seq_run, rep(basename(dirname(dirname(dirname(dirname(dirname(inDirs))))))[i], ncol(exLDF[[i]])-1))
}
# Metadata
metDF <- data.frame(CELL = colnames(exDF)
  , NEXTERA = nexIdx
  , SEQ_RUN = seq_run, stringsAsFactors = FALSE)
# SeqRun
metDF$SEQ_RUN[metDF$SEQ_RUN == "DS-005-006-007-008" & metDF$NEXTERA == "N710"] <- "DS-007-008"
metDF$SEQ_RUN[metDF$SEQ_RUN == "DS-005-006-007-008" & metDF$NEXTERA == "N711"] <- "DS-007-008"
metDF$SEQ_RUN[metDF$SEQ_RUN == "DS-005-006-007-008"] <- "DS-005-006"
# Brain
metDF$BRAIN <- 2
metDF$BRAIN[metDF$SEQ_RUN == "DS-005-006" | metDF$SEQ_RUN == "DS-007-008"] <- 3
metDF$BRAIN[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N712"] <- 4
metDF$BRAIN[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N714"] <- 4
metDF$BRAIN[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N715"] <- 5
metDF$BRAIN[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N716"] <- 5
metDF$BRAIN[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N722"] <- 4
metDF$BRAIN[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N723"] <- 4
metDF$BRAIN[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N724"] <- 5
metDF$BRAIN[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N726"] <- 5
# Region
metDF$REGION <- "GZ"
metDF$REGION[metDF$SEQ_RUN == "DS-007-008" & metDF$NEXTERA == "N711"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N705"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N702"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-002-011" & metDF$NEXTERA == "N701"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-002-011" & metDF$NEXTERA == "N704"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-003-004" & metDF$NEXTERA == "N706"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-003-004" & metDF$NEXTERA == "N709"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N712"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-008-011" & metDF$NEXTERA == "N715"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N722"] <- "CP"
metDF$REGION[metDF$SEQ_RUN == "DS-009" & metDF$NEXTERA == "N724"] <- "CP"
# Library prep
metDF$LIBRARY <- "Plath"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N702"] <- "Geschwind"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-005-006" & metDF$NEXTERA == "N701"] <- "Geschwind"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-007-008" & metDF$NEXTERA == "N710"] <- "Geschwind"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-007-008" & metDF$NEXTERA == "N711"] <- "Geschwind"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-008-011"] <- "Geschwind"
metDF$LIBRARY[metDF$SEQ_RUN == "DS-009"] <- "Geschwind"

save(exDF, metDF
  , file = "../analysis/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7-8-9-11.Rdata")
################################################################################

# ### Digital gene expression - each lane separately
# 
# # Input directory paths
# laneDirs <- list.files("../data/digital_gene_expression/GRCh37_75_assembly_NoERCC/")
# laneDirs <- laneDirs[grep("Sxa", laneDirs)]
# inDirs <- list.dirs(paste0("../data/digital_gene_expression/GRCh37_75_assembly_NoERCC/", laneDirs))
# inDirs <- inDirs[grep("N70[0-9]$", inDirs, perl = TRUE)]
# 
# # Read in each lane as list of dataframes
# exLDF <- lapply(inDirs, function(inDirs) {
#   read.table(paste0(inDirs, "/out_gene_exon_tagged_dge_FtMm250.txt"))
# })
# 
# ## Combine samples into 1 dataframe
# for (i in 1:length(exLDF)) {
#   df <- exLDF[[i]]
#   if (i == 1) {
#     exDF <- df}
#   else {
#     exDF <- merge(exDF, df, by.x = "GENE", by.y = "GENE", all = TRUE)  
#   }
# }
# str(exDF)
# exDF[is.na(exDF)] <- 0
# row.names(exDF) <- exDF$GENE
# exDF <- exDF[ ,-1]
# print("Number of cells input:")
# print(ncol(exDF))
# 
# ## Metadata
# # Samples
# samples <- c()
# for (i in 1:length(inDirs)) {
#   samples <- c(samples, rep(basename(inDirs)[i], ncol(exLDF[[i]])-1))
# }
# # Lanes
# lanes <- c()
# for (i in 1:length(inDirs)) {
#   lanes <- c(lanes, rep(basename(dirname(inDirs))[i], ncol(exLDF[[i]])-1))
# }
# # Metadata
# metDF <- data.frame(CELL = colnames(exDF), LANE = lanes, SAMPLE = samples)
# 
# save(exDF, metDF, file = "../analysis/Expression_Matrix_Compile_dge_FtMm250_LanesSeparate.Rdata")
# ################################################################################