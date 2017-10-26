# Damon Polioudakis
# 2017-04-13
# Frequency of gene classes for each drop-seq sample as well as length and GC
################################################################################

rm(list=ls())
sessionInfo()
set.seed(27)

require(biomaRt)
require(ggplot2)
require(gProfileR)
require(reshape2)

# exDF, metDF
load("../analysis/Expression_Matrix_Compile_dge_FtMm250_DS-2-3-4-5-6-7.Rdata")

## Variables
graphCodeTitle <- "Classes_Of_Genes.R"
graphTitle <- "Drop-seq SeqRuns 2,3,4,5,6"
outGraphPfx <- "../analysis/graphs/Classes_Of_Genes_DS-2-3-4-5-6-7"

## Output Directories
outGraphDir <- dirname(outGraphPfx)
dir.create(outGraphDir, recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size = 18)))
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

# Function to load digital expression matrix and sum genes for each cell with
# >= 1 count
# Input is path to digital expession matrix
# Outputs dataframe with cell barcodes as row names and number genes detected
# column
Calc_Genes_Detected <- function (inDeDF) {
  deDF <- read.table(inDeDF, header = TRUE, stringsAsFactors = FALSE
    , row.names = 1)
  # Remove gene name column
  deDF <- deDF[ ,-1]
  data.frame(colSums(deDF >= 1))
}
################################################################################

## Look up gene info on biomart

# Convert ensembl ID to gene symbol and subset to ensembl IDs found in biomaRt
df <- QueryBiomaRt(row.names(exDF), "hgnc_symbol"
  , c("ensembl_gene_id", "hgnc_symbol"))
geneInfoDF <- df
# Add gene information
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "gene_biotype", "percentage_gc_content"))
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "cds_length"))
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
str(geneInfoDF)
head(geneInfoDF)
tail(geneInfoDF)

metDFL <- split(metDF, list(metDF$SEQ_RUN, metDF$NEXTERA))
# Remove empty dataframes
metDFL <- metDFL[sapply(metDFL, function(x) dim(x)[1]) > 0]

## Frequency of gene biotypes for genes detected >= 1 count

pctDFL <- lapply(names(metDFL), function(sampleID) {
  metDF <- metDFL[[sampleID]]
  ssExDF <- exDF[ ,colnames(exDF) %in% metDF$CELL]
  detDF <- data.frame(NUMBER_DETECTED = rowSums(ssExDF >= 1))
  detDF <- merge(detDF, geneInfoDF, by.x = "row.names", by.y = "hgnc_symbol")
  pctDF <- aggregate(detDF$NUMBER_DETECTED, list(detDF$gene_biotype), sum)
  pctDF$x <- (pctDF$x / sum(pctDF$x)) * 100
  pctDF$SAMPLE <- sampleID
  return(pctDF)
})
pctDF <- do.call("rbind", pctDFL)

names(pctDF)[1:3] <- c("BIOTYPE", "PERCENT", "SAMPLE")

# ggplot
ggplot(pctDF, aes(x = BIOTYPE, y = PERCENT)) +
  geom_bar(aes(fill = SAMPLE), position = "dodge", stat="identity") +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0(graphCodeTitle
    , "\n\n", graphTitle
    , "\nFrequency of gene biotypes detected >=1 count"
    , "\n"))
ggsave(paste0(outGraphPfx, "_Biotype_1.pdf"))

## Frequency of gene biotypes for genes detected >= 3 counts

pctDFL <- lapply(names(metDFL), function(sampleID) {
  metDF <- metDFL[[sampleID]]
  ssExDF <- exDF[ ,colnames(exDF) %in% metDF$CELL]
  detDF <- data.frame(NUMBER_DETECTED = rowSums(ssExDF >= 3))
  detDF <- merge(detDF, geneInfoDF, by.x = "row.names", by.y = "hgnc_symbol")
  pctDF <- aggregate(detDF$NUMBER_DETECTED, list(detDF$gene_biotype), sum)
  pctDF$x <- (pctDF$x / sum(pctDF$x)) * 100
  pctDF$SAMPLE <- sampleID
  return(pctDF)
})
pctDF <- do.call("rbind", pctDFL)

names(pctDF)[1:3] <- c("BIOTYPE", "PERCENT", "SAMPLE")

# ggplot
ggplot(pctDF, aes(x = BIOTYPE, y = PERCENT)) +
  geom_bar(aes(fill = SAMPLE), position = "dodge", stat="identity") +
  coord_flip() +
  theme_bw() +
  ggtitle(paste0(graphCodeTitle
    , "\n\n", graphTitle
    , "\nFrequency of gene biotypes detected >=3 count"
    , "\n"))
ggsave(paste0(outGraphPfx, "_Biotype_3.pdf"))

## Gene length

lengthDFL <- lapply(names(metDFL), function(sampleID) {
  metDF <- metDFL[[sampleID]]
  ssExDF <- exDF[ ,colnames(exDF) %in% metDF$CELL]
  df <- merge(geneInfoDF, ssExDF, by.x = "hgnc_symbol", by.y = "row.names")
  df <- na.omit(df)
  # Length
  lengthDF <- data.frame(apply(df[6:ncol(df)], 2, function(col) {
    sum((col * df$cds_length)) / sum(col)
  }))
  lengthDF$SAMPLE <- sampleID
  return(lengthDF)
})
lengthDF <- do.call("rbind", lengthDFL)

names(lengthDF) <- c("LENGTH", "SAMPLE")

# Annotate library
lengthDF$LIBRARY <- "Plath"
lengthDF$LIBRARY[grep("DS-005-006-007.N70[1,2]", lengthDF$SAMPLE, perl = TRUE)] <- "Geschwind"
lengthDF$LIBRARY[grep("DS-005-006-007.N71[0,1]", lengthDF$SAMPLE, perl = TRUE)] <- "Geschwind"
# Annotate Region
lengthDF$REGION <- metDF$REGION[match(row.names(lengthDF), as.character(metDF$CELL))]
# Label for graph
lengthDF$LABEL <- paste(lengthDF$SAMPLE, lengthDF$LIBRARY, lengthDF$REGION)
# Mean
mnDF <- aggregate(LENGTH ~ LABEL, lengthDF, mean)
mnDF$LENGTH <- round(mnDF$LENGTH, 1)

# ggplot
ggplot(lengthDF, aes(x = LABEL, y = LENGTH)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 4000)) +
  geom_text(data = mnDF, aes(label = LENGTH, x = LABEL, y = 3000)) +
  # geom_bar(aes(fill = SAMPLE), position = "dodge", stat="identity") +
  # coord_flip() +
  # theme_bw() +
  ggtitle(paste0(graphCodeTitle
    , "\n\n", graphTitle
    , "\nCDS length for each cell"
    , "\nText indicates mean"
    , "\nY-axis limit set at 4000"
    , "\n"))
ggsave(paste0(outGraphPfx, "_Length.pdf"), width = 12, height = 8)

## Gene GC

gcDFL <- lapply(names(metDFL), function(sampleID) {
  metDF <- metDFL[[sampleID]]
  ssExDF <- exDF[ ,colnames(exDF) %in% metDF$CELL]
  df <- merge(geneInfoDF, ssExDF, by.x = "hgnc_symbol", by.y = "row.names")
  df <- na.omit(df)
  # Length
  lengthDF <- data.frame(apply(df[6:ncol(df)], 2, function(col) {
    sum((col * df$percentage_gc_content)) / sum(col)
  }))
  lengthDF$SAMPLE <- sampleID
  return(lengthDF)
})
gcDF <- do.call("rbind", gcDFL)

names(gcDF) <- c("GC", "SAMPLE")

# Annotate library
gcDF$LIBRARY <- "Plath"
gcDF$LIBRARY[grep("DS-005-006-007.N70[1,2]", gcDF$SAMPLE, perl = TRUE)] <- "Geschwind"
gcDF$LIBRARY[grep("DS-005-006-007.N71[0,1]", gcDF$SAMPLE, perl = TRUE)] <- "Geschwind"
# Annotate Region
gcDF$REGION <- metDF$REGION[match(row.names(gcDF), as.character(metDF$CELL))]
# Label for graph
gcDF$LABEL <- paste(gcDF$SAMPLE, gcDF$LIBRARY, gcDF$REGION)
# Mean
mnDF <- aggregate(GC ~ LABEL, gcDF, mean)
mnDF$GC <- round(mnDF$GC, 1)

# ggplot
ggplot(gcDF, aes(x = LABEL, y = GC)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = mnDF, aes(label = GC, x = LABEL, y = 50)) +
# geom_bar(aes(fill = SAMPLE), position = "dodge", stat="identity") +
# coord_flip() +
# theme_bw() +
  ggtitle(paste0(graphCodeTitle
    , "\n\n", graphTitle
    , "\nGC length for each cell"
    , "\nText indicates mean"
    , "\n"))
ggsave(paste0(outGraphPfx, "_GC.pdf"), width = 12, height = 8)