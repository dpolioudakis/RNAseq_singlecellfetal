# Damon Polioudakis
# 2015-05-04
# Plot gene length (CDS length) vs fragment length (RNA STAR stats)
################################################################################

require(biomaRt)
require(ggplot2)

# load("../analysis/Expression_Matrix_Compile_dge_FtMm250_LanesSeparate.Rdata")

load("../analysis/Expression_Matrix_Compile_dge_FtMm250.Rdata")

ssDF <- read.csv("../analysis/tables/QC_RNAstar_Stats_HsMm.csv")

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



# Format for ggplot
df <- cbind(nGenesDF, rdepDF)
names(df) <- c("NUMBER_OF_GENES", "READ_DEPTH")
# Add sample info
df$SAMPLE <- paste0(metDF$LANE, "_", metDF$SAMPLE)

# ggplot
ggplot(df, aes(y = NUMBER_OF_GENES, x = READ_DEPTH, color = SAMPLE)) +
  geom_line() +
  coord_cartesian(xlim = c(0, 20000))+
  scale_colour_brewer(type = "qual", palette = 6) +
  ggtitle(paste0(
    "Number_Genes_Detected_Vs_Read_Depth.R"
    , "\n"
    , "\nGenes detected per cell versus number of reads mapping to exons after UMI collapse"
    , "\nSegRuns 5,6,7"
    , "\nMerged lanes"
    , "\nX-axis limit set to 20,000"
    , "\n"
  ))
ggsave("../analysis/graphs/Number_Genes_Detected_Vs_Read_Depth.pdf"
  , width = 8, height = 6)
################################################################################