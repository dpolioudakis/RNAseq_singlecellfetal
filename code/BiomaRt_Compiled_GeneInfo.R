# Damon Polioudakis
# 2017-05-05
# Make table of biomaRt gene IDs and info
################################################################################

rm(list=ls())

require(biomaRt)
################################################################################

### Functions

## Function: Query biomaRt
QueryBiomaRt <- function (genesList, filters, attributes, host) {
  genesDF <- data.frame(genesList)
  # Attribute: values to retrieve - c("hgnc_symbol", "ensembl_gene_id")
  # Filters: input query type
  # Values: input query
  ensembl= useMart(biomart = "ENSEMBL_MART_ENSEMBL"
    , host = host, path = "/biomart/martservice"
    , dataset = "hsapiens_gene_ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
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

### Pull gene info from biomaRt

## feb2014.archive.ensembl.org
# Ensembl 75

# Get all ensembl IDs
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL"
  , host = "feb2014.archive.ensembl.org"
  , path = "/biomart/martservice"
  , dataset = "hsapiens_gene_ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
geneInfoDF <- getBM(attributes = c("ensembl_gene_id")
  , values = "*"
  , mart = ensembl)

# Add gene information
# host = "feb2014.archive.ensembl.org"
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "hgnc_symbol", "entrezgene")
  , host = "feb2014.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "gene_biotype", "percentage_gc_content")
  , host = "feb2014.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "cds_length")
  , host = "feb2014.archive.ensembl.org")
# Take longest CDS
df <- aggregate(cds_length ~ ensembl_gene_id, data = df, max)
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id", all = TRUE)
# Chromosome
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "chromosome_name")
  , host = "feb2014.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id", all.x = TRUE)
str(geneInfoDF)
head(geneInfoDF)
tail(geneInfoDF)

# Write csv
write.csv(geneInfoDF
  , file = "../source/BiomaRt_Compile_GeneInfo_GRCh37_Ensembl75.csv"
  , quote = FALSE)

## dec2016.archive.ensembl.org
# Ensembl 87
# refFlat from Chris is Gencode 25 which corresponds to Ensembl 85, 86, 87
# according to Gencode website

# Get all ensembl IDs
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL"
  , host = "dec2016.archive.ensembl.org"
  , path = "/biomart/martservice"
  , dataset = "hsapiens_gene_ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
geneInfoDF <- getBM(attributes = c("ensembl_gene_id")
  , values = "*"
  , mart = ensembl)

# Add gene information
# host = "dec2016.archive.ensembl.org"
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "hgnc_symbol", "entrezgene")
  , host = "dec2016.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "gene_biotype", "percentage_gc_content")
  , host = "dec2016.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id")
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "cds_length")
  , host = "dec2016.archive.ensembl.org")
# Take longest CDS
df <- aggregate(cds_length ~ ensembl_gene_id, data = df, max)
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id", all = TRUE)
# Chromosome
df <- QueryBiomaRt(geneInfoDF$ensembl_gene_id, "ensembl_gene_id"
  , c("ensembl_gene_id", "chromosome_name")
  , host = "dec2016.archive.ensembl.org")
geneInfoDF <- merge(geneInfoDF, df, by = "ensembl_gene_id", all.x = TRUE)
str(geneInfoDF)
head(geneInfoDF)
tail(geneInfoDF)

# Write csv
write.csv(geneInfoDF
  , file = "../source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"
  , quote = FALSE)
################################################################################