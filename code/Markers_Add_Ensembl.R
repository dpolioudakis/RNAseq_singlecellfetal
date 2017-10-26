# Damon Polioudakis
# 2016-12-06
# Add ensembl IDs to marker gene lists output by Seurat
################################################################################

require(biomaRt)

# In marker gene lists
mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , header = TRUE)
################################################################################

### Functions

## Function: Convert HGNC symbol to Ensembl IDs
# Requires biomaRt package: require(biomaRt)
ConvertGeneSymbolToEnsembl <- function (genesList) {
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genesDF <- data.frame(genesList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Outputs data frame of attributes parameter from biomart
  bmDF <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id")
    , filters = "hgnc_symbol"
    , values = genesDF
    , mart = ensembl
  )
  bmDF
}
################################################################################
# Convert gene symbols to ensembl IDs
ensDF <- ConvertGeneSymbolToEnsembl(mkDF$gene)

# Add ensembl IDs to marker gene data frame
df <- merge(mkDF, ensDF, by.x = "gene", by.y ="hgnc_symbol", all.x = TRUE)
# A few gene symbols have multiple ensembl IDs

# Order by cluster
df <- df[order(df$cluster), ]

# Write out
write.table(df, "../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All_ens.txt"
  , quote = FALSE, row.names = FALSE, sep = "\t")
################################################################################
