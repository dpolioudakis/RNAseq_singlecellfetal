# Damon Polioudakis
# 2017-09-11
# Convert Molyneaux layer marker mouse entrez IDs to human homologs

require(biomaRt)

# Molyneaux entrez IDs
# etz <- c(
#   319582
#   ,68169
#   ,226438
#   ,59058
#   ,18993
#   ,18992
#   ,77220
#   ,27220
#   ,53870
#   ,50766
#   ,14219
#   ,58208
#   ,12971
#   ,56050
#   ,13047
#   ,13048
#   ,19049
#   ,50781
#   ,56419
#   ,207521
#   ,13603
#   ,14009
#   ,56458
#   ,114142
#   ,54713
#   ,140741
#   ,15902
#   ,16010
#   ,16019
#   ,80906
#   ,17311
#   ,17035
#   ,16826
#   ,16870
#   ,16873
#   ,66643
#   ,109593
#   ,16911
#   ,74762
#   ,17260
#   ,387184
#   ,17357
#   ,380684
#   ,18124
#   ,18227
#   ,67013
#   ,18423
#   ,18546
#   ,67784
#   ,170758
#   ,19699
#   ,67792
#   ,225998
#   ,20194
#   ,20349
#   ,18991
#   ,76965
#   ,20678
#   ,266615
#   ,60510
#   ,21375
#   ,21577
#   ,21885
#   ,21887
#   ,21888
#   ,210801
# )

# Molyneaux markers ordered by figure 3
etz <- c(16873, 68169, 19699, 226438, 67792, 27220, 21885, 21887, 74762, 387184, 80906
  , 225998, 56050, 16870, 210801, 140741, 17260, 207521, 13047, 13048, 17311, 266615
  , 59058, 67784, 18991, 17357, 18993, 18992, 14009, 319582, 77220, 380684, 53870
  , 56458, 13603, 66643, 60510, 20194, 67013, 16826, 50766, 18546, 170758, 56419
  , 58208, 12971, 18423, 54713, 16010, 21577, 20678, 50781, 21888, 20349, 18124
  , 17035, 114142, 19049, 16019, 15902, 76965, 21375, 109593, 14219, 18227)
etz[etz == 266615] <- 137970
etz[etz == 21577] <- 28639

# bioMart manual:
# http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
# Attribute: values to retrieve
# Filters: input query type
# Values: input query

# Look up mouse ensembl ID and MGI symbol
genesDF <- data.frame(etz)
ensemblMm = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
ensemblMm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensemblMm <- useDataset("mmusculus_gene_ensembl", mart=ensemblMm)
# Outputs data frame of attributes parameter from biomart
bm1DF <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene", "mgi_symbol")
  , filters = "entrezgene"
  , values = genesDF
  , mart = ensemblMm
)

# Look up human homolog
bm2DF <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene")
  , filters = "ensembl_gene_id"
  , values = bm1DF$ensembl_gene_id
  , mart = ensemblMm
)

# Look up human gene symbol
ensemblHs <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensemblHs <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensemblHs <- useDataset("hsapiens_gene_ensembl", mart=ensemblHs)
bm3DF <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol")
  , filters = "ensembl_gene_id"
  , values = bm2DF$hsapiens_homolog_ensembl_gene
  , mart = ensemblHs
)

# Combine biomart query data frames
bmDF <- merge(bm1DF, bm2DF, by = "ensembl_gene_id")
bmDF <- merge(bmDF, bm3DF, by.x = "hsapiens_homolog_ensembl_gene", by.y = "ensembl_gene_id")

# Order table by Molyneaux figure 3
bmDF <- bmDF[match(etz, bmDF$entrezgene), ]
# Remove NAs from IDs not found in biomart
bmDF <- bmDF[! is.na(bmDF$entrezgene), ]

# Save as csv
write.csv(bmDF, file = "../source/Molyneaux_LayerMarkers_Format.csv"
  , quote = FALSE, row.names = FALSE)