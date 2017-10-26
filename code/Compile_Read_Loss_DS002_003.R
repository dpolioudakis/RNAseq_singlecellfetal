






df <- read.table("../DS-003_DP/data/bam/human_mouse/SxaQSEQsVAP048L8/N706/unaligned_tagged_Cellular.bam_summary.txt",
  header = TRUE)

df <- read.table("../DS-003_DP/data/bam/human_mouse/SxaQSEQsVAP048L8/N706/unaligned_tagged_Molecular.bam_summary.txt",
  header = TRUE)


df$num_barcodes[1]/sum(df$num_barcodes)
sum(df$num_barcodes[-1])/sum(df$num_barcodes)





df <- read.table("../DS-002_DP/data/bam/human_mouse/SxaQSEQsXbp083L1/N701/out_cell_readcounts_readQ0.txt.gz", header = FALSE)





read.table("../DS-002_DP/data/QC/human_mouse/SxaQSEQsXbp083L1/Number_Of_Tagged_Exons.txt", header = TRUE)






inParent <- "../data/digital_gene_expression/human_mouse/SxaQSEQsXbp083L1"
dirsL <- list.files(inParent)

pctLDF <- list()
for (dir in dirsL) {
  dgeDF <- read.table(paste0(inParent, "/", dir
    , "/out_gene_exon_tagged_human.dge.summary.txt")
    , header = TRUE, stringsAsFactors=F)
  cgeDF <- read.table(paste0(inParent, "/", dir
    , "/out_gene_exon_tagged_human.counts.summary.txt")
    , header = TRUE, stringsAsFactors=F)
  # cgeDF <- cgeDF[1:400, ]
  df <- merge(dgeDF, cgeDF, by.x = "CELL_BARCODE", by.y = "CELL_BARCODE")
  df$PERCENT_REMAIN <- df$NUM_TRANSCRIPTS.x / df$NUM_TRANSCRIPTS.y * 100
  df$SAMPLE <- dir
  print(str(df))
  pctLDF[[dir]] <- df
}
ggDF <- do.call("rbind", pctLDF)
mnDF <- aggregate(PERCENT_REMAIN ~ SAMPLE, ggDF, mean)
mnDF$PERCENT_REMAIN <- round(mnDF$PERCENT_REMAIN, 2)


aggregate(NUM_TRANSCRIPTS.y  ~ SAMPLE, ggDF, sum)