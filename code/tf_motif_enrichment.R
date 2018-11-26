# Damon Polioudakis
# 2018-11-20
# TF motif enrichment of cell type mapped REs

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(motifmatchr)
require(GenomicRanges)
require(SummarizedExperiment)
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)
require(chromVAR)
require(tidyverse)
source("Function_Library.R")

## inputs
re_tb <- read_csv(
  "../analysis/tables/ATAC_REs/ATAC_REs_ATAConly_CellTypeEnrichment.csv")

## other variables
script_name <- "tf_motif_enrichment.R"
date <- format(Sys.Date(), "%Y%m%d")
cluster_annot_tb <- tribble(
    ~cluster_number, ~cluster_annot
    , "9",  "vRG"
    , "7",  "oRG"
    , "8",  "PgS"
    , "10", "PgG2M"
    , "2",  "IP"
    , "0",  "ExN"
    , "1",  "ExM"
    , "4",  "ExCal"
    , "3",  "ExDp1"
    , "13", "ExDp2"
    # , "17", "SP"
    , "5",  "InSST"
    , "6",  "InCALB2"
    , "11", "OPC"
    , "12", "End"
    , "14", "Per"
    , "15", "Mic"
    , "16", "NA"
  )

## outputs
out_table <- paste0(
  "../analysis/tables/tf_motif_enrichment/", date, "/tf_motif_enrichment_")
out_scratch <- paste0(
  "/u/flashscratch/d/dpolioud/tf_motif_enrichment/", date
  , "/tf_motif_enrichment_")

# make output directories
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_scratch), recursive = TRUE)
################################################################################

### motif enrichment

# cell type mapped REs
re_gr <- with(re_tb
  , GRanges(
    seqnames = c(paste0("chr", enhancerchr), paste0("chr", promoterchr))
      , ranges = IRanges(
        , start = c(enhancerstart, promoterstart)
        , end = c(enhancerend, promoterend)
    # strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  )))

# background
# sample sequences randomly from genome
# lengths of sequences are mean enhancer length
# mean enhancer length
seq_length <- as.integer(with(re_tb, mean(enhancerend - enhancerstart)))
# sample
background_seqs <-
  # sample 1000 from each chromosome
  tibble(chr = rep(1:22, 1000) %>% sort) %>%
  mutate(seq_start = pmap_chr(., .f = function(chr){
    sample(1:seqlengths(Hsapiens)[chr] - seq_length, 1)
  }) %>% as.integer) %>%
  mutate(seq_end = seq_start + seq_length) %>%
  mutate(chr = paste0("chr", chr))
# genomic ranges
hs_genome_gr <- with(background_seqs
  , GRanges(
    seqnames = chr,
    ranges = IRanges(
      start = seq_start
      , end = seq_end
  )))

## motif matches
motifs <- chromVAR::getJasparMotifs()
# get motif matches for motifs in enhancers
motif_ix <- matchMotifs(motifs, re_gr, genome = "hg19")
# background motif matches
motif_hs_genome_ix <- matchMotifs(motifs, hs_genome_gr, genome = "hg19")
# save to scratch
save(motif_ix, motif_hs_genome_ix
  , file = paste0(out_scratch, "motif_matches.rdata"))

## calculate enrichment with fishers test
tf_cluster_tb <- expand.grid(
  colnames(motifMatches(motif_ix))
  , re_tb$Cluster %>% unique
  ) %>% as_tibble %>% mutate(Var1 = as.character(Var1))
tf_enrich_tb <- pmap(tf_cluster_tb, .f = function(Var1, Var2){
  tf <- Var1
  cluster <- Var2
  cluster_idx <- re_tb$Cluster == cluster
  tf_idx <- colnames(motifMatches(motif_ix)) == tf
  contg_tbl_m <- rbind(
    # number of cell type enriched elements with motif
    c(motifMatches(motif_ix)[cluster_idx, tf_idx] %>% sum
    # number of cell type enriched elements without motif
    , (! motifMatches(motif_ix)[cluster_idx, tf_idx]) %>% sum
    )
    # background number of cell type enriched elements with motif
    , c(motifMatches(motif_hs_genome_ix)[ , tf_idx] %>% sum
    # background  number of cell type enriched elements without motif
    , (! motifMatches(motif_hs_genome_ix)[ , tf_idx]) %>% sum
    )
  )
  ft <- fisher.test(contg_tbl_m)
  tf_enrich_tb <- tibble(
    cluster_number = as.character(cluster), tf = tf, pvalue = ft$p.value
    , odds_ratio = ft$estimate)
  return(tf_enrich_tb)
}) %>% bind_rows

# add cluster annotations
tf_enrich_tb <- tf_enrich_tb %>%
  left_join(cluster_annot_tb) %>%
  # order
  mutate(cluster_annot = factor(cluster_annot
    , levels = cluster_annot_tb$cluster_annot)) %>%
  arrange(cluster_annot, pvalue)

write.csv(tf_enrich_tb, file = paste0(out_table, "vs_hs_genome.csv")
  , quote = FALSE, row.names = FALSE)
################################################################################
