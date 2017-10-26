# Damon Polioudakis
# 2017-04-05
# scRNA-seq marker enrichment in Neel's ASD modules

rm(list=ls())
options(stringsAsFactors = FALSE)

library(biomaRt)
library(reshape2)
library(ggplot2);

# Neel modules from Mike direct message on slack
load("../reference/Neel_neuron_mods.Rdata")

# Drop-seq markers
scMksDF <- read.table("../analysis/tables/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , header = TRUE)

## Cell types in fetal brain (Pollen et al S4) Genes exclusive to each Cell Type from Luis metaMat.R
load(file="../reference/PollenCellTypesGenes-Unique.R")

outGraph <- "../analysis/graphs/ASD_Module_Enrichment_DS002003_"

## Function: Fischers
OR <- function(q,k,m,t) {
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  output <- c(OR,pval,upCI,downCI);   names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

## Function: Count overlaps and run the Fischers analysis
ORA <- function(testpath,refpath,testbackground,refbackground) {
  testpath = testpath[testpath %in% testbackground]
  refpath = refpath[refpath %in% refbackground]
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}

# Background gene list
## First, we get all the genes from Gencode v19
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
geneAnno1 <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)

# ## Calculated by hand to confirm function output
# # Drop-seq
# for(i in 1:length(testset)) {
#   cls <- unique(scMksDF$cluster)
#   print(names(testset)[[i]])
#   for (j in cls) {
#     A <- testset[[i]]$Symbol
#     B <- scMksDF$gene[scMksDF$cluster == j]
#     n <- 57241
#     ctM <- matrix(c(n - length(union(A,B)), length(setdiff(A,B)), length(setdiff(B,A)), length(intersect(A,B))), nrow=2)
#     print(fisher.test(ctM, alternative="greater")$p.value)
#   }
# }
# 
# 
# # Pollen
# for(i in 1:length(testset)) {
#   print(names(testset)[[i]])
#   for (j in 1:length(CellTypesGenes)) {
#     print(names(CellTypesGenes)[[j]])
#     A <- testset[[i]]$ENSG
#     B <- CellTypesGenes[[j]]
#     n <- 57241
#     ctM <- matrix(c(n - length(union(A,B)), length(setdiff(A,B)), length(setdiff(B,A)), length(intersect(A,B))), nrow=2)
#     print(fisher.test(ctM, alternative="greater")$p.value)
#   }
# }

# Pollen
plOrM <- matrix(nrow = length(testset), ncol = length(CellTypesGenes)
  , dimnames = list(names(testset), names(CellTypesGenes)))
for(i in 1:length(testset)) {
  print(names(testset)[[i]])
  for (j in 1:length(CellTypesGenes)) {
    plOrM[i, j] <- as.numeric(ORA(testset[[i]]$ENSG, CellTypesGenes[[j]]
      , geneAnno1$ensembl_gene_id, geneAnno1$ensembl_gene_id)[2])
  }
}
plOrM <- -log10(plOrM)
plOrM <- round(plOrM, 2)

# Drop-seq
dsOrM <- matrix(nrow = length(testset), ncol = length(unique(scMksDF$cluster))
  , dimnames = list(names(testset), c(
    "0 Excitatory Upper Layer Neuron"
    , "1 Excitatory Neuron"
    , "2 Excitatory Upper Layer Neuron"
    , "3 Excitatory Deep Layer Neuron"
    , "4 Intermediate Progenitors"
    , "5 Interneuron"
    , "6 Mitotic Progenitors"
    , "7 oRG"
    , "8 Oligodendrocyte Precursor"
    , "9 Endothelial")))
for(i in 1:length(testset)) {
  print(names(testset)[[i]])
  for (j in 1:length(unique(scMksDF$cluster))) {
    dsOrM[i, j] <- as.numeric(ORA(testset[[i]]$Symbol, scMksDF$gene[scMksDF$cluster == j]
      , geneAnno1$hgnc_symbol, geneAnno1$hgnc_symbol)[2])
  }
}
dsOrM <- -log10(dsOrM)
dsOrM <- round(dsOrM, 2)

# Drop-seq
ggDF <- melt(dsOrM)
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_continuous(low = "white", high = "#e31a1c") +
  geom_text(data = ggDF, aes(x = Var2, y = Var1, label = value)
    , vjust = -0.35, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(legend.position = "none") + 
  ylab("Neel ASD modules") +
  xlab("Geschwind Drop-seq cluster marker lists") +
  ggtitle(paste0("Geschwind Drop-seq Cell-type Enrichment"
    , "\n-log10(p-value)"
    , "\n"))
ggsave(paste0(outGraph, "Heatmap_DropSeq.png"), height = 6, width = 6)

# Pollen
ggDF <- melt(plOrM)
ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_continuous(low = "white", high = "#e31a1c") +
  geom_text(data = ggDF, aes(x = Var2, y = Var1, label = value)
    , vjust = -0.35, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(legend.position = "none") + 
  ylab("Neel ASD modules") +
  xlab("Pollen cluster marker lists") +
  ggtitle(paste0("Pollen Cell-type Enrichment"
    , "\n-log10(p-value)"
    , "\n"))
ggsave(paste0(outGraph, "Heatmap_Pollen.png"), height = 6, width = 6)