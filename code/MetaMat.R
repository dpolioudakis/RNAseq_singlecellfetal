
# Damon Polioudakis
# 2016-10-25
# Derived from Luis' script:
#   /geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/Enrichment/Metamat-CorticalLayers.R
################################################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

require(reshape2)
require(ggplot2)
# require(xlsx)
library(biomaRt) ## First, construct a matrix containing all possible protein coding genes. We will use human gene symbols. Although other identifiers, such as ENSG IDs are preferred, many studies only report human gene symbols, so we will use this as the common identifier.
library(WGCNA) ## For plotting heatmaps
source("Function_Library.R")
################################################################################

### Functions

Subset_OR_And_Log2_Odds_Ratio_Matrix <- function(
  target_list_cols, interpretation_list_cols, ORmat){
  print("Subset_OR_And_Log2_Odds_Ratio_Matrix")
  subORmat <- ORmat[target_list_cols, interpretation_list_cols]
  # Tranform to log2 odds ratio
  dispMat <- log2(exp(subORmat))
  rownames(dispMat) <- rownames(subORmat)
  return(dispMat)
}

Subset_OR_Pval_And_FDR_Correct_Matrix <- function(
  target_list_cols, interpretation_list_cols, Pmat){
  print("Subset_OR_Pval_And_FDR_Correct_Matrix")
  subPmat <- Pmat[target_list_cols, interpretation_list_cols]
  # FDR correction
  # Correcting for all columns and rows
  subPmat <- p.adjust(subPmat, method="BH")
  return(subPmat)
}

Format_OR_Pval_Matrices_For_Barplot <- function(
  dispMat, subPmat, cluster_order){
  print("Format_OR_Pval_Matrices_For_Barplot")
  # Order clusters
  new_cluster_order <- cluster_order
  gg_DF <- melt(dispMat)
  gg_DF$Cluster <- gsub("DS2-11 ", "", gg_DF$Var2)
  gg_DF$Cluster <- factor(gg_DF$Cluster, levels = c(0:15))
  gg_DF$Cluster2 <- factor(gg_DF$Cluster, levels = new_cluster_order)
  # Add p-values
  gg_DF$Pvalue_FDR <- subPmat
  pvals_DF <- melt(Pmat[target_list_cols, interpretation_list_cols, drop=FALSE])
  pvals_DF$Index <- paste0(pvals_DF$Var1, "_", pvals_DF$Var2)
  gg_DF$Index <- paste0(gg_DF$Var1, "_", gg_DF$Var2)
  idx <- match(gg_DF$Index, pvals_DF$Index)
  gg_DF$Pvalue <- pvals_DF$value[idx]
  # Remove clusters not in list
  gg_DF <- gg_DF[gg_DF$Cluster %in% cluster_order, ]
  # Label for coloring bars
  gg_DF <- gg_DF[order(gg_DF$Cluster2), ]
  class_cluster_idx <- c(
    "Radial glia" = 7
    , "Radial glia" = 9
    , "Cycling progenitor" = 8
    , "Cycling progenitor" = 10
    , "Intermediate progenitor" = 2
    , "Excitatory Neuron" = 0
    , "Excitatory Neuron" = 1
    , "Excitatory Neuron" = 4
    , "Excitatory Neuron" = 3
    , "Excitatory Neuron" = 13
    , "Interneuron" = 5
    , "Interneuron" = 6
    , "Oligodendrocyte precursor" = 11
    , "Endothelial" = 12
    , "Pericyte" = 14
    , "Microglia" = 15
  )
  gg_DF$Class <- names(class_cluster_idx)[match(gg_DF$Cluster, class_cluster_idx)]
  gg_DF$Class <- factor(gg_DF$Class
    , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
      , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
      , "Endothelial", "Pericyte", "Microglia")
  )
  # Remove NAs
  gg_DF <- gg_DF[! is.na(gg_DF$Cluster), ]
}
################################################################################

## First, we get all the genes from Gencode v19
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
geneAnno1 <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)
geneAnno1 <- geneAnno1[geneAnno1[,"gene_biotype"]=="protein_coding",]; geneAnno1 <- geneAnno1[!duplicated(geneAnno1[,"hgnc_symbol"]) & geneAnno1[,"hgnc_symbol"]!="",] ## 19163 with hgnc symbols after removing duplicates. We use gene symbols here because the RDNV data comes in gene symbol format. Otherwise, I prefer to use ENSG IDs.

## Construct a matrix to store the gene lists - this starts as a matrix with gene information and list membership information, and new lists are appended. If the new list contains the gene, it gets a 1, if not 0, if the gene is not even in the study (not in the background) we mark an NA
metaMat <- matrix(NA,ncol=1,nrow=nrow(geneAnno1)) ## We can add columns as necessary
rownames(metaMat) <- geneAnno1[,"hgnc_symbol"]
ENSGIDs <- geneAnno1[,"ensembl_gene_id"]
geneLength <- geneAnno1[,"end_position"]-geneAnno1[,"start_position"]
names(geneLength) <- geneAnno1[,"hgnc_symbol"] ## For use later in making the covariate matrix
metaX = geneAnno1[match(rownames(metaMat), geneAnno1$hgnc_symbol),"chromosome_name"]

################################################
## Single-cell RNA-seq Datasets
################################################

## Cell types in fetal brain (Pollen et al S3)
load(file="../source/Gene_Lists/PollenS3GenesENSG.Rda");
for(i in 1:length(PollenS3GenesENSG)){
  thislist = PollenS3GenesENSG[[i]]
  gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
  gename = gename[gename!=""]
  matchlist = match(rownames(metaMat), gename)
  matchlist = ifelse(is.na(matchlist),0,1)
  metaMat = cbind(metaMat, matchlist)
  colnames(metaMat)[dim(metaMat)[2]] = names(PollenS3GenesENSG)[i]
}

## Cell types in fetal brain (Pollen et al S4)
load(file="../source/Gene_Lists/PollenCellTypesGenes.R");
for(i in 1:length(CellTypesGenes)){
  thislist = CellTypesGenes[[i]]
  gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
  gename = gename[gename!=""]
  matchlist = match(rownames(metaMat), gename)
  matchlist = ifelse(is.na(matchlist),0,1)
  metaMat = cbind(metaMat, matchlist)
  colnames(metaMat)[dim(metaMat)[2]] = names(CellTypesGenes)[i]
}

## Cell types in fetal brain (Pollen et al S4) Genes exclusive to each Cell Type
load(file="../source/Gene_Lists/PollenCellTypesGenes-Unique.R")
names(CellTypesGenes) = c("IPC-Uniq", "Neuron-Uniq","Interneuron-Uniq", "vRG-Uniq", "oRG-Uniq");
for(i in 1:length(CellTypesGenes)){
  thislist = CellTypesGenes[[i]]
  gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
  gename = gename[gename!=""]
  matchlist = match(rownames(metaMat), gename)
  matchlist = ifelse(is.na(matchlist),0,1)
  metaMat = cbind(metaMat, matchlist)
  colnames(metaMat)[dim(metaMat)[2]] = names(CellTypesGenes)[i]
}

## Nowakowski 2017
Classify_Cells_By_Type <- function(
  cluster_DE_DF, class_cluster_idx, cluster_col_name){
  print("Classify_Cells_By_Type")
  # Example class_cluster_idx input
  # $`Glia or support cells`
  # [1] "End"     "Per"     "Ast"     "Ast_Cer" "Oli"     "OPC"     "OPC_Cer"
  # [8] "Mic"
  # $`NA`
  # [1] "Gran"
  class_cluster_idx <- melt(class_cluster_idx)
  cluster_DE_DF$Class <- class_cluster_idx$L1[
    match(cluster_DE_DF[[cluster_col_name]], class_cluster_idx$value)]
  return(cluster_DE_DF)
}
# Nowakowski cluster DE table
nowakowski_DE_DF <- read.csv(
  "../nowakowski_2017/Nowakowski_Table_S5_Clustermarkers.csv", header = TRUE
)
# Nowakowski DE not very cell type specifc, set DE filter higher
nowakowski_DE_DF <- nowakowski_DE_DF[nowakowski_DE_DF$avg_diff > 0.75, ]
# Classify cells by type
class_cluster_idx <- list(
  "Astrocyte" = "Astrocyte"
  , "Choroid" = "Choroid"
  , "Mitotic progenitor" = c("IPC-div1", "RG-div1", "RG-div2", "IPC-div2")
  , "RG" = "RG-early"
  , "Endothelial" = "Endothelial"
  , "IP" = c("IPC-nEN1", "IPC-nEN2", "IPC-nEN3")
  , "Microglia" = "Microglia"
  , "OPC" = "OPC"
  , "oRG" = "oRG"
  , "tRG" = "tRG"
  , "vRG" = "vRG"
  , "MGE IPC" = c("MGE-IPC1", "MGE-IPC2")
  , "MGE mitotic progentior" = "MGE-div"
  , "MGE RG" = c("MGE-RG1", "MGE-RG2")
  , "MGE newborn neurons" = c("nIN1", "nIN2", "nIN3", "nIN4", "nIN5")
  , "Mural" = "Mural"
  , "NA" = c("Glyc", "U1", "U2", "U3", "U4")
  , "Excitatory" = c("EN-PFC2", "EN-PFC3", "EN-V1-2", "EN-PFC1", "EN-V1-1", "EN-V1-3")
  , "Excitatory Early" = c("nEN-early2", "nEN-early1")
  , "Excitatory Late" = "nEN-late"
  , "Interneuron" = c("IN-CTX-CGE", "IN-CTX-CGE2", "IN-CTX-MGE2", "IN-CTX-MGE1")
  , "Striatal neurons" = "IN-STR"
)
nowakowski_DE_DF <- Classify_Cells_By_Type(
  cluster_DE_DF = nowakowski_DE_DF
  , class_cluster_idx = class_cluster_idx
  , cluster_col_name = "cluster"
)
Add_Cluster_DE_Data_To_MetaMat <- function(
  cluster_DE_DF, gene_col_name, cluster_col_name, metaMat){
  cluster_DE_DFL <- split(cluster_DE_DF[[gene_col_name]]
    , f = cluster_DE_DF[[cluster_col_name]])
  for(i in 1:length(cluster_DE_DFL)){
    thislist = cluster_DE_DFL[[i]]
    gename = geneAnno1[geneAnno1$hgnc_symbol %in% thislist, 2]
    gename = gename[gename!=""]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = names(cluster_DE_DFL)[i]
  }
  return(metaMat)
}
metaMat <- Add_Cluster_DE_Data_To_MetaMat(
  cluster_DE_DF = nowakowski_DE_DF
  , gene_col_name = "gene"
  , cluster_col_name = "Class"
  , metaMat = metaMat
)

################################################
## Disease gene datasets
################################################

## Sanders et al., 2015, PIMD 26402605
listMat <- read.csv("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/Sanders_2015_TADA.csv")
listMat <- listMat[,1:21]

thisgene = unique(listMat[listMat$tadaFdrAscSscExomeSscAgpSmallDel<0.1, "RefSeqGeneName"])
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "Sanders_TADA0.1_exomedel"

thisgene = unique(listMat[listMat$tadaFdrAscSscExome<0.1, "RefSeqGeneName"])
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "Sanders_TADA0.1_exome"

## Iossifov ASD - load for exomeLength
listMat <- read.csv("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/AllDeNovoByGene.csv") ## From the supplmental tables of Iossifov et al., 2014 - note that some gene symbols have the excel conversion error in this data, but it will get removed when we overlap things
listMat <- listMat[!duplicated(listMat[,1]),]
rownames(listMat) <- listMat[,"gene"]
exomeLength <- listMat[,"codingLenInTarget"] ## For later use in making the covariate matrix
names(exomeLength) <- rownames(listMat)

## iHART
listMat <- read.csv("../source/Gene_Lists/ASD.risk-genes.ForDamon.SingleCellExp_2018-04-18.csv")
listMat$HGNC.gene.symbol <- gsub("\"", "", listMat$HGNC.gene.symbol)
idx <- match(row.names(metaMat), listMat$HGNC.gene.symbol)
for(i in 4:ncol(listMat)){
  thislist <- listMat[,i][idx]
  thislist[is.na(thislist)] <- 0
  metaMat <- cbind(metaMat, thislist)
  colnames(metaMat)[dim(metaMat)[2]] <- colnames(listMat)[i]
}

## Epilepsy Ruzzo
listMat <- read.table(
  "../source/Gene_Lists/High-Confidence_Epilepsy_Risk_Genes_Ruzzo_2018-05-11.txt", fill = TRUE, header = TRUE, sep = "\t"
)
listMat$Epilepsy_All <- 1
listMat$Epilepsy_High_Conf <- 0
listMat$Epilepsy_High_Conf[listMat$Classification == "High-confidence"] <- 1
listMat$Epilepsy_High_Conf_Early <- 0
listMat$Epilepsy_High_Conf_Early[
  listMat$Classification == "High-confidence" &
  listMat$Early.Infantile == "Yes"] <- 1
listMat$Epilepsy_High_Conf_Late <- 0
listMat$Epilepsy_High_Conf_Late[
  listMat$Classification == "High-confidence" &
  listMat$Early.Infantile == "No"] <- 1
idx <- match(row.names(metaMat), listMat$Gene)
for(i in 4:ncol(listMat)){
  thislist <- listMat[,i][idx]
  thislist[is.na(thislist)] <- 0
  metaMat <- cbind(metaMat, thislist)
  colnames(metaMat)[dim(metaMat)[2]] <- colnames(listMat)[i]
}

## de novo ID: NEJM + Lancet
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_deLigt_NEJM.txt", header=T, sep="\t", fill=T)
nejmmat = listMat[listMat$nature_of_mutation=="D",1]
listMat <- read.table("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_Rauch_Lancet.txt", header=T, sep="\t", fill=T)
lancetmat = listMat[listMat$Type %in% c("frameshift", "nonsense", "splice"),1]
thislist = unique(c(lancetmat, nejmmat))
listname <- "DeNovoLGDsInID"
metaMat <- cbind(metaMat,rep(NA,nrow(metaMat)))
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat), thislist)
metaMat[!is.na(matchlist), listname] <- 1
metaMat[is.na(matchlist), listname] <- 0

## Biological Processes Gene lists from Luis (subset to ID)
Add_Gene_List_To_Binary_Matrix <- function(
  gene_binary_M, gene_list, gene_list_name){
  binary_gene_list <- as.numeric(rownames(gene_binary_M) %in% gene_list)
  gene_binary_M <- cbind(gene_binary_M, binary_gene_list)
  colnames(gene_binary_M)[dim(gene_binary_M)[2]] <- gene_list_name
  return(gene_binary_M)
}
Add_Gene_Lists_To_Binary_Matrix <- function(gene_binary_M, gene_lists_LL){
  for(i in 1:length(gene_lists_LL)){
    gene_binary_M <- Add_Gene_List_To_Binary_Matrix(
      gene_binary_M = gene_binary_M
      , gene_list = gene_lists_LL[[i]]
      , gene_list_name = names(gene_lists_LL)[i]
    )
  }
  return(gene_binary_M)
}
ORAgenelist <- read.csv("../source/Gene_Lists/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/GO_ORAGeneLists.csv", fill=T, header=T)
cluster_genes_LL <- split(
  ORAgenelist$hgnc_symbol
  , f = ORAgenelist$Category
)
cluster_genes_LL <- cluster_genes_LL[
  names(cluster_genes_LL) %in% c("ID.All", "ID.only")]
metaMat <- Add_Gene_Lists_To_Binary_Matrix(
  metaMat, gene_lists_LL = cluster_genes_LL
)

############
#Allen Human-Macaque-Specific Gene Expression Trajectories
############

inputdir = "../source/Gene_Lists/geschwindlabshares/atacseq/Enrichment_Lists/Bakken-AllenMacaque-Nature/";
files = dir(inputdir)
category = c("Primate-specific","Conserved","Human-specific","Rhesus-specific", "Not Conserved");

for(i in 1:length(files)){
  for(j in 1:length(category)){
    thislist = read.xlsx(paste(inputdir,files[[i]],sep=""),1, startRow=2, header=TRUE);
    ind = which(thislist$Set==category[j]);
    thislist = thislist[ind,1];
    matchlist = match(rownames(metaMat), thislist);
    matchlist = ifelse(is.na(matchlist),0,1);
    metaMat = cbind(metaMat, matchlist);
    colnames(metaMat)[dim(metaMat)[2]] = category[j];
  }
}

################################################
## Regulatory region datasets
################################################

## Luis ATAC enhancers
# All enhancers
listMat <- read.csv("../source/Gene_Lists/TorreUbieta_2018_TS2.csv"
  , header = TRUE)
thisgene = unique(listMat$hgnc_symbol)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "ATAC_enhancers"
# Enhancers CP>GZ
listMat <- read.csv("../source/Gene_Lists/TorreUbieta_2018_TS2.csv"
  , header = TRUE)
listMat <- listMat[listMat$log2FoldChange < 0 & listMat$padj.diffacc < 0.05, ]
thisgene = unique(listMat$hgnc_symbol)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "ATAC_enhancers_CP"
# Enhancers GZ>CP
listMat <- read.csv("../source/Gene_Lists/TorreUbieta_2018_TS2.csv"
  , header = TRUE)
listMat <- listMat[listMat$log2FoldChange > 0 & listMat$padj.diffacc < 0.05, ]
thisgene = unique(listMat$hgnc_symbol)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "ATAC_enhancers_GZ"

## Human gained enhancers
# All enhancers
listMat <- read.csv("../source/Gene_Lists/HGEs_Genes.csv"
  , header = TRUE)
thisgene = unique(listMat$HGNC)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "Human_Gained_Enhancers"
# Enhancers CP>GZ
listMat <- read.csv("../source/Gene_Lists/HGEs_Genes.csv"
  , header = TRUE)
listMat <- listMat[listMat$GZ.or.CP == "CP>GZ", ]
thisgene = unique(listMat$HGNC)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "Human_Gained_Enhancers_CP"
# Enhancers GZ>CP
listMat <- read.csv("../source/Gene_Lists/HGEs_Genes.csv"
  , header = TRUE)
listMat <- listMat[listMat$GZ.or.CP == "GZ>CP", ]
thisgene = unique(listMat$HGNC)
matchlist = match(rownames(metaMat), thisgene)
matchlist = ifelse(is.na(matchlist),0,1)
metaMat = cbind(metaMat, matchlist)
colnames(metaMat)[ncol(metaMat)] <- "Human_Gained_Enhancers_GZ"

################################################
## Fetal drop-seq
################################################

## Drop-seq Seurat cluster enriched genes
Load_Dropseq <- function(inDropseq, alignmentName) {
  ds2DF <- read.table(inDropseq, header = TRUE)
  ds2L <- split(ds2DF$Gene, f = ds2DF$Cluster)
  for(i in 1:length(ds2L)){
    thislist = ds2L[[i]]
    gename = geneAnno1[geneAnno1$hgnc_symbol %in% thislist, 2]
    gename = gename[gename!=""]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = paste0(alignmentName, names(ds2L)[i])
  }
  return(metaMat)
}
metaMat <- Load_Dropseq(
   inDropseq = "../analysis/tables/Seurat_ClusterDE/DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/res054/Seurat_ClusterDE_ClusterX_Vs_All_Clusters.txt"
  , alignmentName = "DS2-11 ")

### 2) Run the overlaps using gene length as a covariate
## These enrichments will be performed in a logistic regression framework

intGenes <- intersect(names(geneLength),names(exomeLength))
cor(geneLength[match(intGenes,names(geneLength))],
    exomeLength[match(intGenes,names(exomeLength))],
    method="spearman") ## Only 0.51 correlation (lower with Pearson) between coding and total gene length - this is interesting and good to keep in mind

exomeLength <- exomeLength[match(names(geneLength),names(exomeLength))]
dnvcols <- unique(c(grep("LGD",colnames(metaMat)), grep("denovo", colnames(metaMat)), grep("dnv", colnames(metaMat)),
                    grep("LOF", colnames(metaMat)), grep("mis", colnames(metaMat)), grep("disrupting", colnames(metaMat)),
                    grep("sibling", colnames(metaMat)), grep("proband", colnames(metaMat))))
genecols <- unique(c(grep("_gene", colnames(metaMat)), grep("gene", colnames(metaMat)),
grep("MAGMA", colnames(metaMat)), grep("clst", colnames(metaMat)), grep("Hi-C", colnames(metaMat))))

metaMat[metaMat>0] <- 1
table(as.numeric(data.matrix(metaMat))) ## Confirm all values are 1s and 0s

enrichOR <- enrichP <- matrix(NA,nrow=ncol(metaMat),ncol=ncol(metaMat))
colnames(enrichOR) <- colnames(enrichP) <- rownames(enrichOR) <- rownames(enrichP) <- colnames(metaMat)

for (i in 2:ncol(metaMat)) { ## Query Lists; Loop through columns of the data
  for (j in 2:ncol(metaMat)) { ## Biology/Interpretation Lists
    if (!is.na(match(j,dnvcols))) { ## if using de novo gene sets, use exome length
      thiscovar <- exomeLength
      print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
      dat1 <- as.numeric(metaMat[,i])
      dat2 <- as.numeric(metaMat[,j])
      glm.out <- glm(dat1~dat2+thiscovar,family=binomial)
      keepdat <- !is.na(dat1)&!is.na(dat2)

    } else if (!is.na(match(j, genecols)) | !is.na(match(i, genecols)) & j<120) {
    thiscovar <- geneLength ## For other gene sets, use gene length
    print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
    dat1 <- as.numeric(metaMat[,i])
    dat2 <- as.numeric(metaMat[,j])
    glm.out <- glm(dat1~dat2+thiscovar,family=binomial)
    keepdat <- !is.na(dat1)&!is.na(dat2)

    } else{
      print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
      dat1 <- as.numeric(metaMat[,i])
      dat2 <- as.numeric(metaMat[,j])
      glm.out <- glm(dat1~dat2,family=binomial)
      keepdat <- !is.na(dat1)&!is.na(dat2)
    }
    ## Compute glm that asks how well dat2 predicts dat1 when controlling for thiscovar

    enrichP[i,j] <- summary(glm.out)$coefficients[2,4]
    enrichOR[i,j] <- summary(glm.out)$coefficients[2,1]
  }
}

## Save the output in case you don't want to run all of the above again
dir.create("../analysis/analyzed_data/MetaMat", recursive = TRUE)
save(enrichP, enrichOR, metaMat
  , file="../analysis/analyzed_data/MetaMat/MetaMat.Rdata"
)
# load("../analysis/analyzed_data/MetaMat/MetaMat.Rdata")

Pmat <- enrichP
diag(Pmat) <- 0
Pmat[upper.tri(Pmat)] <- t(Pmat)[upper.tri(Pmat)] ## Fill in the other half, keep a symmetric matrix

ORmat <- enrichOR
diag(ORmat) <- Inf
ORmat[upper.tri(Pmat)] <- t(ORmat)[upper.tri(ORmat)] ## Fill in the other half, keep a symmetric matrix

#moduleSizes <- apply(metaMat,2,sum,na.rm=TRUE)


### 3) Take a subset of the relationships for plotting

###################################### Subset Target Lists here
keepcols2 = c(25:41);        #Interpretation Lists
keepcols1 = c(18:18);  #Target Lists

subORmat <- ORmat[keepcols1,keepcols2]
subPmat <- Pmat[keepcols1,keepcols2]

## FDR correction
subPmat <- p.adjust(subPmat, method="BH")   ##Correcting for all columns and rows

dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
rmval <- subPmat > 0.05 ## Display only values with p > 0.05
dispMat[rmval] <- 0
rownames(dispMat) <- rownames(subORmat)

#Use this with single query list vs all bio lists

#load("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/Enrichment/AllEnrichments-ATAC-HiC-Evo.Rdata")
#setwd("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/Enrichment")

#pdf("AllEnrichments-ATAC-HiC-Evo.pdf",height=5,width=20)
#FDRmat = subPmat; FDRfilt <- FDRmat<0.05

#textMat <- paste(signif(dispMat,2),"\n(",signif(FDRmat,2),")",sep="")
#textMat[!FDRfilt|subORmat==Inf] <- ""

#dispMat = matrix(dispMat, ncol=length(dispMat), nrow=1) #This is the color gradient for OR
#textMat = matrix(textMat, ncol=length(textMat), nrow=1) #This is the actual OR and p-values
#colnames(dispMat)=colnames(textMat)=colnames(metaMat)[-dim(metaMat)[2]]
#rownames(dispMat)=rownames(textMat)="ASDWGS"

#labeledHeatmap(Matrix=dispMat,
 #              textMatrix=textMat,
  #             xLabels=colnames(dispMat),
   #            yLabels=rownames(dispMat),
    #           colors=blueWhiteRed(100),
     #          cex.lab.x=0.4,
      #         cex.lab.y=0.4,
       #        zlim=c(-2,2),
        #       cex.text=0.2,
         #      setStdMargins=FALSE)
#dev.off()

#Use this with multiple query lists vs all bio lists

dir.create("../analysis/graphs/MetaMat", recursive = TRUE)
# pdf("../analysis/graphs/MetaMat/MetaMat_TADA_CorrelationHeatmap.pdf", height = 8, width = 12)
# FDRmat = subPmat; FDRfilt <- FDRmat<0.05
#
# # textMat <- matrix(paste(signif(dispMat,2),"\n(",signif(FDRmat,2),")",sep=""),nrow=nrow(dispMat),ncol=ncol(dispMat))
# textMat <- matrix(signif(dispMat,2),nrow=nrow(dispMat),ncol=ncol(dispMat))
# textMat[!FDRfilt|subORmat==Inf] <- ""
#
# # Transform FDR pvalues to -log10 and save as matrix for plotting
# # (dispMat is for plotting log2 Odds Ratio)
# dispMat2 <- matrix(FDRmat, nrow=nrow(dispMat), ncol=ncol(dispMat))
# dispMat2 <- -1*log(dispMat2, 10)
# dispMat2[is.infinite(dispMat2)] <- NA
#
# # Plot
# labeledHeatmap(Matrix=dispMat2,
#                textMatrix=textMat,
#                xLabels=colnames(dispMat),
#                yLabels=rownames(dispMat),
#                colors=colorRampPalette(c("white", "blue"))(100),
#                cex.lab.x = 1,
#                cex.lab.y = 1,
#                zlim=c(0, 100),
#                cex.text = 1,
#                setStdMargins = TRUE)
# dev.off();


## ASD TADA

keepcols2 = c(25:41);        #Interpretation Lists
keepcols1 = c(18:18);  #Target Lists

subORmat <- ORmat[keepcols1,keepcols2]
subPmat <- Pmat[keepcols1,keepcols2]

## FDR correction
subPmat <- p.adjust(subPmat, method="BH")   ##Correcting for all columns and rows

dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
rmval <- subPmat > 0.05 ## Display only values with p > 0.05
dispMat[rmval] <- 0
rownames(dispMat) <- rownames(subORmat)

## Barplot of pvalues
new_cluster_order <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
df1 <- melt(dispMat)
df1$Cluster <- as.factor(c(0:16))
df1$Cluster2 <- factor(df1$Cluster, levels = new_cluster_order)
# df1$Cluster2 <- factor(df1$Cluster2, levels = new_cluster_order)
df1$Pvalue_FDR <- subPmat
df1$Pvalue <- c(Pmat[keepcols1,keepcols2,drop=FALSE][1,])
df1 <- df1[df1$Cluster != 16, ]
# Label for coloring bars
df1 <- df1[order(df1$Cluster2), ]
df1$Class = factor(c(
      rep("Radial glia", 2)
      , rep("Cycling progenitor", 2)
      , "Intermediate progenitor"
      , rep("Excitatory Neuron", 5)
      , rep("Interneuron", 2)
      , "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia"
    )
  , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
    , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia")
)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)

ggplot(df1, aes(x = Cluster, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90) +
  ylim(0,11) +
  geom_hline(yintercept = -log(0.05, 10), color = "red")
ggsave("../analysis/graphs/MetaMat/MetaMat_TADA_BarPlot.pdf"
  , width = 8, height = 5)

# Reorder clusters
ggplot(df1, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90) +
  ylim(0,11) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_TADA_BarPlot_ReorderClusters.pdf"
  , width = 8, height = 5)

## Human specific

keepcols2 = c(25:42);        #Interpretation Lists
keepcols1 = c(22:22);  #Target Lists

subORmat <- ORmat[keepcols1,keepcols2,drop=FALSE]
subPmat <- Pmat[keepcols1,keepcols2,drop=FALSE]

# FDR correction
subPmat <- p.adjust(subPmat, method="BH")   ##Correcting for all columns and rows

dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
rmval <- subPmat > 0.05 ## Display only values with p > 0.05
dispMat[rmval] <- 0
rownames(dispMat) <- rownames(subORmat)

# Barplot of pvalues
new_cluster_order <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
df1 <- melt(dispMat)
df1$Cluster <- as.factor(c(0:16))
df1$Cluster2 <- factor(df1$Cluster, levels = new_cluster_order)
# df1$Cluster2 <- factor(df1$Cluster2, levels = new_cluster_order)
df1$Pvalue_FDR <- subPmat
df1$Pvalue <- c(Pmat[keepcols1,keepcols2,drop=FALSE][1,])
df1 <- df1[df1$Cluster != 16, ]
# Label for coloring bars
df1 <- df1[order(df1$Cluster2), ]
df1$Class = factor(c(
      rep("Radial glia", 2)
      , rep("Cycling progenitor", 2)
      , "Intermediate progenitor"
      , rep("Excitatory Neuron", 5)
      , rep("Interneuron", 2)
      , "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia"
    )
  , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
    , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia")
)

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)

# Plot
ggplot(df1, aes(x = Cluster, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90) +
  ylim(0,7) +
  geom_hline(yintercept = -log(0.05, 10), color = "red")
ggsave("../analysis/graphs/MetaMat/MetaMat_AllenHumanSpecific_BarPlot.pdf"
  , width = 8, height = 5)

# Plot with reordered clusters
ggplot(df1, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90) +
  ylim(0,7) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ggplot_set_theme_publication
ggsave(
  "../analysis/graphs/MetaMat/MetaMat_AllenHumanSpecific_ReorderClusters_BarPlot.pdf"
  , width = 8, height = 5)


## ASD iHART

keepcols2 = c(25:41);        #Interpretation Lists
keepcols1 = c(18:24);  #Target Lists

subORmat <- ORmat[keepcols1,keepcols2]
subPmat <- Pmat[keepcols1,keepcols2]

## FDR correction
subPmat <- p.adjust(subPmat, method="BH")   ##Correcting for all columns and rows
subPmat

dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
# rmval <- subPmat > 0.05 ## Display only values with p > 0.05
# dispMat[rmval] <- 0
rownames(dispMat) <- rownames(subORmat)

## Barplot of pvalues
new_cluster_order <- c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
df1 <- melt(dispMat)
df1$Pvalue_FDR <- subPmat
df1$Cluster <- gsub("DS2-11 ", "", df1$Var2)
df1$Cluster <- factor(df1$Cluster, levels = c(0:15))
df1$Cluster2 <- factor(df1$Cluster, levels = new_cluster_order)
# df1$Cluster2 <- factor(df1$Cluster2, levels = new_cluster_order)

pvals_DF <- melt(Pmat[keepcols1,keepcols2,drop=FALSE])
pvals_DF$Index <- paste0(pvals_DF$Var1, "_", pvals_DF$Var2)
df1$Index <- paste0(df1$Var1, "_", df1$Var2)
idx <- match(df1$Index, pvals_DF$Index)
df1$Pvalue <- pvals_DF$value[idx]

df1 <- df1[df1$Cluster != 16, ]
# Label for coloring bars
df1 <- df1[order(df1$Cluster2), ]
class_cluster_idx <- c(
  "Radial glia" = 7
  , "Radial glia" = 9
  , "Cycling progenitor" = 8
  , "Cycling progenitor" = 10
  , "Intermediate progenitor" = 2
  , "Excitatory Neuron" = 0
  , "Excitatory Neuron" = 1
  , "Excitatory Neuron" = 4
  , "Excitatory Neuron" = 3
  , "Excitatory Neuron" = 13
  , "Interneuron" = 5
  , "Interneuron" = 6
  , "Oligodendrocyte precursor" = 11
  , "Endothelial" = 12
  , "Pericyte" = 14
  , "Microglia" = 15
)
df1$Class <- names(class_cluster_idx)[match(df1$Cluster, class_cluster_idx)]
df1$Class <- factor(df1$Class
  , levels = c("Radial glia", "Cycling progenitor", "Intermediate progenitor"
    , "Excitatory Neuron", "Interneuron", "Oligodendrocyte precursor"
    , "Endothelial", "Pericyte", "Microglia")
)
# Remove NAs
df1 <- df1[! is.na(df1$Cluster), ]

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)

# iHART + Sanders
ggplot(df1, aes(x = Cluster, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Var1) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red")
ggsave("../analysis/graphs/MetaMat/MetaMat_iHART_Sanders_BarPlot.pdf"
  , width = 12, height = 7)
# Reorder clusters
ggplot(df1, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Var1) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_iHART_Sanders_BarPlot_ReorderClusters.pdf"
  , width = 10, height = 6)
# iHART 69
ggDF <- df1[df1$Var1 == "iHART.69", ]
ggplot(ggDF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_iHART69_BarPlot_ReorderClusters.pdf"
  , width = 8, height = 5)
# iHART
ggDF <- df1[df1$Var1 %in% c(
  "Sanders.65", "iHART.17novel", "iHART.69", "Sanders.not.replicated.13"
  , "PreviouslyKnownASDgene.52"), ]
ggplot(ggDF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Var1) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_iHART_BarPlot_ReorderClusters.pdf"
  , width = 10, height = 6)

## ATAC enhancers

# Interpretation Lists
interpretation_list_cols = c(31:47)
# Target Lists
target_list_cols = c(25:27)

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.border = element_blank()
)

dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- Format_OR_Pval_Matrices_For_Barplot(
  dispMat, subPmat, cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)

ggplot(gg_DF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  facet_wrap(~Var1) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_ATACenhancers_BarPlot_ReorderClusters.png"
  , width = 10, height = 4)


## Human gained enhancers

# Interpretation Lists
interpretation_list_cols = c(31:47)
# Target Lists
target_list_cols = c(28:30)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- Format_OR_Pval_Matrices_For_Barplot(
  dispMat, subPmat, cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
ggplot(gg_DF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  facet_wrap(~Var1) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_HumanGainedEnhancers_BarPlot_ReorderClusters.png"
  , width = 10, height = 4)

# Interpretation Lists
interpretation_list_cols = c(18:39)
# Target Lists
target_list_cols = c(50:52)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- melt(dispMat)
gg_DF$Cluster <- gg_DF$Var2
# Add p-values
gg_DF$Pvalue_FDR <- subPmat
pvals_DF <- melt(Pmat[target_list_cols, interpretation_list_cols, drop=FALSE])
pvals_DF$Index <- paste0(pvals_DF$Var1, "_", pvals_DF$Var2)
gg_DF$Index <- paste0(gg_DF$Var1, "_", gg_DF$Var2)
idx <- match(gg_DF$Index, pvals_DF$Index)
gg_DF$Pvalue <- pvals_DF$value[idx]
ggplot(gg_DF, aes(x = Cluster, y = -log(Pvalue, 10), fill = Cluster)) +
  facet_wrap(~Var1) +
  geom_bar(stat = "identity", position = position_dodge()) +
  # scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("../analysis/graphs/MetaMat/MetaMat_HumanGainedEnhancers_Nowakowski_BarPlot_ReorderClusters.png"
  , width = 13, height = 6)


## Nowakowski

# Interpretation Lists
interpretation_list_cols = c(18:39)
# Target Lists
target_list_cols = c(53:69)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- melt(dispMat)
gg_DF$Cluster <- gg_DF$Var2
# Add p-values
gg_DF$Pvalue_FDR <- subPmat
pvals_DF <- melt(Pmat[target_list_cols, interpretation_list_cols, drop=FALSE])
pvals_DF$Index <- paste0(pvals_DF$Var1, "_", pvals_DF$Var2)
gg_DF$Index <- paste0(gg_DF$Var1, "_", gg_DF$Var2)
idx <- match(gg_DF$Index, pvals_DF$Index)
gg_DF$Pvalue <- pvals_DF$value[idx]
# Make 0 p-values smallest possible value for log transform
gg_DF$Pvalue[gg_DF$Pvalue == 0] <- 1.333258e-300
gg_DF$Var1 <- gsub("DS2-11 ", "", gg_DF$Var1)
gg_DF$Var1 <- factor(gg_DF$Var1
  , levels = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15))
ggplot(gg_DF, aes(x = Var1, y = Var2, fill = -log(Pvalue, 10))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", space = "Lab", name = "-log10 p-value") +
  geom_text(aes(x = Var1, y = Var2, label = signif(value, 2))
    , color = "black", size = 3) +
  ylab("Nowakowski et al. 2017") +
  xlab("Drop-seq") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_Nowakowski_Heatmap.png"
  , width = 9, height = 7)
################################################################################

### Epilepsy

# Interpretation Lists
interpretation_list_cols = c(60:76)
# Target Lists
target_list_cols = c(47:50)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- Format_OR_Pval_Matrices_For_Barplot(
  dispMat, subPmat, cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
ggplot(gg_DF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  facet_wrap(~Var1, scales = "free_x") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_Epilepsy_BarPlot_ReorderClusters.png"
  , width = 10, height = 8)
ggsave("../analysis/graphs/MetaMat/MetaMat_Epilepsy_BarPlot_ReorderClusters.pdf"
  , width = 10, height = 8)

# For paper
# Interpretation Lists
interpretation_list_cols = c(60:76)
# Target Lists
target_list_cols = c(48:50)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- Format_OR_Pval_Matrices_For_Barplot(
  dispMat, subPmat, cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
ggplot(gg_DF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  facet_wrap(~Var1, scales = "free_x") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave(
  "../analysis/graphs/MetaMat/MetaMat_Epilepsy_BarPlot_ReorderClusters_paper.pdf"
  , width = 12, height = 4)
################################################################################

### ID

# Interpretation Lists
interpretation_list_cols = c(60:76)
# Target Lists
target_list_cols = c(51:53)
dispMat <- Subset_OR_And_Log2_Odds_Ratio_Matrix(
  target_list_cols, interpretation_list_cols, ORmat
)
subPmat <- Subset_OR_Pval_And_FDR_Correct_Matrix(
  target_list_cols, interpretation_list_cols, Pmat
)
gg_DF <- Format_OR_Pval_Matrices_For_Barplot(
  dispMat, subPmat, cluster_order = c(9,7,8,10,2,0,1,4,3,13,5,6,11,12,14,15)
)
ggplot(gg_DF, aes(x = Cluster2, y = -log(Pvalue, 10), fill = Class)) +
  facet_wrap(~Var1, scales = "free_x") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(type = "qual", palette = "Set3", direction = 1) +
  geom_text(aes(x = Cluster2, y = -log(Pvalue, 10), label = round(value, 1))
    , position = position_dodge(width = 1), hjust = -0.3, color = "black"
    , angle = 90, size = 3) +
  # ylim(0,15) +
  geom_hline(yintercept = -log(0.05, 10), color = "red") +
  ylab("-Log10 p-value") +
  xlab("Reordered clusters") +
  ggplot_set_theme_publication
ggsave("../analysis/graphs/MetaMat/MetaMat_ID_BarPlot_ReorderClusters.png"
  , width = 10, height = 5)
ggsave("../analysis/graphs/MetaMat/MetaMat_ID_BarPlot_ReorderClusters.pdf"
  , width = 10, height = 5)
################################################################################
