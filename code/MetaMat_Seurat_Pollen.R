
# Damon Polioudakis
# 2016-10-25
# Derived from Luis' script:
#   /geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/Enrichment/Metamat-CorticalLayers.R
################################################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

library(biomaRt) ## First, construct a matrix containing all possible protein coding genes. We will use human gene symbols. Although other identifiers, such as ENSG IDs are preferred, many studies only report human gene symbols, so we will use this as the common identifier.
library(WGCNA) ## For plotting heatmaps

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
## b) Single-cell RNA-seq Datasets
################################################
## Cell types in fetal brain (Pollen et al S3)
load(file="../source/metaMat/PollenS3GenesENSG.Rda");
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
load(file="../source/metaMat/PollenCellTypesGenes.R");
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
load(file="../source/metaMat/PollenCellTypesGenes-Unique.R")
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

## Iossifov ASD - load for exomeLength
listMat <- read.csv("../source/metaMat/AllDeNovoByGene.csv") ## From the supplmental tables of Iossifov et al., 2014 - note that some gene symbols have the excel conversion error in this data, but it will get removed when we overlap things
listMat <- listMat[!duplicated(listMat[,1]),]
rownames(listMat) <- listMat[,"gene"]
exomeLength <- listMat[,"codingLenInTarget"] ## For later use in making the covariate matrix
names(exomeLength) <- rownames(listMat)

## Drop-seq Seurat cluster marker genes
Load_Dropseq <- function (inDropseq, alignmentName) {
  ds2DF <- read.table(inDropseq, header = TRUE)
  ds2L <- split(ds2DF$gene, f = ds2DF$cluster)
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
metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , "Exon ")
metaMat <- Load_Dropseq("../analysis/tables/Seurat_Markers_DS002_003_exon_FtMm250_200gd_2500dg_Mt5_Pc110_Res3_Marker_Genes_Clusters_Vs_All.txt"
  , "PC1-10 Res 3.0 ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample50p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon ds50p ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample25p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon ds25p ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_downSample10p_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon ds10p ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_ExonIntron_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon + Intron ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_ExonIntron_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon + Intron ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_DsCl50p_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon DsCell50p ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_DsCl25p_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon DsCell25p ")
# metaMat <- Load_Dropseq("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_DsCl10p_Marker_Genes_Clusters_Vs_All.txt"
#   , "Exon DsCell10p ")

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
colnames(enrichOR) <- rownames(enrichP) <- rownames(enrichOR) <- rownames(enrichP) <- colnames(metaMat)


for (i in 18:44) { ## Query Lists; Loop through columns of the data
  for (j in 2:17) { ## Biology/Interpretation Lists
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
save(enrichP, enrichOR, metaMat, file="../analysis/MetaMat_Seurat_Pollen.Rdata")

Pmat <- enrichP
diag(Pmat) <- 0
Pmat[upper.tri(Pmat)] <- t(Pmat)[upper.tri(Pmat)] ## Fill in the other half, keep a symmetric matrix

ORmat <- enrichOR
diag(ORmat) <- Inf
ORmat[upper.tri(Pmat)] <- t(ORmat)[upper.tri(ORmat)] ## Fill in the other half, keep a symmetric matrix

#moduleSizes <- apply(metaMat,2,sum,na.rm=TRUE)


### 3) Take a subset of the relationships for plotting

###################################### Subset Target Lists here
keepcols2 = c(2:17);        #Interpretation Lists
keepcols1 = c(18:44);  #Target Lists

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

pdf("../analysis/graphs/MetaMat_Seurat_Pollen.pdf", height = 8, width = 12)
FDRmat = subPmat; FDRfilt <- FDRmat<0.05

# textMat <- matrix(paste(signif(dispMat,2),"\n(",signif(FDRmat,2),")",sep=""),nrow=nrow(dispMat),ncol=ncol(dispMat))
textMat <- matrix(signif(dispMat,2),nrow=nrow(dispMat),ncol=ncol(dispMat))
textMat[!FDRfilt|subORmat==Inf] <- ""

# Transform FDR pvalues to -log10 and save as matrix for plotting
# (dispMat is for plotting log2 Odds Ratio)
dispMat2 <- matrix(FDRmat, nrow=nrow(dispMat), ncol=ncol(dispMat))
dispMat2 <- -1*log(dispMat2, 10)
dispMat2[is.infinite(dispMat2)] <- NA

# Plot
labeledHeatmap(Matrix=dispMat2,
               textMatrix=textMat,
               xLabels=colnames(dispMat),
               yLabels=rownames(dispMat),
               colors=colorRampPalette(c("white", "blue"))(100),
               cex.lab.x = 1,
               cex.lab.y = 1,
               zlim=c(0, 100),
               cex.text = 1,
               setStdMargins = TRUE)
dev.off();
