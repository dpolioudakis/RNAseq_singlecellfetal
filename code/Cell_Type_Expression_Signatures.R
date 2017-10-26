# Damon Polioudakis
# 2017-05-09
# 
################################################################################

rm(list = ls())

require(Seurat)
require(WGCNA)
require(GGally)
require(ggplot2)

# load("../analysis/Cell_Type_Expression_Signatures_Workspace.RData")

# Drop-seq
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
# Cluster markers
mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , header = TRUE)

# Pollen
pnExDF <- read.csv("../pollen_2015/data/Pollen_2015_ST2.csv", header = TRUE)
pnMtDF <- read.csv("../pollen_2015/metadata/Cell paper - updated attributes.csv"
  , header = TRUE)

# Kang
load("../neurogenesis/analysis/Kang_Expression_RgCv.RData")
# nmExDF
# sampleKey is ProcessedKangMetaData.csv
KangAnnotRAW <- read.csv("../neurogenesis/orig.data/InVivoData/annot.csv"
  , row.names = 1)

# Miller
load("../neurogenesis/orig.data/LCMDE/AllenLCM.Rdata")
mlExDF = AllenLCM$datExpr
mlMetDF = AllenLCM$datTraits
Zones = read.csv("../neurogenesis/orig.data/LCMDE/LCM_Zones_CPio.csv")
MillerAnnotRAW = read.csv("../neurogenesis/orig.data/LCMDE/annot.csv", row.names = 1)

# Gene info from biomaRt
geneinfoDF <- read.csv("../source/BiomaRt_Compile_GeneInfo_GRCh37_Ensembl75.csv"
  , header = TRUE)

## Variables
outGraph <- "../analysis/graphs/Cell_Type_Expression_Signatures_"
outAnalysis <- "../analysis/Cell_Type_Expression_Signatures_"
################################################################################

### Format

## Pollen

row.names(pnExDF) <- pnExDF$X
pnExDF <- pnExDF[-1, ]

## Kang

# Kang probe IDs? to entrez
row.names(nmExDF) <- KangAnnotRAW$ENTREZ_ID[
  match(row.names(nmExDF), row.names(KangAnnotRAW))]
# Entrez to hgnc_symbol
row.names(nmExDF) <- geneinfoDF$hgnc_symbol[
  match(row.names(nmExDF), geneinfoDF$entrezgene)]

# Subset to Stage 1-6
nmExDF <- nmExDF[ ,colnames(nmExDF) %in% sampleKey$Sample.ID[sampleKey$Stage %in% c(1:6)]]

## Miller

# Select cortical layers only

ZonesMiller = mlMetDF$structure_acronym

mlMetDF$Labels = rep(NA, nrow(mlMetDF))
for (Label in Zones$Label)
  mlMetDF[grep(Label, ZonesMiller), "Labels"] = Label

Zoneindex = which(! is.na(mlMetDF$Labels))

mlExDF = mlExDF[, Zoneindex]
mlMetDF = mlMetDF[Zoneindex, ]

# Metadata character class to factor
mlMetDF[sapply(mlMetDF, is.character)] <- lapply(mlMetDF[sapply(mlMetDF, is.character)]
  , as.factor)

# Format zone codes
millerZones <- mlMetDF[match(mlMetDF$well_id, colnames(mlExDF)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)

mlMetDF$Labels <- millerZones

# Only ENTREZ annotated probes
mlExDF = mlExDF[! is.na(MillerAnnotRAW$ENTREZ_ID), ]
MillerAnnot = MillerAnnotRAW[! is.na(MillerAnnotRAW$ENTREZ_ID), ]

# Get maximum expression probe
Genes = unique(MillerAnnot$ENTREZ_ID)
keepind = matrix(nrow = 0, ncol = 0);
for (ii in 1:length(Genes)) {
  genematchind = which(MillerAnnot$ENTREZ_ID == Genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(mlExDF[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}
mlExDF = mlExDF[keepind, ];
rownames(mlExDF) = Genes;
MillerAnnot = MillerAnnot[keepind, ];
################################################################################

# Miller
genes <- mkDF$gene[mkDF$cluster == 7][1:20]
genes <- geneinfoDF$entrez[match(genes, geneinfoDF$hgnc_symbol)]
colors <- rep(NA, nrow(mlExDF))
colors[row.names(mlExDF) %in% genes] <- 7
mlME = moduleEigengenes(t(mlExDF), colors = colors)
cmDF <- data.frame(cor(t(mlExDF), mlME$eigengenes[1]))
cmDF$GENE <- row.names(cmDF)
cmDF <- cmDF[order(-cmDF$ME7), ]
# Entrez to hgnc_symbol
cmDF$GENE <- geneinfoDF$hgnc_symbol[match(cmDF$GENE, geneinfoDF$entrezgene)]

# Kang
genes <- mkDF$gene[mkDF$cluster == 7][1:20]
colors <- rep(NA, nrow(nmExDF))
colors[row.names(nmExDF) %in% genes] <- 7
kaME = moduleEigengenes(t(nmExDF), colors = colors)
ckDF <- cor(t(nmExDF), kaME$eigengenes[1])
ckDF <- data.frame(GENE = row.names(ckDF), ME7 = ckDF)
ckDF <- ckDF[order(-ckDF$ME7), ]

# Drop-seq
genes <- mkDF$gene[mkDF$cluster == 7][1:20]
colors <- rep(NA, nrow(fetb@scale.data))
colors[row.names(fetb@scale.data) %in% genes] <- 7
# exDF <- fetb@scale.data[ ,fetb@ident == 7]
exDF <- fetb@scale.data
dsME = moduleEigengenes(t(exDF), colors = colors)
cdDF <- cor(t(exDF), dsME$eigengenes[1])
cdDF <- data.frame(GENE = row.names(cdDF), ME7 = cdDF)
cdDF <- cdDF[order(-cdDF$ME7), ]

# Pollen
genes <- mkDF$gene[mkDF$cluster == 7][1:20]
colors <- rep(NA, nrow(pnExDF))
colors[row.names(pnExDF) %in% genes] <- 7
# exDF <- pnExDF[ ,colnames(pnExDF) %in% gsub("-", ".", pnMtDF$Cell[pnMtDF$Inferred.Cell.Type == "RG"])]
exDF <- pnExDF
pnME = moduleEigengenes(t(exDF), colors = colors)
cpDF <- cor(t(exDF), pnME$eigengenes[1])
cpDF <- data.frame(GENE = row.names(cpDF), ME7 = cpDF)
cpDF <- cpDF[order(-cpDF$ME7), ]

save(list = ls(), file = paste0(outAnalysis, "Workspace.RData"))

length(intersect(as.character(cmDF$GENE[1:1000]), as.character(ckDF$GENE)[1:1000]))
length(intersect(as.character(cpDF$GENE[1:1000]), as.character(ckDF$GENE)[1:1000]))
ckDF$GENE <- as.character(ckDF$GENE)
row.names(ckDF) <- 1:nrow(ckDF)
cmDF$GENE <- as.character(cmDF$GENE)
row.names(cmDF) <- 1:nrow(cmDF)
df <- merge(ckDF, cmDF, by.x = "GENE1", by.y = "GENE", all = FALSE)
ggplot(df, aes(x = ME7.x, y = ME7.y)) +
  geom_point(alpha = 0.25, size = 0.25)
ggsave(paste0(outGraph, "test.png"), dpi = 100)
cor(df[ ,2:3], method = "spearman")
cor(df[ ,2:3], method = "pearson")

df <- merge(cdDF, cpDF, by = "row.names")
df <- merge(df, cmDF, by.x = "Row.names", by.y = "GENE")
df <- merge(df, ckDF, by.x = "Row.names", by.y = "GENE")
df <- df[ ,-c(2,4)]
colnames(df) <- c("GENE", "DROPSEQ", "POLLEN", "MILLER", "KANG")
png(paste0(outGraph, "Correlation_Matrix.png"))
ggpairs(df, columns = 2:5)
dev.off()

df <- merge(cdDF, ckDF, by.x = "row.names", by.y = "GENE")
cor(df[ ,3:4], use = 'pairwise.complete.obs', method = "spearman")

