# Damon Polioudakis
# 2016-12-08
# Plot modules defined in Kang dataset across Kang, Miller, Luis VZ/CP, hNPCs,
# and Primary cultures
################################################################################

rm(list = ls())

require(WGCNA)
require(ggplot2)
require(reshape2)
require(biomaRt)
require(gridExtra)

options(stringsAsFactors = F)

## Load data and assign variables

# Drop-seq cluster marker gene lists
mkDF <- read.table("../analysis/tables/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_Marker_Genes_Clusters_Vs_All.txt"
  , header = TRUE)

# Kang
# Non regressed data
load("../neurogenesis/orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
KangExpr = t(datExpr)
rm(datExpr)
# Regressed data
load("../neurogenesis/analysis/Kang_Expression_RgCv.RData")
KangExpr <- nmExDF
# sampleKey is ProcessedKangMetaData.csv
metKangDF <- sampleKey
# metKangDF <- read.csv("../neurogenesis/orig.data/InVivoData/ProcessedKangMetaData.csv")
KangAnnotRAW = read.csv("../neurogenesis/orig.data/InVivoData/annot.csv", row.names = 1)
KangAnnotnet2 = KangAnnotRAW[which(rownames(KangAnnotRAW) %in% rownames(KangExpr)), ] 
load("../neurogenesis/orig.data/InVivoData/net2_power16_cutHeigt0.15.RData")

# Miller
load("../neurogenesis/orig.data/LCMDE/AllenLCM.Rdata")
MillerExprRAW = AllenLCM$datExpr
MillerMetaRAW = AllenLCM$datTraits
Zones = read.csv("../neurogenesis/orig.data/LCMDE/LCM_Zones_CPio.csv")
MillerAnnotRAW = read.csv("../neurogenesis/orig.data/LCMDE/annot.csv", row.names = 1)

## Variables
graphCodeTitle <- "Kang_ModuleEG_Expression.R"
outGraphPfx <- "../analysis/graphs/Kang_ModuleEG_Expression_"

## Output Directory
dir.create("../analysis/graphs", recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

## Function: Convert HGNC symbol to Entrez and Ensembl IDs
# Requires biomaRt package: require(biomaRt)
ConvertGeneSymbolToEntrez <- function (genesList) {
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
    attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene")
    , filters = "hgnc_symbol"
    , values = genesDF
    , mart = ensembl
  )
  bmDF
}
################################################################################

### Format and subset of data

## Cluster marker genes gene symbols to Entrez

enzDF <- ConvertGeneSymbolToEntrez(mkDF$gene)
mkDF <- merge(mkDF, enzDF, by.x = "gene", by.y = "hgnc_symbol")

## Format and subset Miller

ZonesMiller = MillerMetaRAW$structure_acronym

# Change name of Kang metadata sample column from X to SampleID
names(metKangDF)[1] = "SampleID"

## Select cortical layers only

MillerMetaRAW$Labels = rep(NA, nrow(MillerMetaRAW))
for (Label in Zones$Label)
  MillerMetaRAW[grep(Label, ZonesMiller), "Labels"] = Label

Zoneindex = which(! is.na(MillerMetaRAW$Labels))

MillerExpr = MillerExprRAW[, Zoneindex]
MillerMeta = MillerMetaRAW[Zoneindex, ]


## Only ENTREZ annotated probes

MillerExpr = MillerExpr[! is.na(MillerAnnotRAW$ENTREZ_ID), ]
MillerAnnot = MillerAnnotRAW[! is.na(MillerAnnotRAW$ENTREZ_ID), ]


## Get maximum expression probe

# Miller Dataset

Genes = unique(MillerAnnot$ENTREZ_ID)

keepind = matrix(nrow = 0, ncol = 0);

for (ii in 1:length(Genes)) {
  genematchind = which(MillerAnnot$ENTREZ_ID == Genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(MillerExpr[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}

MillerExpr = MillerExpr[keepind, ];
rownames(MillerExpr) = Genes;
MillerAnnot = MillerAnnot[keepind, ];


## Select only ENTREZ genes evaluated in the Kang datasets

KangAnnot = KangAnnotnet2

Index = MillerAnnot$ENTREZ_ID %in% KangAnnot$ENTREZ_ID

MillerAnnot = MillerAnnot[Index, ]
MillerExpr = MillerExpr[Index, ]

nrow(MillerExpr)
# 16534

## Limit to ENTREZ IDs present in all datasets
## Arbitrary for Module preservation analysis, but mandatory for follow up
## procedures

if (TRUE) {
  ENTREZset = intersect(
    intersect(KangAnnotnet2$ENTREZ_ID, MillerAnnot$ENTREZ_ID)
        , mkDF$entrezgene)
  
  Index = match(ENTREZset, KangAnnotnet2$ENTREZ_ID)
  KangExpr = KangExpr[Index, ]
  KangAnnot = KangAnnotnet2[Index, ]
  
  Index = match(ENTREZset, MillerAnnot$ENTREZ_ID)
  MillerExpr = MillerExpr[Index, ]
  MillerAnnot = MillerAnnot[Index, ]
  
  # Use %in% because gene can be marker for more than 1 cluster
  # (match only returns position of 1 match)
  Index = mkDF$entrezgene %in% ENTREZset
  mkDF = mkDF[Index, ]
}
################################################################################

### Plot mean expression of cluster markers across Kang and Miller

# Remove genes in multiple cluster marker lists (should probably get rid of this step)
mkDF <- mkDF[! duplicated(mkDF$entrezgene), ]

## Calculate MEs
# Kang
clusters <- mkDF$cluster[match(KangAnnot$ENTREZ_ID, mkDF$entrezgene)]
df <- data.frame(KangExpr)
df$CLUSTERS <- clusters
kangMns <- data.frame(lapply(split(df, df$CLUSTERS), colMeans))
kangMns <- kangMns[-nrow(kangMns), ]
# Miller
clusters <- mkDF$cluster[match(row.names(MillerExpr), mkDF$entrezgene)]
df <- data.frame(MillerExpr)
df$CLUSTERS <- clusters
millerMns <- data.frame(lapply(split(df, df$CLUSTERS), colMeans))
millerMns <- millerMns[-nrow(millerMns), ]

# Subset Kang metadata to Kang samples present in expression DF
sbMetKangDF <- metKangDF[match(colnames(KangExpr), metKangDF$SampleID), ]
sbMetKangDF <- sbMetKangDF[! is.na(sbMetKangDF[1]), ]

## Extract stage for each sample, ordered by expression DF - same order as
## module eigengene list
# Kang
kangStage <- sbMetKangDF[match(sbMetKangDF$SampleID, colnames(KangExpr)), ]$Stage
# Miller
millerZones <- MillerMeta[match(MillerMeta$well_id, colnames(MillerExpr)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)

# Kang
df <- data.frame(Stage = kangStage, kangMns)
kangGgDF <- melt(df)
kangGgDF <- melt(df, id.vars = "Stage"
  , variable.name = "Mn", value.name = "Expression")
kangGgLDF <- split(kangGgDF, kangGgDF$Mn)
kangGGL <- lapply(names(kangGgLDF), function(name) {
  mnGgDF <- aggregate(kangGgLDF[[name]], list(kangGgLDF[[name]]$Stage), mean)
  ggplot(kangGgLDF[[name]], aes(x = Stage, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") +
    geom_line(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") + 
    xlab("Stage\n") +
    ylab("Mean expression") +
    ggtitle(paste0(name, "\n"))
})

# Miller
df <- data.frame(Zone = millerZones, millerMns)
df$Zone <- factor(df$Zone, levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG"))
millerGgDF <- melt(df)
millerGgDF <- melt(df, id.vars = "Zone"
  , variable.name = "Mn", value.name = "Expression")
millerGgLDF <- split(millerGgDF, millerGgDF$Mn)
millerGGL <- lapply(names(millerGgLDF), function(name) {
  mnGgDF <- aggregate(millerGgLDF[[name]], list(millerGgLDF[[name]]$Zone), mean)
  mnGgDF$Zone <- mnGgDF$Group.1
  ggplot(millerGgLDF[[name]], aes(x = Zone, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") +
    geom_line(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") + 
    xlab("Zone\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Combine to set layout of grid.arrange
ggLL <- mapply(list, kangGGL, millerGGL)

pdf(paste0(outGraphPfx, "Split_Line_Graphs_KangRegCv_Means.pdf"), height = 24, width = 6)
do.call("grid.arrange", c(ggLL, ncol = 2))
dev.off()
################################################################################
### Plot ME expression from Kang modules in datasets

mkDF <- mkDF[! duplicated(mkDF$entrezgene), ]

## Calculate MEs
# Kang
clusters <- mkDF$cluster[match(KangAnnot$ENTREZ_ID, mkDF$entrezgene)]
kangMEs <- moduleEigengenes(t(KangExpr), clusters)$eigengenes
# Miller
clusters <- mkDF$cluster[match(row.names(MillerExpr), mkDF$entrezgene)]
millerMEs <- moduleEigengenes(t(MillerExpr), clusters)$eigengenes

# Subset Kang metadata to Kang samples present in expression DF
sbMetKangDF <- metKangDF[match(colnames(KangExpr), metKangDF$SampleID), ]
sbMetKangDF <- sbMetKangDF[! is.na(sbMetKangDF[1]), ]

## Extract stage for each sample, ordered by expression DF - same order as
## module eigengene list
# Kang
kangStage <- sbMetKangDF[match(sbMetKangDF$SampleID, colnames(KangExpr)), ]$Stage
# Miller
millerZones <- MillerMeta[match(MillerMeta$well_id, colnames(MillerExpr)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)

# Kang
df <- data.frame(Stage = kangStage, kangMEs)
kangGgDF <- melt(df)
kangGgDF <- melt(df, id.vars = "Stage"
  , variable.name = "ME", value.name = "Expression")
kangGgLDF <- split(kangGgDF, kangGgDF$ME)
kangGGL <- lapply(names(kangGgLDF), function(name) {
  mnGgDF <- aggregate(kangGgLDF[[name]], list(kangGgLDF[[name]]$Stage), mean)
  ggplot(kangGgLDF[[name]], aes(x = Stage, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") +
    geom_line(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") + 
    xlab("Stage\n") +
    ylab("ME expression") +
    ggtitle(paste0(name, "\n"))
})

# Miller
df <- data.frame(Zone = millerZones, millerMEs)
df$Zone <- factor(df$Zone, levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG"))
millerGgDF <- melt(df)
millerGgDF <- melt(df, id.vars = "Zone"
  , variable.name = "ME", value.name = "Expression")
millerGgLDF <- split(millerGgDF, millerGgDF$ME)
millerGGL <- lapply(names(millerGgLDF), function(name) {
  mnGgDF <- aggregate(millerGgLDF[[name]], list(millerGgLDF[[name]]$Zone), mean)
  mnGgDF$Zone <- mnGgDF$Group.1
  ggplot(millerGgLDF[[name]], aes(x = Zone, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") +
    geom_line(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") + 
    xlab("Zone\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Combine to set layout of grid.arrange
ggLL <- mapply(list, kangGGL, millerGGL)

pdf(paste0(outGraphPfx, "Split_Line_Graphs_KangRegCv_MEs.pdf"), height = 24, width = 6)
do.call("grid.arrange", c(ggLL, ncol = 2))
dev.off()
################################################################################


# # Calculate mean ME expression at each time point or brain region
# Dev_Time_Point_Mean_MEexpr <- function (timePoint, MEs) {
#   df <- cbind(data.frame(Time_Point = timePoint), MEs)
#   tmpL <- split(df, df$Time_Point)
#   # Remove Time_Point column
#   tmpL <- lapply(tmpL, function(df) df[ ,-1])
#   outM <- sapply(tmpL, function(timePoint) {
#     apply(timePoint, 2, mean)
#   })
#   outM
# }
# 
# # Kang
# ggM <- Dev_Time_Point_Mean_MEexpr(kangStage, kangMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# kangGgDF <- melt(ggDF, id.vars = "Time_Point"
#                  , variable.name = "ME", value.name = "Expression")
# kangGgDF$Data_Set <- "Kang"
# 
# # hNPs
# ggM <- Dev_Time_Point_Mean_MEexpr(hNPwK, hNPMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# gghNPdF <- melt(ggDF, id.vars = "Time_Point"
#                 , variable.name = "ME", value.name = "Expression")
# gghNPdF$Data_Set <- "hNP"
# 
# # Miller
# ggM <- Dev_Time_Point_Mean_MEexpr(millerZones, millerMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# # Convert zones to numeric for plotting on Kang time course
# ggDF$Time_Point <- as.numeric(factor(ggDF$Time_Point
#                                      , levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG")))
# ggMillerdF <- melt(ggDF, id.vars = "Time_Point"
#                    , variable.name = "ME", value.name = "Expression")
# ggMillerdF$Data_Set <- "Miller"
# 
# # Combine for ggplot2 plotting
# ggDF <- rbind(kangGgDF, gghNPdF, ggMillerdF)
# ggDF$Time_Point <- as.numeric(ggDF$Time_Point)
# 
# ggplot(ggDF, aes(x = Time_Point, y = Expression, color = Data_Set)) +
#   facet_wrap(~ME) +
#   geom_line() +
#   xlab("Stage / Week / Zone") +
#   ylab("ME Expression") +
#   theme(axis.text.x = element_blank()) +
#   ggtitle(paste0(graphCodeTitle
#                  ,"\nKang Module Eigengene Expression"
#                  ,"\nX-axis:"
#                  ,"\nKang: Stage 1, 2, 3, 4, 5, 6, 7, 8"
#                  ,"\nhNPs: 1, 4, 8, 12 weeks"
#                  ,"\nMiller: VZ, SZ, IZ, SP, CP, MZ, SG"
#                  ,"\n"))
# ggsave(paste0(outGraphPfx, "Line_Graphs.pdf"), width = 14, height = 12)
################################################################################