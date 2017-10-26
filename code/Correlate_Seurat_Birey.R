

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(dplyr)
require(Matrix)


load("../analysis/DS002003_exon_FtMm250_Seurat_NoScale.Robj")

biExDF <- read.csv("../source/Birey_2017_DataTable.csv")
biNmExDF <- read.csv("../source/Birey_2017_NormalizedDataTable.csv")
biMtDF <- read.csv("../source/Birey_2017_AnnotationTable.csv")

### Number of each cell class

df <- data.frame(table(biMtDF$X9.final.clusters..one.empty.))
df <- df[-9, ]
df$Percent <- round(df$Freq/sum(df$Freq) * 100, 1)
df$Var1 <- c("1 Glutamatergic Neurons", "2 Intermediate Progenitors"
  , "3 Radial Glia", "4 Astroglia", "5 Ventral Progenitors"
  ,  "6 GABAergic Neurons", "7 OPC", "8 Choroid Plexus")

df2 <- data.frame(table(seuratO@ident))
df2$Percent <- round(df2$Freq/sum(df2$Freq) * 100, 1)

df <- rbind(df, df2)

## Subset Pasca to pallium
biExDF <- biExDF[biMtDF$sample == "P", ]
biNmExDF <- biNmExDF[biMtDF$sample == "P", ]
biMtDF <- biMtDF[biMtDF$sample == "P", ]

df3 <- data.frame(table(biMtDF$X9.final.clusters..one.empty.))
df3 <- df3[-9, ]
df3$Percent <- round(df3$Freq/sum(df3$Freq) * 100, 1)
df3$Var1 <- c("1 Glutamatergic Neurons", "2 Intermediate Progenitors"
  , "3 Radial Glia", "4 Astroglia", "5 Ventral Progenitors"
  ,  "6 GABAergic Neurons", "7 OPC", "8 Choroid Plexus")

df <- rbind(df, df3)

write.csv(df, "../analysis/tables/Correlate_Seurat_Birey_CellsPerCluster.csv"
  , quote = FALSE)

### Correlation

## Raw counts

# Geschwind
exDF <- as.matrix(seuratO@raw.data)
exDF <- exDF[ ,colnames(exDF) %in% colnames(as.matrix(seuratO@scale.data))]
df <- data.frame(seuratO@ident, t(exDF))
mnDF <- sapply(split(df[,-1], df$seuratO.ident), function(df) {colMeans(df)})

# Pasca
df <- data.frame(biMtDF["X9.final.clusters..one.empty."], biExDF)
mnBiDF <- sapply(split(df[,-1], df$X9.final.clusters..one.empty.), function(df) {colMeans(df)})
mnBiDF <- mnBiDF[ ,-9]
colnames(mnBiDF) <- c("1 Glutamatergic Neurons", "2 Intermediate Progenitors"
  , "3 Radial Glia", "4 Astroglia", "5 Ventral Progenitors"
  ,  "6 GABAergic Neurons", "7 OPC", "8 Choroid Plexus")

table(row.names(mnDF) %in% row.names(mnBiDF))

df <- merge(mnDF, mnBiDF, by = "row.names")

corDF <- cor(df[ ,-1], method = "spearman")

round(corDF, 2)


## Geschwind and Pasca normalizations

# Geschwind
exDF <- as.matrix(seuratO@scale.data)
df <- data.frame(seuratO@ident, t(exDF))
mnDF <- sapply(split(df[,-1], df$seuratO.ident), function(df) {rank(colMeans(df))})

# Pasca
df <- data.frame(biMtDF["X9.final.clusters..one.empty."], biNmExDF)
mnBiDF <- sapply(split(df[,-1], df$X9.final.clusters..one.empty.), function(df) {rank(colMeans(df))})
mnBiDF <- mnBiDF[ ,-9]
colnames(mnBiDF) <- c("1 Glutamatergic Neurons", "2 Intermediate Progenitors"
  , "3 Radial Glia", "4 Astroglia", "5 Ventral Progenitors"
  ,  "6 GABAergic Neurons", "7 OPC", "8 Choroid Plexus")

table(row.names(mnDF) %in% row.names(mnBiDF))

df <- merge(mnDF, mnBiDF, by = "row.names")

corDF <- cor(df[ ,-1], method = "spearman")

round(corDF, 2)


## Read depth normalizing log2 transform Geschwind and Pasca raw counts using seurat

dgSO <- new("seurat", raw.data = as.matrix(seuratO@raw.data))

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes 
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules
# (as in Macosko et al. Cell 2015)
dgSO <- Setup(dgSO, min.cells = 3, min.genes = 200, do.logNormalize = T
  , total.expr = 1e4, project = "10X_fetb")


df <- t(biExDF)
colnames(df) <- paste0("Cell_", df[1, ])
df <- df[-1, ]
biSO <- new("seurat", raw.data = df)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes 
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules
# (as in Macosko et al. Cell 2015)
biSO <- Setup(biSO, min.cells = 3, min.genes = 200, do.logNormalize = T
  , total.expr = 1e4, project = "10X_fetb")

# Geschwind
exDF <- as.matrix(dgSO@raw.data)
exDF <- exDF[ ,colnames(exDF) %in% colnames(as.matrix(seuratO@scale.data))]
df <- data.frame(seuratO@ident, t(exDF))
mnDF <- sapply(split(df[,-1], df$seuratO.ident), function(df) {colMeans(df)})

# Pasca
biExDF <- as.matrix(biSO@raw.data)
biExDF <- biExDF[ ,colnames(biExDF) %in% colnames(as.matrix(seuratO@scale.data))]
df <- data.frame(biMtDF["X9.final.clusters..one.empty."], biExDF)
mnBiDF <- sapply(split(df[,-1], df$X9.final.clusters..one.empty.), function(df) {colMeans(df)})
mnBiDF <- mnBiDF[ ,-9]
colnames(mnBiDF) <- paste0("Pasca ", c("1 Glutamatergic Neurons", "2 Intermediate Progenitors"
  , "3 Radial Glia", "4 Astroglia", "5 Ventral Progenitors"
  ,  "6 GABAergic Neurons", "7 OPC", "8 Choroid Plexus"))

table(row.names(mnDF) %in% row.names(mnBiDF))

df <- merge(mnDF, mnBiDF, by = "row.names")

corDF <- cor(df[ ,-1], method = "spearman")

round(corDF, 2)

write.csv(round(corDF, 2), "../analysis/tables/Correlate_Seurat_Birey.csv", quote = FALSE)


