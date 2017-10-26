# Damon Polioudakis
# 2016-12-08
# Recluster the clusters output by Cluster_Seurat.R

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3.0
################################################################################

rm(list = ls())

require(Seurat)

# Cluster identities
load("../analysis/Cluster_Seurat/Cluster_Seurat_exon_FtMm250_fetb_seurat.Robj")
clIdDF <- data.frame(fetb@ident)
################################################################################

### Cells per cluster

# Total cells per cluster
cpcDF <- data.frame(table(fetb@ident))
cpcDF$PERCENTAGE <- round((cpcDF$Freq / sum(cpcDF$Freq)) * 100, 2)

# VZ cells per cluster
clIdDF$CELL_ID <- row.names(clIdDF)
df <- clIdDF[grep("_VZ", row.names(clIdDF)), ]
df <- data.frame(table(df$fetb.ident))
df$PERCENTAGE <- round((df$Freq / sum(df$Freq)) * 100, 2)
cpcDF <- merge(cpcDF, df, by = "Var1")

# CP cells per cluster
clIdDF$CELL_ID <- row.names(clIdDF)
df <- clIdDF[grep("_CP", row.names(clIdDF)), ]
df <- data.frame(table(df$fetb.ident))
df$PERCENTAGE <- round((df$Freq / sum(df$Freq)) * 100, 2)
cpcDF <- merge(cpcDF, df, by = "Var1")

colnames(cpcDF) <- c("CLUSTER", "NUMBER_OF_CELLS", "PERCENT"
  , "NUMBER_OF_CELLS_VZ", "PERCENT_VZ", "NUMBER_OF_CELLS_CP", "PERCENT_CP")

write.table(cpcDF, "../analysis/tables/Cluster_Seurat_Cells_Per_Cluster.txt"
  , sep = "\t", quote = FALSE, row.names = FALSE)
################################################################################
