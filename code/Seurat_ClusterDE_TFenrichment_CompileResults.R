# Damon Polioudakis
# 2017-10-03
# Compile results from Vijays TF binding site enrichment pipeline
################################################################################

rm(list = ls())

require(xlsx)
source('GeneSymbol2EnsgID.R')
################################################################################

### Compile list TFs in JASPAR database and their gene symbol aliases

# TFs in JASPAR
fileName <- "../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/TF_Results_Seurat_ClusterDE_TFenrichment_logFC04_Cluster0/JASPAR_CLOVER.txt"
conn <- file(fileName, open = "r")
jaspar <- readLines(conn)
close(conn)
# 204
length(grep(">", jaspar))

v1 <- jaspar[grep(">", jaspar)]
v1 <- gsub(">.*_", "", v1)
v1 <- toupper(v1)
# Split genes of the form:
# SMAD2::SMAD3::SMAD4
v1 <- unlist(strsplit(v1, "::"))

# Look up other TF gene symbol aliases
jasperAliasDF <- FetchGeneAliases(v1)

# Save as csv
write.csv(jasperAliasDF, "../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/Compiled_TFs_in_JASPAR_Database.csv"
  , quote = FALSE)
################################################################################

### Compile list TFs in TRANSFAC database and their gene symbol aliases

# TFs in TRANSFAC
fileName <- fileName <- "~/R/x86_64-pc-linux-gnu-library/3.4/RegFacEnc/data/Clover_TRANSFAC"
conn <- file(fileName, open = "r")
transfac <- readLines(conn)
close(conn)
# 2208
length(transfac[grep(">", transfac)])

v1 <- transfac[grep(">", transfac)]
v1 <- gsub(">.*_._", "", v1)
v1 <- gsub("_.*", "", v1)

v1 <- gsub("ALPHA", "A", v1)
v1 <- gsub("BETA", "B", v1)
v1 <- gsub("GAMMA", "G", v1)
v1 <- gsub("KAPPA", "K", v1)

# Split genes of the form:
# SMAD2::SMAD3::SMAD4
v1 <- unlist(strsplit(v1, "::"))

# Look up other TF gene symbol aliases
transfacAliasDF <- FetchGeneAliases(v1)

# Save as csv
write.csv(transfacAliasDF, "../analysis/Seurat_ClusterDE_TFenrichment/TRANSFAC/Compiled_TFs_in_TRANSFAC_Database.csv"
  , quote = FALSE)
################################################################################

### Compile JASPAR results

# List of file paths
inFiles <- list.files("../analysis/Seurat_ClusterDE_TFenrichment/JASPAR", full.names = TRUE)
inFiles <- inFiles[grep("TF_Results_Seurat_ClusterDE", inFiles)]
inFiles <- list.files(inFiles, full.names = TRUE)
inFiles <- inFiles[grep("xls", inFiles)]

# Load Vijay TF output into list of data frames
ldf <- lapply(inFiles, function(inFile) read.table(inFile, skip = 1))
# Load header
header <- read.table(inFiles[1], nrow = 1)

# Add cluster name to list
names <- basename(inFiles)
names <- gsub("TF_Enrichment_Results Seurat_ClusterDE_TFenrichment_logFC04_Cluster", "", names)
names <- gsub(" .xls", "", names)
names(ldf) <- names
# Order
nameOrder <- as.character(sort(as.numeric(names)))
ldf <- ldf[nameOrder]

# Format and combine data frames into one
ldf1 <- lapply(names(ldf), function(name) {
  df <- data.frame(ldf[[name]])
  df <- df[ ,-1]
  colnames(df) <- unlist(header[1, ])
  
  # Split RXRA::VDR to separate rows
  df$GENE <- df$GENENAME
  doubleRows <- df[grep("\\:\\:", df$GENENAME), ]
  if (nrow(doubleRows) > 0) {
    
    doublesDF <- apply(doubleRows, 1, function(row){
      # Split genes of the form:
      # SMAD2::SMAD3::SMAD4
      v1 <- unlist(strsplit(row[["GENE"]], "::"))
      df <- data.frame(lapply(v1, function(gene) c(row, gene)))
      df <- as.data.frame(t(df))
      df$GENE <- v1
      row.names(df) <- v1
      df <- df[ ,! colnames(df) %in% "V8"]
      return(df)
      
    })
    
    doublesDF <- data.frame(doublesDF[[1]])
    doublesDF$GENE <- as.character(doublesDF$GENE)
    df$GENE <- as.character(df$GENE)
    
    colnames(doublesDF) <- colnames(df)
    
    df <- rbind(df, doublesDF)
  }
  
  # Add column with cluster name
  df$CLUSTER <- name
  
  # Add gene symbol aliases
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  # Look up other TF gene symbol aliases
  aliasDF <- FetchGeneAliases(df$GENE)
  # Duplicate rows that have multiple aliases
  df <- df[match(aliasDF$symbol, df$GENE), ]
  df$ALIAS <- aliasDF$alias
  
  return(df)
})
# Combine into data frame
df <- do.call("rbind", ldf1)
df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)

# Save as csv
write.csv(df, "../analysis/Seurat_ClusterDE_TFenrichment/JASPAR/Seurat_ClusterDE_TFenrichment_logFC04_TF_Enrichment_Results.csv")
################################################################################

### Compile TRANSFAC results

# List of file paths
inFiles <- list.files("../analysis/Seurat_ClusterDE_TFenrichment/TRANSFAC", full.names = TRUE)
inFiles <- inFiles[grep("TF_Results_Seurat_ClusterDE", inFiles)]
inFiles <- list.files(inFiles, full.names = TRUE)
inFiles <- inFiles[grep("xls", inFiles)]

# Load Vijay TF output into list of data frames
ldf <- lapply(inFiles, function(inFile) read.table(inFile, skip = 1))
# Load header
header <- read.table(inFiles[1], nrow = 1)

# Add cluster name to list
names <- basename(inFiles)
names <- gsub("TF_Enrichment_Results Seurat_ClusterDE_TFenrichment_logFC04_Cluster", "", names)
names <- gsub(" .xls", "", names)
names(ldf) <- names
# Order
nameOrder <- as.character(sort(as.numeric(names)))
ldf <- ldf[nameOrder]

# Format and combine data frames into one
ldf1 <- lapply(names(ldf), function(name) {
  df <- data.frame(ldf[[name]])
  df <- df[ ,-1]
  colnames(df) <- unlist(header[1, ])
  
  # Split RXRA::VDR to separate rows
  df$GENE <- df$Genes
  doubleRows <- df[grep("\\:\\:", df$Genes), ]
  if (nrow(doubleRows) > 0) {
    
    doublesDF <- apply(doubleRows, 1, function(row){
      # Split genes of the form:
      # SMAD2::SMAD3::SMAD4
      print(row)
      v1 <- unlist(strsplit(row[["Genes"]], "::"))
      df <- data.frame(lapply(v1, function(gene) c(row, gene)))
      df <- as.data.frame(t(df))
      df$GENE <- v1
      row.names(df) <- v1
      df <- df[ ,! colnames(df) %in% "V8"]
      return(df)
    })
    
    doublesDF <- data.frame(doublesDF[[1]])
    doublesDF$GENE <- as.character(doublesDF$GENE)
    df$GENE <- as.character(df$GENE)
    
    colnames(doublesDF) <- colnames(df)
    
    df <- rbind(df, doublesDF)
  }
  
  # Split HLH-2:HLH-15 to separate rows
  df$GENE <- df$Genes
  doubleRows <- df[grep("\\:", df$Genes), ]
  if (nrow(doubleRows) > 0) {
    
    doublesDF <- apply(doubleRows, 1, function(row){
      # Split genes of the form:
      # HLH-2:HLH-15
      print(row)
      v1 <- unlist(strsplit(row[["Genes"]], ":"))
      df <- data.frame(lapply(v1, function(gene) c(row, gene)))
      df <- as.data.frame(t(df))
      df$GENE <- v1
      # Deal with TFs of the form P50:P50 (can't have duplicate rownames)
      if (any(duplicated(v1))) {
        v2 <- paste0(v1, seq(1, length(v1)))
        row.names(df) <- v2
      }
      else {
        row.names(df) <- v1  
      }
      df <- df[ ,! colnames(df) %in% "V8"]
      return(df)
    })
    
    doublesDF <- data.frame(doublesDF[[1]])
    doublesDF$GENE <- as.character(doublesDF$GENE)
    df$GENE <- as.character(df$GENE)
    
    colnames(doublesDF) <- colnames(df)
    
    df <- rbind(df, doublesDF)
  }
  
  # Add column with cluster name
  df$CLUSTER <- name
  
  # Add gene symbol aliases
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  # Look up other TF gene symbol aliases
  aliasDF <- FetchGeneAliases(df$GENE)
  # Duplicate rows that have multiple aliases
  df <- df[match(aliasDF$symbol, df$GENE), ]
  df$ALIAS <- aliasDF$alias
  
  return(df)
})
# Combine into data frame
df <- do.call("rbind", ldf1)
df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

# Save as csv
write.csv(df, "../analysis/Seurat_ClusterDE_TFenrichment/TRANSFAC/Seurat_ClusterDE_TFenrichment_logFC04_TF_Enrichment_Results.csv")
################################################################################