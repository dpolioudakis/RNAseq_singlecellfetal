# Damon Polioudakis
# 2017-04-11
# Format RNA STAR stats into table from Compile_RNAstar_Stats_HsMm.sh
################################################################################

rm(list=ls())


# Input directory paths
inSeqDir <- c("DS-002-011", "DS-003", "DS-004", "DS-005-006-007-008"
  , "DS-008-011", "DS-009")
inDirs <- c(
  list.dirs(paste0("../", inSeqDir, "/metadata/human_mouse"))
  # , list.dirs(paste0("../", inSeqDir, "/metadata/GRCh37_75_assembly_NoERCC"))
  , list.dirs(paste0("../", inSeqDir, "/metadata/GRCh38_Gencode25"))
)
starStatsFilesL <- list.files(inDirs, full.names = TRUE)
#
# starStatsFilesL <- list.files("../metadata/human_mouse/")
# starStatsFilesL <- starStatsFilesL[grep("RNAstar_Stats*", starStatsFilesL)]

starStatsL <- list()
for (i in 1:length(starStatsFilesL)) {

  starStatsFile <- starStatsFilesL[[i]]

  starDF <- read.table(starStatsFile, header = TRUE, fill = TRUE, sep = "\t")

  df <- data.frame(t(starDF[ ,c(1, 2, 3, 5, 6, 7, 20, 21, 22, 23, 25, 26, 27)]))

  rNames <- gsub("X", "%", row.names(df))
  rNames <- gsub("\\.", " ", rNames)

  outDF <- cbind(rNames, df)

  lane <- gsub(".*Stats_", "", starStatsFile)
  lane <- gsub(".txt", "", lane)
  lane <- c("Lane", rep(lane, (ncol(outDF)-1)))
  lane <- t(data.frame(lane))
  print(lane)

  outDF <- rbind(lane, as.matrix(outDF))

  starStatsL[[i]] <- outDF

}
starStatsL <- lapply(starStatsL, function(df) df[, -1])
starStatsDF <- do.call(cbind, starStatsL)
starStatsDF <- t(starStatsDF)

dir.create("../analysis/tables", recursive = TRUE)
write.csv(starStatsDF, "../analysis/tables/QC_RNAstar_Stats_HsMm.csv"
  , quote = FALSE, row.names = FALSE)
################################################################################
