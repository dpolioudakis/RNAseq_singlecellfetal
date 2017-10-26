# Damon Polioudakis
# 2017-06-20
# Calculate read depth per cell
##################################################################################

rm(list=ls())


# Input directory paths
ssDF <- read.csv("../analysis/tables/QC_RNAstar_Stats_HsMm.csv", header = TRUE
  , stringsAsFactors = FALSE)

nCells <- c(
  # Brain 2
  N701 = 800
  , N702 = 800
  , N703 = 800
  , N704 = 800
  , N705 = 800
  , N706 = 800
  , N707 = 800
  , N708 = 1600
  , N709 = 1600
  # Brain 3
  , N701 = 2500
  , N702 = 2500
  , N704 = 2000
  , N705 = 2000
  , N710 = 1500
  , N711 = 1500
  # Brain 4 and 5
  , N712 = 2900
  , N714 = 3100
  , N715 = 2625
  , N716 = 2750
  , N722 = 2200
  , N723 = 2800
  , N724 = 2375
  , N726 = 2625
)

# Split some of brain 2 cells to aggregate separately because have same indeces
# as brain 3
ss1DF <- ssDF[1:5, ]
ss2DF <- ssDF[6:53, ]

rdDF <- data.frame(READ_DEPTH = c(tapply(ss1DF$Number.of.input.reads, ss1DF$SampleID, sum)
  , tapply(ss2DF$Number.of.input.reads, ss2DF$SampleID, sum))
, NUMBER_OF_CELLS = nCells, INDEX = names(nCells)
)

rdDF$BRAIN <- c(rep(2, 9), rep(3, 6), 4, 4, 5, 5, 4, 4, 5, 5)

rdDF$READS_PER_CELL <- round(rdDF$READ_DEPTH / rdDF$NUMBER_OF_CELLS, 1)

write.csv(rdDF, "../analysis/tables/Read_Depth_Per_Cell.csv", quote = FALSE
  , row.names = FALSE)