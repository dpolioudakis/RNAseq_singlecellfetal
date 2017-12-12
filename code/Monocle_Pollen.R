# Damon Polioudakis
# 2017-12-05
# Run monocle on each Seurat cluster
################################################################################

rm(list = ls())

require(monocle)
require(Seurat)
require(cowplot)
require(ggplot2)
require(viridis)
require(reshape2)
require(fdrtool)
source("Function_Library.R")

# qsub task ID for selecting Seurat cluster ID to run through Monocle
args <- commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(args)

## Inputs

exDF <- read.csv("../pollen_2015/data/htseq/Exprs_HTSCexon.csv")
mtDF <- read.csv("../pollen_2015/metadata/Cell paper - updated attributes.csv", header = TRUE)

## Variables
graphCodeTitle <- "Monocle_Pollen.R"
outGraph <- "../analysis/graphs/Monocle_Pollen/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/Monocle_Pollen_"
outTable <- "../analysis/tables/Monocle_Pollen/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/Monocle_Pollen_"
outRdat <- "../analysis/Monocle_Pollen/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Comp1-10/Monocle_Pollen_"

## Output Directories
dir.create(dirname(outGraph), recursive = TRUE)
dir.create(dirname(outTable), recursive = TRUE)
dir.create(dirname(outRdat), recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 12)))
theme_update(plot.title = element_text(size = 12))
theme_update(axis.line = element_line(colour = "black")
  , plot.background = element_blank()
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , panel.border = element_blank()
)
################################################################################

### Format

# Remove ERCCs and STAR stats from tail
tail(exDF, 10)[1:10, 1:5]
exDF <- head(exDF, -97)
tail(exDF, 5)[1:10, 1:5]

# Move gene names
row.names(exDF) <- exDF$X
exDF <- exDF[ ,-1]

# Clean chr number of ens IDs
row.names(exDF) <- gsub("\\.[0-9]*", "", row.names(exDF))

# Format metadata cell IDs to match expression matrix column names
mtDF$Cell <- gsub("-", ".", mtDF$Cell)

# Subset metadata to cells in expression data
mtDF <- mtDF[as.character(mtDF$Cell) %in% colnames(exDF), ]

# Metadata row names to cell IDs
row.names(mtDF) <- mtDF$Cell

# Order metadata as columns in expression data
mtDF <- mtDF[match(colnames(exDF), row.names(mtDF)), ]
################################################################################

### Run Monocle

cell_type_combos <- list(
  c("Neuron", "IPC")
  , c("Neuron")
  , c("Neuron", "RG", "IPC")
  , c("Neuron", "Interneuron", "RG", "IPC")
  )

lapply(cell_type_combos, function(cellTypes) {
  
  print(paste0("Running monocle for Pollen cell types: ", cellTypes))
  
  # Subset metadata to cell types of interest
  subMtDF <- mtDF[mtDF$Inferred.Cell.Type %in% cellTypes, ]
  # Subset expression data to cells in metadata
  subExDF <- exDF[ ,colnames(exDF) %in% as.character(subMtDF$Cell)]
  
  # Setup monocle object
  print("Initializing Monocle object...")
  feature_data = data.frame(gene_short_name = rownames(subExDF))
  rownames(feature_data) = feature_data$gene_short_name
  pd <- new("AnnotatedDataFrame", data = subMtDF)
  fd <- new("AnnotatedDataFrame", data = feature_data) 
  # Note metadata row order must match expression data column order
  mo <- newCellDataSet(cellData = as(as.matrix(subExDF), "sparseMatrix"), 
    phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, 
    expressionFamily = negbinomial.size())
  
  # Filter genes and cells
  # retains all 
  mo <- detectGenes(mo, min_expr = 0.1)
  print(head(fData(mo)))
  # genes expressed in at least 10 cells
  expressed_genes <- row.names(subset(fData(mo), num_cells_expressed >= 10))
  print(head(pData(mo)))
  
  # If you are using RPC values to measure expression, as we are in this vignette,
  # it's also good to look at the distribution of mRNA totals across the cells:
  pData(mo)$Total_mRNAs <- Matrix::colSums(exprs(mo))
  
  mo <- mo[,pData(mo)$Total_mRNAs < 1e6]
  
  pdf(paste0(outGraph, "nUMI_density_"
    , paste0(cellTypes, collapse = '-'), ".pdf"))
  qplot(Total_mRNAs, data = pData(mo), color = as.factor(Anatomical.Source), geom = "density")
  qplot(Total_mRNAs, data = pData(mo), color = Brain, geom = "density")
  dev.off()
  
  # Filter genes and cells
  # retains all
  mo <- detectGenes(mo, min_expr = 0.1)
  print(head(fData(mo)))
  # genes expressed in at least 10 cells
  expressed_genes <- row.names(subset(fData(mo), num_cells_expressed >= 10))
  print(head(pData(mo)))
  
  print("Size factors and dispersion")
  # Size factors and dispersion
  mo <- estimateSizeFactors(mo) 
  mo <- estimateDispersions(mo)
  
  print("Dispersed genes to use for pseudotime ordering")
  # Dispersed genes to use for pseudotime ordering
  disp_table <- dispersionTable(mo) 
  ordering_genes <- subset(disp_table
    , mean_expression >= 0.01 & dispersion_empirical >= 0.75 * dispersion_fit)$gene_id
  mo_filtered <- setOrderingFilter(mo, ordering_genes) 
  # Plot
  png(paste0(outGraph, "OrderingGenesDispersion_"
    , paste0(cellTypes, collapse = '-'), ".png"))
  plot_ordering_genes(mo_filtered)
  dev.off()
  
  # print("Variance explained by each PC")
  # # Variance explained by each PC
  # png(paste0(outGraph, "PCA_VarianceExplained_Cluster_"
  #   , paste0(cellTypes, collapse = '-'), ".png"))
  # plot_pc_variance_explained(mo_filtered, verbose = TRUE
  #   , use_existing_pc_variance = TRUE, return_all = FALSE) 
  # dev.off()
  
  print("Reduce data dimensionality")
  # Reduce data dimensionality
  # Use number of genes expressed or total mRNAs?
  mo_filtered <- reduceDimension(mo_filtered, max_components = 10,
    residualModelFormulaStr = "~Total_mRNAs"
    , verbose = TRUE)
  
  print("Order cells along trajectory")
  # Order cells along trajectory
  mo_filtered <- orderCells(mo_filtered)
  
  # Save Monocle round 2 objects
  save(mo_filtered, file = paste0(
    outRdat, "monocleO_", paste0(cellTypes, collapse = '-'), ".Robj"))
  
  
  ## Differential gene test as function of pseudotime
  
  print("Differential gene test as function of pseudotime")
  
  ptDE <- differentialGeneTest(mo_filtered
    , fullModelFormulaStr = "~sm.ns(Pseudotime) + Total_mRNAs"
    , reducedModelFormulaStr = "~Total_mRNAs")
  save(ptDE, file = paste0(
    outRdat, "DEpseudotime_", paste0(cellTypes, collapse = '-'), ".Robj"))
  
  
  # ## DE as function of state using linear model
  # 
  # print("Differential gene test as function of state (linear model DE)")
  # 
  # # Print states
  # print(unique(pData(mo_filtered)$State))
  # 
  # # Monocle states
  # states <- sort(unique(pData(mo_filtered)$State))
  # 
  # # Loop through states
  # stateLDF <- lapply(states, function(state) {
  #   print(state)
  #   
  #   # Cell IDs in state
  #   cellIDs <- pData(mo_filtered)$CELL[pData(mo_filtered)$State %in% state]
  #   
  #   if (length(cellIDs) < 10) {
  #     print(paste0("State ", state, " is < 10 cells, skipping..."))
  #     NULL
  #     
  #   } else {
  #     print(paste0("State ", state, " DE is running..."))
  #     
  #     # Subset expression matrix to genes
  #     exM <- as.matrix(centSO@data)
  #     exM <- exM[
  #       row.names(exM) %in% c(genesPt) &
  #         row.names(exM) %in% c(genesFC)
  #       , ]
  #     # And cells in monocle object
  #     exM <- exM[ ,colnames(exM) %in% colnames(exprs(mo_filtered))]
  #     
  #     # DE Linear model
  #     termsDF <- pData(mo_filtered)[c("nUMI")]
  #     # Add term TRUE/FALSE cell is in state
  #     termsDF$state <- FALSE
  #     termsDF$state[pData(mo_filtered)$State == state] <- TRUE
  #     mod <- "y ~ state+nUMI"
  #     deLM <- DE_Linear_Model(exDatDF = exM, termsDF = termsDF, mod = mod)
  #     
  #     # Format LM output into data frame
  #     # Combine log2 fold changes, p-values
  #     deDF <- data.frame(GENE = row.names(deLM$coefmat)
  #       , LOG_FC = deLM$coefmat[ ,2]
  #       , PVALUE = deLM$pvalmat[ ,2])
  #     # Order by pvalue
  #     deDF <- deDF[order(deDF$PVALUE), ]
  #     # Add cluster ID
  #     deDF$STATE <- state
  #     print(head(deDF))
  #     # Order by log fold change
  #     deDF <- deDF[order(-deDF$LOG_FC), ]
  #     
  #     # FDR correct
  #     # NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
  #     corrected <- fdrtool(deDF$PVALUE, statistic = "pvalue", plot = FALSE)
  #     deDF$FDR <- corrected$lfdr
  #     # Check
  #     table(deDF$PVALUE < 0.05)
  #     table(deDF$FDR < 0.05)
  #     print(head(deDF))  
  #     
  #     return(deDF)
  #   }  
  # })
  # # Remove NULLs from states < 10 cells and combine list of data frames
  # stateLDF <- stateLDF[! is.null(stateLDF)]
  # stateDE <- do.call("rbind", stateLDF)
  # 
  # # Save DE data frame
  # save(stateDE, file = paste0(
  #   outRdat, "DEstate_cluster", paste0(cellTypes, collapse = '-'), ".Robj"))
  # 
  # # Differential gene expression test by state
  # stateDE <- differentialGeneTest(mo_filtered
  #   , fullModelFormulaStr = "~State + individual + librarylab + Total_mRNAs"
  #   , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")
  # save(stateDE, file = paste0(
  #   outRdat, "DEstate_cluster", paste0(clid, collapse = '-'), ".Robj"))
  
  ## DE BEAM (branch DE)
  
  print("Differential branch expression")
  
  # Not sure how to extract states from monocle object so loop through 1-100
  # and return NULL if state does not exist
  states <- c(1:100)
  beamLDF <- lapply(states, function(branch_point){
    BEAM_res <- tryCatch(
      BEAM(
        mo_filtered
        , fullModelFormulaStr = "~sm.ns(Pseudotime)*Branch + Total_mRNAs"
        , reducedModelFormulaStr = "~sm.ns(Pseudotime) + Total_mRNAs"
        , branch_point = branch_point, cores = 1, verbose = TRUE
      )
      , error = function(e) NULL)
    return(BEAM_res)
  })
  # Remove NULLs and change list names to states
  idx <- sapply(beamLDF, is.null)
  beamLDF <- beamLDF[! idx]
  names(beamLDF) <- states[! idx]
  # Save
  save(beamLDF, file = paste0(
    outRdat, "DEbranch_cluster", paste0(cellTypes, collapse = '-'), ".Robj"))
  
})
print("End of Monocle_Round2.R...")