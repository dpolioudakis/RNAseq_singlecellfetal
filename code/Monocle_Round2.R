# Damon Polioudakis
# 2017-08-15
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

# Monocle
load("../analysis/Monocle/Monocle_PC1-40_monocleO.Robj")
# Seurat clustering
load("../analysis/Seurat_Cluster_DS2-11/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Seurat_Cluster_DS2-11_seuratO.Robj")

## Variables
graphCodeTitle <- "Monocle_Round2.R"
outGraph <- "../analysis/graphs/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Monocle_Round2_Comp1-20_"
outTable <- "../analysis/tables/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Monocle_Round2_Comp1-20_"
outRdat <- "../analysis/Monocle_Round2/FtMm250_200-3sdgd_Mt5_RegNumiLibBrain_KeepCC_PC1to40/Monocle_Round2_Comp1-20_"

## Output Directories
outDir <- dirname(outGraph)
dir.create(outDir, recursive = TRUE)
outTableDir <- dirname(outTable)
dir.create(outTableDir, recursive = TRUE)
outRdatDir <- dirname(outRdat)
dir.create(outRdatDir, recursive = TRUE)

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

### Run Monocle

# Extract Seurat cluster IDs from Seurat object then remove object for memory
clusterIDs <- centSO@ident
centSO@raw.data <- NULL
centSO@scale.data <- NULL
# Add custom cluster ID combinations to ID list
clustersL <- c(list(c(0, 1, 4, 13), c(0, 1), c(3, 14), c(5, 6), c(7, 9))
  , as.list(sort(as.numeric(as.character(unique(clusterIDs))))))
# Select ID from list using qsub task ID
clid <- unlist(clustersL[as.numeric(args[1])])

print(paste0("Running Monocle for Seurat Cluster ID: "
  , paste(clid, collapse = " ")))

# Skip clusters < 100 cells
if (sum(clusterIDs %in% clid) < 100) {
  print("Less than 100 cells in cluster, skipping...")
  NULL
  
} else {
  
  # Subset mo object to cells in cluster
  mo <- mo[ ,row.names(pData(mo))[pData(mo)$cluster %in% c(clid)]]
  
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
  png(paste0(outGraph, "OrderingGenesDispersion_Cluster"
    , paste0(clid, collapse = '-'), ".png"))
  plot_ordering_genes(mo_filtered)
  dev.off()
  
  print("Variance explained by each PC")
  # Variance explained by each PC
  png(paste0(outGraph, "PCA_VarianceExplained_Cluster"
    , paste0(clid, collapse = '-'), ".png"))
  plot_pc_variance_explained(mo_filtered, verbose = TRUE
    , use_existing_pc_variance = TRUE, return_all = FALSE) 
  dev.off()
  
  print("Reduce data dimensionality")
  # Reduce data dimensionality
  # Use number of genes expressed or total mRNAs?
  mo_filtered <- reduceDimension(mo_filtered, max_components = 20,
    residualModelFormulaStr = "~individual + librarylab + Total_mRNAs"
    , verbose = TRUE)
  
  # mo_filtered <- reduceDimension(mo_filtered, max_components = 16,
  #   residualModelFormulaStr = "~method + batch + num_genes_expressed"
  #   , verbose = TRUE)
  
  print("Order cells along trajectory")
  # Order cells along trajectory
  mo_filtered <- orderCells(mo_filtered)

  # Save Monocle round 2 objects
  save(mo_filtered, file = paste0(
    outRdat, "monocleO_cluster", paste0(clid, collapse = '-'), ".Robj"))

  
  ## Differential gene test as function of pseudotime
  
  print("Differential gene test as function of pseudotime")
  
  ptDE <- differentialGeneTest(mo_filtered
    , fullModelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
    , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")
  save(ptDE, file = paste0(
    outRdat, "DEpseudotime_cluster", paste0(clid, collapse = '-'), ".Robj"))
  
  # ## DE as function of state using Monocle DE function
  # 
  # print("Differential gene test as function of state (monocle DE)")
  # 
  # # Print states
  # print(unique(pData(mo_filtered)$State))
  # 
  # # Differential gene expression test by state
  # stateDE <- differentialGeneTest(mo_filtered
  #   , fullModelFormulaStr = "~State + individual + librarylab + Total_mRNAs"
  #   , reducedModelFormulaStr = "~individual + librarylab + Total_mRNAs")
  # save(stateDE, file = paste0(
  #   outRdat, "DEstate_cluster", paste0(clid, collapse = '-'), ".Robj"))
  
  
  ## DE as function of state using linear model
  
  Filter_ExpMatrix_PercentExpressedCluster <- function(
    exprM
    , minPercent
    , cellID) {
    
    # Subset expression matrix to cluster
    cdf <- exprM[ ,colnames(exprM) %in% cellID]
    
    # Expressed > 0 counts in > X% of cells in cluster
    idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (minPercent / 100)
    print(paste0("Genes expressed in > ", minPercent, "% of cells in cluster"))
    print(table(idxp))
    
    pfCellIDs <- row.names(cdf)[idxp]
    
    return(pfCellIDs)
  }
  
  Filter_ExpMatrix_FoldChangeCluster <- function(
    exprM
    , foldChange
    , cellID) {
    
    # Subset expression matrix to cluster
    cdf <- exprM[ ,colnames(exprM) %in% cellID]
    # Subset expression matrix to all other cells
    ndf <- exprM[ , ! colnames(exprM) %in% cellID]
    
    # Fold change
    v1 <- rowMeans(cdf) - rowMeans(ndf)
    idxf <- v1 > foldChange
    print(paste0("Genes > ", foldChange, " fold change in cluster versus all other cells"))
    print(table(idxf))
    
    pfCellIDs <- row.names(cdf)[idxf]
    
    return(pfCellIDs)
  }
  
  print("Differential gene test as function of state (linear model DE)")
  
  # Print states
  print(unique(pData(mo_filtered)$State))
  
  # Monocle states
  states <- sort(unique(pData(mo_filtered)$State))
  
  # Loop through states
  stateLDF <- lapply(states, function(state) {
    print(state)
    
    # Cell IDs in state
    cellIDs <- pData(mo_filtered)$CELL[pData(mo_filtered)$State %in% state]
    
    if (length(cellIDs) < 10) {
      print(paste0("State ", state, " is < 10 cells, skipping..."))
      NULL
    
    } else {
      print(paste0("State ", state, " DE is running..."))
    
      # Filter genes expressed in < 10% of cells in state
      genesPt <- Filter_ExpMatrix_PercentExpressedCluster(
        exprM = as.matrix(exprs(mo_filtered))
        , minPercent = 10
        , cellID = cellIDs)
      
      # Filter genes expressed in < 0.2 log fold change of cells
      genesFC <- Filter_ExpMatrix_FoldChangeCluster(
        exprM = noCentExM
        , foldChange = 0.2
        , cellID = cellIDs)
      
      # Subset expression matrix to genes
      exM <- as.matrix(centSO@data)
      exM <- exM[
        row.names(exM) %in% c(genesPt) &
          row.names(exM) %in% c(genesFC)
        , ]
      # And cells in monocle object
      exM <- exM[ ,colnames(exM) %in% colnames(exprs(mo_filtered))]
      
      # DE Linear model
      termsDF <- pData(mo_filtered)[c("nUMI", "librarylab", "individual")]
      # Add term TRUE/FALSE cell is in state
      termsDF$state <- FALSE
      termsDF$state[pData(mo_filtered)$State == state] <- TRUE
      mod <- "y ~ state+nUMI+librarylab+individual"
      deLM <- DE_Linear_Model(exDatDF = exM, termsDF = termsDF, mod = mod)
      
      # Format LM output into data frame
      # Combine log2 fold changes, p-values
      deDF <- data.frame(GENE = row.names(deLM$coefmat)
        , LOG_FC = deLM$coefmat[ ,2]
        , PVALUE = deLM$pvalmat[ ,2])
      # Order by pvalue
      deDF <- deDF[order(deDF$PVALUE), ]
      # Add cluster ID
      deDF$STATE <- state
      print(head(deDF))
      # Order by log fold change
      deDF <- deDF[order(-deDF$LOG_FC), ]
      
      # FDR correct
      # NOTE: p-values are so low that FDR tool is returning FDR of 1 for everything
      corrected <- fdrtool(deDF$PVALUE, statistic = "pvalue", plot = FALSE)
      deDF$FDR <- corrected$lfdr
      # Check
      table(deDF$PVALUE < 0.05)
      table(deDF$FDR < 0.05)
      print(head(deDF))  
  
      return(deDF)
    }  
  })
  # Remove NULLs from states < 10 cells and combine list of data frames
  stateLDF <- stateLDF[! is.null(stateLDF)]
  stateDE <- do.call("rbind", stateLDF)
  
  # Save DE data frame
  save(stateDE, file = paste0(
    outRdat, "DEstate_cluster", paste0(clid, collapse = '-'), ".Robj"))
  
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
        , fullModelFormulaStr = "~sm.ns(Pseudotime)*Branch + individual + librarylab + Total_mRNAs"
        , reducedModelFormulaStr = "~sm.ns(Pseudotime) + individual + librarylab + Total_mRNAs"
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
    outRdat, "DEbranch_cluster", paste0(clid, collapse = '-'), ".Robj"))
  
}

print("End of Monocle_Round2.R...")
################################################################################