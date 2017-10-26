library(monocle)
library(ksheu.lib)
memory.limit()
memory.limit(20000)

setwd("D:/Rotation3/")
load("DS002003_exon_FtMm250_Seurat_RawData.Rdat")
#write.table(exM, "DS002003_exon_FtMm250_Seurat_RawData.txt", row.names = T, quote = F, sep = "\t")
#load("/u/project/eeskin/geschwind/dpolioud/RNAseq_singlecellfetal/analysis/DS002003_exon_FtMm250_Seurat_RawData.Rdat")
markers = read.csv("MarkersforSingleCell_2017-01-05.csv")
clusters = read.delim("DS002003_exon_FtMm250_Seurat_ClusterIDs.txt")
assignments = data.frame(current.cluster.ids =  c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),new.cluster.ids = c("Excitatory Upper Layer Neuron 1",
                                                                                                    "Excitatory Neuron",
                                                                                                    "Excitatory Upper Layer Neuron 2",
                                                                                                    "Excitatory Deep Layer Neuron",
                                                                                                    "Intermediate Progenitors", 
                                                                                                    "Interneuron", 
                                                                                                    "Mitotic Progenitors",
                                                                                                    "oRG", 
                                                                                                    "Oligodendrocyte Precursor",
                                                                                                    "Endothelial"))
assignments$current.cluster.ids = as.factor(assignments$current.cluster.ids)
clusters$CLUSTER_ID = as.factor(clusters$CLUSTER_ID)
clusters$CLUSTER_NAME = assignments$new.cluster.ids[match(clusters$CLUSTER_ID, assignments$current.cluster.ids)]   

exM = exM[, colnames(exM) %in% clusters$CELL_ID]
#set.seed(1234)
#exM = exM[, sample(colnames(exM), 2500)]


phenoData = data.frame(name = colnames(exM))
phenoData$region = ifelse(grepl("CP", phenoData$name), "CP", "GZ")
phenoData$method = ifelse(grepl("H", phenoData$name), "hard", "soft")
phenoData$batch = ifelse(grepl("1", phenoData$name), "1", "2")
rownames(phenoData) = phenoData$name

feature_data = data.frame(gene_short_name = rownames(exM))
rownames(feature_data) = feature_data$gene_short_name

pd <- new("AnnotatedDataFrame", data = phenoData) 
fd <- new("AnnotatedDataFrame", data = feature_data) 
HSMM <- newCellDataSet(as(as.matrix(exM), "sparseMatrix"), 
                       phenoData = pd, featureData = fd, lowerDetectionLimit=0.5, 
                       expressionFamily=negbinomial.size())

HSMM <- detectGenes(HSMM, min_expr = 0.1) #retains all
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10)) #genes expressed in at least 10 cells
print(head(pData(HSMM)))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- estimateSizeFactors(HSMM) 
HSMM <- estimateDispersions(HSMM)


disp_table <- dispersionTable(HSMM) 
ordering_genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 0.25 * dispersion_fit)$gene_id
HSMM_filtered <- setOrderingFilter(HSMM, ordering_genes) 
plot_ordering_genes(HSMM_filtered)

plot_pc_variance_explained(HSMM_filtered, verbose = T, use_existing_pc_variance = T, return_all = F) 

HSMM_filtered <- reduceDimension(HSMM_filtered, max_components=25,
                                 residualModelFormulaStr = "~method + batch + num_genes_expressed",verbose = T)
HSMM_filtered <- orderCells(HSMM_filtered)
plot_cell_trajectory(HSMM_filtered, 1, 3, color_by="region")

pData(HSMM_filtered)$cluster_name = clusters$CLUSTER_NAME[match(pData(HSMM_filtered)$name, clusters$CELL_ID)]
head(pData(HSMM_filtered))
plot_cell_trajectory(HSMM_filtered, 1, 3, show_tree = T,show_backbone = T, color_by="cluster_name")

GM_state <- function(cds){ if (length(unique(pData(cds)$State)) > 1){ 
  T0_counts <- table(pData(cds)$State, pData(cds)$cluster_name)[,"Mitotic Progenitors"] 
  return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) 
}else { return (1) } 
} 
HSMM_filtered <- orderCells(HSMM_filtered, root_state=GM_state(HSMM_filtered))
plot_cell_trajectory(HSMM_filtered, 1, 2, color_by="Pseudotime")
plot_cell_trajectory(HSMM_filtered, 1,3,color_by="State")


neuron_markers = markers$Gene.Symbol[markers$Notes=="Neuron"]
IP_markers = markers$Gene.Symbol[markers$Notes =="IP"]

temp = (colMeans(exM[neuron_markers,][exM[neuron_markers,]!=0]) )
pData(HSMM_filtered)$neuron_marker = log(temp + 1, 2)
plot_cell_trajectory(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color_by="neuron_marker")
pData(HSMM_filtered)$IP_marker = log(t(colMeans(exM[IP_markers,])) + 1, 2)
plot_cell_trajectory(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color_by="IP_marker" )

pData(HSMM_filtered)$EOMES = log(t(exM["EOMES",]) + 1, 2)
plot_cell_trajectory(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color_by="EOMES" )

 




plot_cell_trajectory(HSMM_filtered, 1, 3, show_backbone = T, markers = "EOMES",color_by="cluster_name")

#look at marker genes
head(fData(HSMM_filtered))
HSMM_expressed_genes <- row.names(subset(fData(HSMM_filtered), num_cells_expressed >= 10)) 
HSMM_expressed <- HSMM_filtered[HSMM_expressed_genes,]
#my_genes <- row.names(subset(fData(HSMM_expressed), gene_short_name %in% c("ASCL1", "TBR1", "PAX6", "OLIG2")))
#cds_subset <- HSMM_expressed[my_genes,] 
#plot_genes_branched_pseudotime(cds_subset, branch_point=1, color_by="State", ncol=1)

marker_genes <- na.omit(markers$Gene.Symbol)
marker_genes = as.character(marker_genes[!duplicated(as.character(marker_genes))])
marker_genes = marker_genes[marker_genes != ""]


to_test = row.names(subset(fData(HSMM_filtered), gene_short_name %in% marker_genes))
cds_subset <- HSMM_filtered[to_test,]
head(fData(cds_subset))
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr="~sm.ns(Pseudotime)")
#diff_test_res[,c("gene_short_name", "pval", "qval")]

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
dev.off()
plot_pseudotime_heatmap(HSMM_filtered[sig_gene_names,], cluster_rows = T, num_clusters = 9, 
                        cores = 1, show_rownames = T)

library(Vennerable)
overlap = intersect_all(sig_gene_names, marker_genes)
VennSet = Venn(SetNames = c("literature markers", "pseudotime markers"), Weight = c('00' = 0, '01' = 0, '11' = 124, '10' = 90))
VennList <- compute.Venn(VennSet, doWeights = TRUE)
plot(VennList)

#analyzing branches


BEAM_res <- BEAM(cds_subset, branch_point=1, cores = 1) 
BEAM_res <- BEAM_res[order(BEAM_res$qval),] 
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.table(BEAM_res, "BEAM_res_5D_all_branch1", sep = "\t", quote = F, row.names = T)

plot_genes_branched_heatmap(HSMM_filtered[row.names(subset(BEAM_res, qval < 1e-20)),], 
                            branch_point = 1, #num_clusters = 6, 
                            cores = 1, use_gene_short_name = T, show_rownames = T)
test_genes <- row.names(subset(fData(HSMM_filtered), gene_short_name %in% c("ASCL1", "TBR1", "PAX6", "OLIG2"))) 
plot_genes_branched_pseudotime(HSMM_filtered[test_genes,], branch_point=1, color_by="cluster_name", ncol=1)


overlap = intersect_all(BEAM_res$gene_short_name, marker_genes)
VennSet = Venn(SetNames = c("literature markers", "pseudotime markers"), Weight = c('00' = 0, '01' = 0, '11' = 124, '10' = 90))
VennList <- compute.Venn(VennSet, doWeights = TRUE)
plot(VennList)

