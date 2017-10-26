
#monocle excluding endothelial, OPCs, interneurons

library(monocle)
library(ksheu.lib)
memory.limit()
memory.limit(20000)
source("E:/Rotation3/plot_cell_trajectory2.R")
source("C:/Users/msheu/Desktop/Katherine/RRHO/intersect_then_rank_RRHO.R")

setwd("E:/Rotation3/fromDamon/")
load("DS002003_exon_FtMm250_Seurat_tSNErot.Rdat")
tsneRot$type = pData(HSMM_filtered)$State[match(rownames(tsneRot), pData(HSMM_filtered)$name)]
tsneRot$seurat = pData(HSMM_filtered)$cluster_name[match(rownames(tsneRot), pData(HSMM_filtered)$name)]
p <- ggplot(tsneRot, aes(x = tSNE_1, y = tSNE_2, color = seurat)) + 
  geom_point()
p

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

library(reshape)
tbr1 = data.frame(t(exM[rownames(exM) =="TBR1", ]))
tbr1$type = clusters$CLUSTER_NAME[match(rownames(tbr1), clusters$CELL_ID)]
library(ggplot2)
p <- ggplot(tbr1, aes(TBR1, type))+
  geom_boxplot() + 
  geom_jitter() 
p

DLX2 = data.frame(t(exM[rownames(exM) =="DLX2", ]))
DLX2$type = clusters$CLUSTER_NAME[match(rownames(DLX2), clusters$CELL_ID)]
library(ggplot2)
p <- ggplot(DLX2, aes(DLX2, type))+
  geom_boxplot() + 
  geom_jitter() 
p


########################
#monocle2 continue
########################

subset = clusters$CELL_ID[!grepl("Endothelial|Oligo|Interneuron|Deep|oRG", clusters$CLUSTER_NAME)]
#subset = clusters$CELL_ID[grepl("Interneuron", clusters$CLUSTER_NAME)]
exM = exM[, colnames(exM) %in% subset]
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

# We can filter genes based on average expression level, and we can additionally
# select genes that are unusually variable across cells.
disp_table <- dispersionTable(HSMM) 
ordering_genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 0.25 * dispersion_fit)$gene_id
# The setOrderingFilter function marks genes that will be used for clustering in
# subsequent calls to clusterCells
HSMM_filtered <- setOrderingFilter(HSMM, ordering_genes)
# The plot_ordering_genes function shows how variability (dispersion) in a
# gene's expression depends on the average expression across cells. The red line
# shows Monocle's expectation of the dispersion based on this relationship. The
# genes we marked for use in clustering are shown as black dots, while the
# others are shown as grey dots.
plot_ordering_genes(HSMM_filtered)

plot_pc_variance_explained(HSMM_filtered, verbose = T, use_existing_pc_variance = T, return_all = F) 

# Monocle allows us to subtract the effects of "uninteresting" sources of
# variation to reduce their impact on the clustering.
HSMM_filtered <- reduceDimension(HSMM_filtered, max_components = 16
  , residualModelFormulaStr = "~method + batch + num_genes_expressed"
  , verbose = TRUE)




HSMM_filtered <- orderCells(HSMM_filtered)
plot_cell_trajectory(HSMM_filtered, 1, 2, color_by="region")

pData(HSMM_filtered)$cluster_name = clusters$CLUSTER_NAME[match(pData(HSMM_filtered)$name, clusters$CELL_ID)]
head(pData(HSMM_filtered))
plot_cell_trajectory(HSMM_filtered, 1, 4, show_tree = T,show_backbone = T, color_by="cluster_name")


GM_state <- function(cds){ if (length(unique(pData(cds)$State)) > 1){ 
  T0_counts <- table(pData(cds)$State, pData(cds)$cluster_name)[,"Mitotic Progenitors"] 
  return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) 
}else { return (1) } 
} 
HSMM_filtered2 <- orderCells(HSMM_filtered, root_state=GM_state(HSMM_filtered))


environment(plot_cell_trajectory2) <- environment(plot_cell_trajectory)
plot_cell_trajectory(HSMM_filtered2, 1, 3, color_by="Pseudotime" )
plot_cell_trajectory(HSMM_filtered, 1,4,color_by="State")
plot_cell_trajectory(HSMM_filtered, 1,3,color_by="batch")
plot_cell_trajectory(HSMM_filtered, 1,2,color_by="method")

#plot specific genes
pData(HSMM_filtered)$EOMES = log(t(exM["EOMES",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 3,show_tree = T,show_backbone = T, color = "red", color_by="EOMES" )
pData(HSMM_filtered)$SOX2 = log(t(exM["SOX2",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 3,show_tree = T,show_backbone = T, color = "blue", color_by="SOX2" )
pData(HSMM_filtered)$VIM = log(t(exM["VIM",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "blue", color_by="VIM" )
pData(HSMM_filtered)$TBR1 = log(t(exM["TBR1",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 3,show_tree = T,show_backbone = T, color = "darkgreen", color_by="TBR1" )
pData(HSMM_filtered)$TBR1 = t(exM["TBR1",])
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "darkgreen", color_by="TBR1" )
pData(HSMM_filtered)$PAX6 = log(t(exM["PAX6",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "orange", color_by="PAX6" )
pData(HSMM_filtered)$OLIG2 = log(t(exM["OLIG2",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "red", color_by="OLIG2" )
pData(HSMM_filtered)$DLX2 = log(t(exM["DLX2",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "blue", color_by="DLX2" )
pData(HSMM_filtered)$DLX1 = log(t(exM["DLX1",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "blue", color_by="DLX1" )
pData(HSMM_filtered)$OLIG1 = log(t(exM["OLIG1",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "red", color_by="OLIG1" )
pData(HSMM_filtered)$NKX2.1 = log(t(exM["TTF1",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "red", color_by="NKX2.1" )
pData(HSMM_filtered)$CXCR2 = log(t(exM["CXCR2",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "red", color_by="CXCR2" )
pData(HSMM_filtered)$MCM5 = log(t(exM["MCM5",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color = "red", color_by="MCM5" )

neuron_markers = markers$Gene.Symbol[markers$Notes=="Neuron"]
IP_markers = markers$Gene.Symbol[markers$Notes =="IP"]
temp = colSums(exM[neuron_markers,])
pData(HSMM_filtered)$neuron_marker = log(temp + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color_by="neuron_marker")
temp = colSums(exM[IP_markers,])
pData(HSMM_filtered)$IP_marker = log(temp + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color_by="IP_marker" )


pData(HSMM_filtered)$mitotic1_marker = log(colSums(exM[mitotic1_markers,]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color_by="mitotic1_marker" )
pData(HSMM_filtered)$mitotic2_marker = log(colSums(exM[mitotic2_markers,]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 4,show_tree = T,show_backbone = T, color_by="mitotic2_marker" )


#marker genes
marker_genes <- na.omit(markers$Gene.Symbol)
marker_genes = as.character(marker_genes[!duplicated(as.character(marker_genes))])
marker_genes = marker_genes[marker_genes != ""]

#all expressed genes
head(fData(HSMM_filtered))
HSMM_expressed_genes <- row.names(subset(fData(HSMM_filtered), num_cells_expressed >= 10)) #14449 genes
HSMM_expressed <- HSMM_filtered2[HSMM_expressed_genes,]

to_test = row.names(subset(fData(HSMM_expressed), gene_short_name %in% marker_genes))
cds_subset <- HSMM_expressed[to_test,]
head(fData(cds_subset))
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr="~sm.ns(Pseudotime)")
#diff_test_res[,c("gene_short_name", "pval", "qval")]

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.001))
dev.off()
plot_pseudotime_heatmap(HSMM_filtered[sig_gene_names,], cluster_rows = T, #num_clusters = 5, 
                        cores = 1, show_rownames = T)


#heatmap of marker genes by monocle stage
info = pData(HSMM_filtered)
exM_topgenes = exM[rownames(exM) %in% marker_genes, ]
exM_topgenes_log = log2(exM_topgenes+1)
library(pheatmap)
pheatmap(as.matrix(exM_topgenes_log), scale = "row", annotation_col = info[, c(2,9)], cluster_rows = T, 
         clustering_method = "complete", show_colnames = F)



#all expressed genes
head(fData(HSMM_filtered))
HSMM_expressed_genes <- row.names(subset(fData(HSMM_filtered), num_cells_expressed >= 10)) #14449 genes
HSMM_expressed <- HSMM_filtered[HSMM_expressed_genes,]

cds_subset <- HSMM_expressed
head(fData(cds_subset))
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr="~sm.ns(Pseudotime)")
#diff_test_res[,c("gene_short_name", "pval", "qval")]

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.000001))
dev.off()
plot_pseudotime_heatmap(HSMM_filtered[sig_gene_names,], cluster_rows = T, #num_clusters = 9, 
                        cores = 1, show_rownames = T)


# saveRDS(HSMM_filtered, "monocle_16D_0.5_neuronal_lineage.rds")
neuronal = readRDS("monocle_16D_0.5_neuronal_lineage.rds")
plot_cell_trajectory(neuronal, 1, 2, show_tree = T,show_backbone = T, color_by="cluster_name")


# saveRDS(HSMM_filtered, "monocle_16D_0.5_all.rds")
HSMM_filtered = readRDS("monocle_16D_0.5_all.rds")

# saveRDS(HSMM_filtered2, "monocle_16D_0.5_mito_IP_excit_upper12.rds")

layers = read.csv("Bakken_2016_ST7_LayerMarkers.csv")
pData(HSMM_filtered)$MYCN = log(t(exM["MYCN",]) + 1, 2)
plot_cell_trajectory2(HSMM_filtered, 1, 2,show_tree = T,show_backbone = T, color = "red", color_by="MYCN" )

