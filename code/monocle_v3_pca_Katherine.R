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




#for mitotic progenitors
head(pData(HSMM_filtered))
test = pData(HSMM_filtered)[pData(HSMM_filtered)$cluster_name =="Mitotic Progenitors", ]

exM = exM[, colnames(exM) %in% clusters$CELL_ID]
subset = clusters$CELL_ID[grepl("Mitotic", clusters$CLUSTER_NAME)]
exM = exM[, colnames(exM) %in% subset]

#exM_orderinggenes = exM[rownames(exM) %in% ordering_genes, 

#temp
data = cbind(gene = rownames(data), data)

write.table(data, "Mitotic_singlecell.txt", row.names = F, quote = F, sep = "\t")


PCA_from_file(file = "Mitotic_singlecell.txt", center = T, scale = T)
plot_pca("Mitotic_singlecell_prcomp_scores.txt", test$name, test$State, title = "Mitotic Progenitors", PCx = "PC1", PCy = "PC2",labels =F)

loadings = read.delim("Mitotic_singlecell_prcomp_loadings.txt")
loadings_PC2 = loadings[, c(1,3)]
write.table(loadings_PC2, "Mitotic_singlecell_prcomp_loadings_PC2.rnk",row.names = F, col.names = F, quote = F, sep = "\t")

top_loadings = loadings_PC2[order((loadings_PC2$PC2), decreasing = T),]
top_loadings = top_loadings$Loading[c(1:50, 17193:17242)]
mitotic1_markers = top_loadings$Loading[c(1:50)]
mitotic2_markers = top_loadings$Loading[c(17193:17242)]

exM_topgeens = exM[rownames(exM) %in% top_loadings, colnames(exM) %in% subset]
exM_topgeens_log = log2(exM_topgeens+1)
library(pheatmap)
library(RColorBrewer)
colors<-colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(255)
pheatmap(as.matrix(exM_topgeens_log), col=colors, scale = "row", annotation_col = test[, c(2,9)], cluster_rows = T, 
         clustering_method = "complete", show_colnames = F)

low_expressing = test$name[test$num_genes_expressed <1000]
exM_topgeens_log_highexp = exM_topgeens_log[, !(colnames(exM_topgeens_log) %in% low_expressing)]
pheatmap(as.matrix(exM_topgeens_log_highexp), col=colors, scale = "row", annotation_col = test[, c(2,9)], cluster_rows = T, 
         clustering_method = "complete", show_colnames = F)


#for excitatory neurons
test = pData(HSMM_filtered)[pData(HSMM_filtered)$cluster_name =="Excitatory Neuron", ]
exM = exM[, colnames(exM) %in% clusters$CELL_ID]
exM_excitneu = exM[, colnames(exM) %in% test$name]
exM_excitneu = cbind(gene = rownames(exM_excitneu), exM_excitneu)

write.table(exM_excitneu, "excitatoryneuron_singlecell.txt", row.names = F, quote = F, sep = "\t")
PCA_from_file(file = "excitatoryneuron_singlecell.txt", center = T, scale = T)
plot_pca("excitatoryneuron_singlecell_prcomp_scores.txt", pData(HSMM_filtered)$name, pData(HSMM_filtered)$State, title = "Excitatory neuron", 
         PCx = "PC1", PCy = "PC2",labels =F)

#for interneurons

