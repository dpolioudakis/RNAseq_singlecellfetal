#!/bin/bash

# Damon Polioudakis
# 2017-05-28
# Code to run Seurat_Cluster_DS2-11.R as array job
# Task IDs should be number of times to run Seurat_Cluster_DS2-11.R
# (e.g. -t 1-4)

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0


# qsub:
# qsub Seurat_Cluster_DS2-11_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Seurat
#$ -o logs/Seurat_Cluster_DS2-11_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Seurat_Cluster_DS2-11_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=256G,h_rt=24:00:00,highp
#$ -t 1-1
################################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run Seurat_Cluster_DS2-11.R
${pathRscript} Seurat_Cluster_DS2-11.R

# /u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript Seurat_Cluster_DS2-11.R ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N701/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N702/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N703/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N704/out_gene_exon_tagged_dge_FtMm250.txt ../analysis/graphs/Seurat_Cluster_DS2-11_exon_FtMm250_
