#!/bin/bash

# Damon Polioudakis
# 2016-12-06
# Code to run Cluster_Seurat.R as array job
# Task ideas should be number of times to run Cluster_Seurat.R
# (e.g. -t 1-4)

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub Cluster_Seurat_QSUB.sh [File of Arguments for Cluster_Seurat.R]
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Seurat
#$ -o logs/Cluster_Seurat_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Cluster_Seurat_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
#$ -t 1-1
################################################################################

# Text file of arguments for Cluster_Seurat.R
inArgs=$1
# Path to R
pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

# Read lines until blank line and assign to elements of array
read -r -a argsA <<< $(awk -v taskID="${SGE_TASK_ID}" 'BEGIN {count=1}
  $0 == "" {count++}
  count == taskID {print $0}' ${inArgs})

# Stdout all Cluster_Seurat.R argsA on separate lines
echo "Arguments for Cluster_Seurat.R:"
printf '%s\n' "${argsA[@]}"

## Run Cluster_Seurat.R
${pathRscript} Cluster_Seurat.R ${argsA[0]} ${argsA[1]} ${argsA[2]} ${argsA[3]} ${argsA[4]} ${argsA[5]} ${argsA[6]} ${argsA[7]} ${argsA[8]} ${argsA[9]} ${argsA[10]} ${argsA[11]}

# /u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript Cluster_Seurat.R ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N701/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N702/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N703/out_gene_exon_tagged_dge_FtMm250.txt ../data/digital_gene_expression/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXbp083L1/N704/out_gene_exon_tagged_dge_FtMm250.txt ../analysis/graphs/Cluster_Seurat_exon_FtMm250_
