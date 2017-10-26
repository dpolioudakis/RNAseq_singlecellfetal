#!/bin/bash

# Damon Polioudakis
# 2017-02-21
# qsub Seurat_Markers.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub Seurat_Markers_QSUB.sh [Seurat object]
#  [Seurat_Markers.R output path graphs]
#  [Seurat_Markers.R output path tables]
#  [Seurat_Markers.R output path Rdata]
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SrtMarks
#$ -o logs/Seurat_Markers_QSUB_$JOB_ID.log
#$ -e logs/Seurat_Markers_QSUB_$JOB_ID.error
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
################################################################################

echo "Starting Seurat_Markers_QSUB.sh ... "$(date)
################################################################################

# Seurat object
inSeuratO=$1
# Seurat_Markers.R output paths
# Graphs
outGraph=$2
# Tables
outTable=$3
# Rdata
outRdata=$4

# Path to R
RscriptPath=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript

# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

echo ""
echo "Running Seurat_Markers.R..."
## Run Seurat_Markers.R
${RscriptPath} Seurat_Markers.R ${inSeuratO} ${outGraph} ${outTable} ${outRdata}
################################################################################

echo ""
echo "End of Seurat_Markers_QSUB.sh ... "$(date)
################################################################################
