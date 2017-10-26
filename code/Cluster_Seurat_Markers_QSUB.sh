#!/bin/bash

# Damon Polioudakis
# 2017-02-13
# Code to run Cluster_Seurat_Markers.R as array job

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub Cluster_Seurat_Markers_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Seurat
#$ -o logs/Cluster_Seurat_Markers_QSUB_$JOB_ID.log
#$ -e logs/Cluster_Seurat_Markers_QSUB_$JOB_ID.error
#$ -l h_data=8G,h_rt=24:00:00
#$ -pe shared 4
################################################################################

echo ""
echo "Starting Cluster_Seurat_Markers.sh ... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

echo "Running Cluster_Seurat_Markers.R..."
## Run Cluster_Seurat_Markers.R
${pathRscript} Cluster_Seurat_Markers.R
################################################################################

echo ""
echo "End of Cluster_Seurat_Markers.sh ... "$(date)
################################################################################
