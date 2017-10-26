#!/bin/bash

# Damon Polioudakis
# 2016-12-10
# Code to run Cluster_Seurat_KnownMarkers_Heatmap.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub Cluster_Seurat_KnownMarkers_Heatmap.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Seurat
#$ -o logs/Cluster_Seurat_KnownMarkers_Heatmap_$JOB_ID.log
#$ -e logs/Cluster_Seurat_KnownMarkers_Heatmap_$JOB_ID.error
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
################################################################################

### Set paths
# Path to R
pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

### Run Cluster_Seurat_KnownMarkers_Heatmap.R
${pathRscript} Cluster_Seurat_KnownMarkers_Heatmap.R
