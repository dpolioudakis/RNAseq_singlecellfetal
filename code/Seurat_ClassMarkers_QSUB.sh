#!/bin/bash

# Damon Polioudakis
# 2018-05-21
# Code to run Seurat_ClassMarkers.R as array job
# Task IDs should be number of cell classes
# length(unique(names(class_cluster_idx)))
# (e.g. -t 1-9)

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+


# qsub:
# qsub Seurat_ClassMarkers_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SrtClasM
#$ -o logs/Seurat_ClassMarkers_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Seurat_ClassMarkers_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=12:00:00,highp
#$ -t 1-9
##########################################################################
echo ""
echo "Starting Seurat_ClassMarkers.sh ${SGE_TASK_ID}... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run Seurat_ClassMarkers.R
${pathRscript} Seurat_ClassMarkers.R ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of Seurat_ClassMarkers.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
