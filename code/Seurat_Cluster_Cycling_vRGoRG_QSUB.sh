#!/bin/bash

# Damon Polioudakis
# 2018-01-16
# Code to run Seurat_Cluster_Cycling_vRGoRG.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.4.0

# qsub:
# qsub Seurat_Cluster_Cycling_vRGoRG_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N vRGoRG
#$ -o logs/Seurat_Cluster_Cycling_vRGoRG_QSUB_$JOB_ID.log
#$ -e logs/Seurat_Cluster_Cycling_vRGoRG_QSUB_$JOB_ID.error
#$ -l h_data=128G,h_rt=24:00:00,highp
################################################################################
echo ""
echo "Starting Seurat_Cluster_Cycling_vRGoRG_QSUB.sh... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/bin/Rscript
# Set path to GCC compiler
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run Seurat_Cluster_Cycling_vRGoRG.R
${pathRscript} Seurat_Cluster_Cycling_vRGoRG_02.R
################################################################################

echo ""
echo "End of Seurat_Cluster_Cycling_vRGoRG_QSUB.sh... "$(date)
################################################################################
