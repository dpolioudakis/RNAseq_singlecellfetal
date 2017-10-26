#!/bin/bash

# Damon Polioudakis
# 2017-03-30
# qsub Seurat_Cluster_Stability.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub Seurat_Cluster_Stability_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N StClStab
#$ -o logs/Seurat_Cluster_Stability_QSUB_$JOB_ID.log
#$ -e logs/Seurat_Cluster_Stability_QSUB_$JOB_ID.error
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 4
################################################################################

echo "Starting Seurat_Cluster_Stability_QSUB.sh ... "$(date)
################################################################################

# Path to R
RscriptPath=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript

# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

echo ""
echo "Running Seurat_Cluster_Stability.R..."
## Run Seurat_Cluster_Stability.R
${RscriptPath} Seurat_Cluster_Stability.R
################################################################################

echo ""
echo "End of Seurat_Cluster_Stability_QSUB.sh ... "$(date)
################################################################################
