#!/bin/bash

# Damon Polioudakis
# 2018-01-02
# Code to run single_cell_clustering.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3.0

# qsub:
# qsub single_cell_clustering_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SC3
#$ -o logs/single_cell_clustering_QSUB_$JOB_ID.log
#$ -e logs/single_cell_clustering_QSUB_$JOB_ID.error
#$ -l h_data=256G,h_rt=24:00:00,highp
################################################################################
echo ""
echo "Starting single_cell_clustering_QSUB.sh ... "$(date)
echo ""
################################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run single_cell_clustering.R
${pathRscript} single_cell_clustering.R ExprMatrix.txt out --input_tpm INPUT_TPM --no_log NO_LOG
################################################################################

echo ""
echo "End of single_cell_clustering_QSUB.sh ... "$(date)
################################################################################
