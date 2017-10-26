#!/bin/bash

# Damon Polioudakis
# 2017-05-16
# Code to run Seurat_InfoMap_Louvain.R as array job

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.4.0

# qsub:
# qsub Seurat_InfoMap_Louvain_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Seurat_InfoMap_Louvain
#$ -o logs/Seurat_InfoMap_Louvain_QSUB_$JOB_ID.log
#$ -e logs/Seurat_InfoMap_Louvain_QSUB_$JOB_ID.error
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/bin/Rscript
# Set path to GCC compiler
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run Seurat_InfoMap_Louvain.R
${pathRscript} Seurat_InfoMap_Louvain.R
################################################################################
