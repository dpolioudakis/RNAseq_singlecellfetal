#!/bin/bash

# Damon Polioudakis
# 2017-11-10
# Code to run Pooled_Vs_Bulk.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub Pooled_Vs_Bulk_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N PoolBulk
#$ -o logs/Pooled_Vs_Bulk_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Pooled_Vs_Bulk_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=256G,h_rt=12:00:00,highp
#$ -t 1-1
################################################################################
echo ""
echo "Starting Pooled_Vs_Bulk_QSUB.sh ${SGE_TASK_ID}... "$(date)
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

## Run Pooled_Vs_Bulk.R
${pathRscript} Pooled_Vs_Bulk.R
##########################################################################

echo ""
echo "End of Pooled_Vs_Bulk_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
