#!/bin/bash

# Damon Polioudakis
# 2017-08-31
# Code to run TF_Expression.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub TF_Expression_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N TFexpr
#$ -o logs/TF_Expression_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/TF_Expression_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=3:00:00,highp
#$ -t 1-1
################################################################################
echo ""
echo "Starting TF_Expression_QSUB.sh ${SGE_TASK_ID}... "$(date)
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

## Run TF_Expression.R
${pathRscript} TF_Expression.R
##########################################################################

echo ""
echo "End of TF_Expression_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
