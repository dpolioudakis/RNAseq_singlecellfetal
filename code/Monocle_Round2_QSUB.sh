#!/bin/bash

# Damon Polioudakis
# 2017-08-15
# Code to run Monocle_Round2.R as array job

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Monocle_Round2_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Monocle_Round2
#$ -o logs/Monocle_Round2_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Monocle_Round2_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=24:00:00,highp
#$ -t 1-23
################################################################################

echo ""
echo "Starting Monocle_Round2_QSUB.sh ${SGE_TASK_ID}... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

# Run Monocle_Round2.R
echo "Running Monocle_Round2.R..."
${pathRscript} Monocle_Round2.R ${SGE_TASK_ID}
################################################################################

echo ""
echo "End of Monocle_Round2_QSUB.sh ${SGE_TASK_ID}... "$(date)
################################################################################
