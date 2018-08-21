#!/bin/bash

# Damon Polioudakis
# 2017-08-15
# Code to run Monocle.R as array job

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Monocle_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Monocle
#$ -o logs/Monocle_QSUB_$JOB_ID.log
#$ -e logs/Monocle_QSUB_$JOB_ID.error
#$ -l h_data=256G,h_rt=24:00:00,highp
################################################################################

echo ""
echo "Starting Monocle.sh... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

# Run Monocle.R
echo "Running Monocle.R..."
${pathRscript} Monocle.R
################################################################################

echo ""
echo "End of Monocle.sh... "$(date)
################################################################################
