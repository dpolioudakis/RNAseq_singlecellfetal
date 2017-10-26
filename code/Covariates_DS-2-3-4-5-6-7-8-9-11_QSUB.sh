#!/bin/bash

# Damon Polioudakis
# 2017-04-17
# Code to run Covariates_DS-2-3-4-5-6-7-8-9-11.R as array job

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Covariates_DS-2-3-4-5-6-7-8-9-11_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Cov_Ds
#$ -o logs/Covariates_DS-2-3-4-5-6-7-8-9-11_QSUB_$JOB_ID.log
#$ -e logs/Covariates_DS-2-3-4-5-6-7-8-9-11_QSUB_$JOB_ID.error
#$ -l h_data=128G,h_rt=24:00:00,highp
################################################################################

echo ""
echo "Starting Covariates_DS-2-3-4-5-6-7-8-9-11.sh... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

# Run Covariates_DS-2-3-4-5-6-7-8-9-11.R
echo "Running Covariates_DS-2-3-4-5-6-7-8-9-11.R..."
${pathRscript} Covariates_DS-2-3-4-5-6-7-8-9-11.R
################################################################################

echo ""
echo "End of Covariates_DS-2-3-4-5-6-7-8-9-11.sh... "$(date)
################################################################################
