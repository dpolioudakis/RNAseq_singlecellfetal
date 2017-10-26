#!/bin/bash

# Damon Polioudakis
# 2017-02-13
# Code to run Covariates_DS-002-003.R as array job

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Covariates_DS-002-003_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N Cov_Ds23
#$ -o logs/Covariates_DS-002-003_QSUB_$JOB_ID.log
#$ -e logs/Covariates_DS-002-003_QSUB_$JOB_ID.error
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
################################################################################

echo ""
echo "Starting Covariates_DS-002-003.sh ${SGE_TASK_ID}... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/current/bin/Rscript

# Run Covariates_DS-002-003.R
echo "Running Covariates_DS-002-003.R..."
${pathRscript} Covariates_DS-002-003.R
################################################################################

echo ""
echo "End of Covariates_DS-002-003.sh ${SGE_TASK_ID}... "$(date)
################################################################################
