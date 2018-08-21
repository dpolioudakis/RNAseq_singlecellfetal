#!/bin/bash

# Damon Polioudakis
# 2018-06-20
# Code to run Subplate.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub Subplate_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N subplate
#$ -o logs/Subplate_QSUB_$JOB_ID.log
#$ -e logs/Subplate_QSUB_$JOB_ID.error
#$ -l h_data=128G,h_rt=24:00:00,highp
################################################################################
echo ""
echo "Starting Subplate_QSUB.sh ... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript

## Run Subplate.R
${pathRscript} Subplate.R
##########################################################################

echo ""
echo "End of Subplate_QSUB.sh ... "$(date)
##########################################################################
