#!/bin/bash

# Damon Polioudakis
# 2017-08-31
# Code to run Expression_GenesList.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub Expression_GenesList_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SC_Plot
#$ -o logs/Expression_GenesList_QSUB_$JOB_ID.log
#$ -e logs/Expression_GenesList_QSUB_$JOB_ID.error
#$ -l h_data=128G,h_rt=2:00:00,highp
################################################################################
echo ""
echo "Starting Expression_GenesList_QSUB.sh... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript

## Run Expression_GenesList.R
${pathRscript} Expression_GenesList.R
##########################################################################

echo ""
echo "End of Expression_GenesList_QSUB.sh... "$(date)
##########################################################################
