#!/bin/bash

# Damon Polioudakis
# 2017-12-10
# Code to run Monocle_Round2_Analysis_Telley.R

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Monocle_Round2_Analysis_Telley_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N MnTelley
#$ -o logs/Monocle_Round2_Analysis_Telley_QSUB_$JOB_ID.log
#$ -e logs/Monocle_Round2_Analysis_Telley_QSUB_$JOB_ID.error
#$ -l h_data=128G,h_rt=8:00:00,highp
################################################################################

echo ""
echo "Starting Monocle_Round2_Analysis_Telley_QSUB.sh ... "$(date)
echo ""
################################################################################

# Path to R
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

# Run Monocle_Round2_Analysis_Telley.R
echo "Running Monocle_Round2_Analysis_Telley.R..."
${pathRscript} Monocle_Round2_Analysis_Telley.R
################################################################################

echo ""
echo "End of Monocle_Round2_Analysis_Telley_QSUB.sh ... "$(date)
################################################################################
