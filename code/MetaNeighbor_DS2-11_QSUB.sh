#!/bin/bash

# Damon Polioudakis
# 2017-08-31
# Code to run MetaNeighbor_DS2-11.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub MetaNeighbor_DS2-11_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N MetNeig
#$ -o logs/MetaNeighbor_DS2-11_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/MetaNeighbor_DS2-11_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=2:00:00,highp
#$ -t 2-2066
################################################################################
echo ""
echo "Starting MetaNeighbor_DS2-11_QSUB.sh ${SGE_TASK_ID}... "$(date)
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

## Run MetaNeighbor_DS2-11.R
${pathRscript} MetaNeighbor_DS2-11.R ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of MetaNeighbor_DS2-11_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
