#!/bin/bash

# Will Connell
# 2018-06-26
# Code to run SC3_vs_seurat.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.4.0

# qsub:
# qsub SC3_vs_seurat_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SC3_vs_seurat
#$ -o logs/SC3_vs_seurat_QSUB_$JOB_ID.log
#$ -e logs/SC3_vs_seurat_QSUB_$JOB_ID.error
#$ -l h_data=256G,h_rt=24:00:00,highp
#$ -m bea
################################################################################
echo ""
echo "SC3_vs_seurat_QSUB.sh ... "$(date)
echo ""
################################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.5.0/gcc-6.3.0_MKL-2017/lib64/R/bin/Rscript

## Run single_cell_clustering.R
${pathRscript} SC3_vs_seurat.R
################################################################################

echo ""
echo "End of SC3_vs_seurat_QSUB.sh ... "$(date)
################################################################################
