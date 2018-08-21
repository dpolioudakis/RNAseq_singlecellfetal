#!/bin/bash

# Damon Polioudakis
# 2018-05-31
# Code to run Seurat_ClusterRound2_MarkerPlotsIndv.R

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub Seurat_ClusterRound2_MarkerPlotsIndv_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N StR2_KMi
#$ -o logs/Seurat_ClusterRound2_MarkerPlotsIndv_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Seurat_ClusterRound2_MarkerPlotsIndv_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=64G,h_rt=24:00:00
#$ -t 1-16
#$ -tc 8
################################################################################
echo ""
echo "Starting Seurat_ClusterRound2_MarkerPlotsIndv_QSUB.sh ${SGE_TASK_ID}... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript

## Run Seurat_ClusterRound2_MarkerPlotsIndv.R
${pathRscript} Seurat_ClusterRound2_MarkerPlotsIndv.R ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of Seurat_ClusterRound2_MarkerPlotsIndv_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
