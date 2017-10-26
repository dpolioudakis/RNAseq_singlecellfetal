#!/bin/bash

# Damon Polioudakis
# 2017-10-04
# Code to run Seurat_ClusterDE_TFenrichment.R as array job
# Task IDs should be number of clusters
# (e.g. -t 1-4)

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+


# qsub:
# qsub Seurat_ClusterDE_TFenrichment_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N SrtDETF
#$ -o logs/Seurat_ClusterDE_TFenrichment_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Seurat_ClusterDE_TFenrichment_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=3:00:00,highp
#$ -t 1-18
##########################################################################
echo ""
echo "Starting Seurat_ClusterDE_TFenrichment_QSUB.sh ${SGE_TASK_ID}... "$(date)
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

## Run Seurat_ClusterDE_TFenrichment.R
${pathRscript} Seurat_ClusterDE_TFenrichment_JASPAR.R ${SGE_TASK_ID};
${pathRscript} Seurat_ClusterDE_TFenrichment_TRANSFAC.R ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of Seurat_ClusterDE_TFenrichment_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
