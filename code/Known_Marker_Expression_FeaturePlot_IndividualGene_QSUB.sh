#!/bin/bash

# Damon Polioudakis
# 2018-06-21
# Code to run Known_Marker_Expression_FeaturePlot_IndividualGene.R
# Array length:
#   > length(km_DFL)
#   [1] 36

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.3+

# qsub:
# qsub Known_Marker_Expression_FeaturePlot_IndividualGene_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N KM_indv
#$ -o logs/Known_Marker_Expression_FeaturePlot_IndividualGene_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/Known_Marker_Expression_FeaturePlot_IndividualGene_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=64G,h_rt=12:00:00,highp
#$ -t 1-36
#$ -tc 6
################################################################################
echo ""
echo "Starting Known_Marker_Expression_FeaturePlot_IndividualGene_QSUB.sh ${SGE_TASK_ID}... "$(date)
echo ""
##########################################################################

### Run

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript

## Run Known_Marker_Expression_FeaturePlot_IndividualGene.R
${pathRscript} Known_Marker_Expression_FeaturePlot_IndividualGene.R ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of Known_Marker_Expression_FeaturePlot_IndividualGene_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
