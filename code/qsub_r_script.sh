#!/bin/bash

# Damon Polioudakis
# 2018-11-15
# job submission script for R scripts

# to submit this script:
# qsub qsub_r_script.sh -p [path to R script]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N rscript
#$ -o logs/qsub_r_script_$JOB_ID.log
#$ -e logs/qsub_r_script_$JOB_ID.error
#$ -l h_data=64G,h_rt=12:00:00
# #$ -pe shared 8
#$ -t 1-16
# #$ -m bea

# reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting qsub_r_script.sh ${SGE_TASK_ID}... "$(date)
echo ""
################################################################################

### functions

usage () {
  message=$1
  if [ -n "${message}" ]
    then
      echo "${message}"
  fi
  echo "usage: $0"
  echo "  [-p path to R script]"
}
################################################################################

### main

## variables

# path to Rscript
rscript_path=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# arguments
# required
# path to script to run
in_r_script=

## handle arguments

while [ $# -gt 0 ]; do
  case "$1" in
  	-p) in_r_script="$2"; shift;;
  	--)	shift; break;;
    # usage message and terminates if an unknown command line flag starting with
    # a dash was specified
  	-*)
      echo >&2 usage "Invalid command line flag"
      exit 1;;
    *) break;;	# terminate while loop
  esac
  shift
done
# all command line switches are processed,
# "$@" contains all file names

# check arguments
# -z switch will test if the expansion of "$1" is a null string or not. If it is
# a null string then the body is executed.
if  [ -z "${in_r_script}" ]
  then
    usage "Missing required arguments"
    exit 1
fi

echo "Arguments:"
echo "Path to R script to run: ${in_r_script}"

## run Rscript
${rscript_path} ${in_r_script} ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of qsub_r_script.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
