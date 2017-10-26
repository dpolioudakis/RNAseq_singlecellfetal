#!/bin/bash

# Damon Polioudakis
# 2017-07-17
# Compile statistics for Illumina PF filter
################################################################################

echo ""
echo "Starting 1_Process_And_Align.sh..."$(date)
echo ""


inQSEQlog=../DS-009/code/logs/0_QSEQ_to_fastq_QSUB_51109_demultiplex.log

grep "Working on directory:.*" ${inQSEQlog}
grep "Number of unfiltered reads in the lane:.*" ${inQSEQlog}
grep "Number of filtered reads in the lane:.*" ${inQSEQlog}
grep "Percentage of useable reads:.*" ${inQSEQlog}

################################################################################

echo ""
echo "End of 1_Process_And_Align.sh... "$(date)
################################################################################
