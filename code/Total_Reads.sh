#!/bin/bash

# Damon Polioudakis
# 2017-03-28
# Calculate total lines in fastq
################################################################################

# Reminder:
#  make /logs directory in code directory

# qsub:
# qsub Total_Reads.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N TotReads
#$ -o logs/Total_Reads_$JOB_ID.log
#$ -e logs/Total_Reads_$JOB_ID.error
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 4
################################################################################

echo "Starting Total_Reads.sh ... "$(date)
################################################################################

outFile="../metadata/Total_Reads.txt"

echo "Lines in fastq R1:" > ${outFile}

echo "DS-001 SxaQSEQsXap108L1" >> ${outFile}
gzip -cd ../DS-001_DP/data/fastq/SxaQSEQsXap108L1/*R1.fastq.gz | wc -l >> ${outFile}

echo "DS-001 SxaQSEQsXap110L1" >> ${outFile}
gzip -cd ../DS-001_DP/data/fastq/SxaQSEQsXap110L1/*R1.fastq.gz | wc -l >> ${outFile}

echo DS-002 >> ${outFile}
gzip -cd ../DS-002_DP/data/fastq/SxaQSEQsXbp083L1/*R1.fastq.gz | wc -l >> ${outFile}

echo DS-003 >> ${outFile}
gzip -cd ../DS-003_DP/data/fastq/SxaQSEQsVAP048L8/*R1.fastq.gz | wc -l >> ${outFile}

echo DS-003v2 >> ${outFile}
gzip -cd ../DS-003v2_DP/data/fastq/Sxa*/*R1.fastq.gz | wc -l >> ${outFile}
# Adding up picard output alignment_stats.txt came out to 357,567,664
# Lines in fastq divided by 4 = 357,674,855
################################################################################

echo "End of Total_Reads.sh ... "$(date)
################################################################################
