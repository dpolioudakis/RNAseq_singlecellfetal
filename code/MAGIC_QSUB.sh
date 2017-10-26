#!/bin/bash

# Damon Polioudakis
# 2017-08-11
# Code to run MAGIC.py imputation

# Reminder:
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load python/3.6.1

# qsub:
# qsub MAGIC_QSUB.sh [file with magic.py commands]
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N MAGIC_QSUB
#$ -o logs/MAGIC_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/MAGIC_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=12:00:00,highp
#$ -t 1-4
################################################################################

# File of magic.py commands
inCmds=$1

# Loop through file of magic.py commands and store in array
i=0
while read line
do
  # Skip blank lines
  # The -z test means true if empty. ! negates (i.e. true if not empty).
  if [ ! -z "$line" ]; then
        arr[$i]="$line"
        i=$((i+1))
      fi
done < ${inCmds}
echo "Commands to execute:"
printf '%s\n' "${arr[@]}"

# Run each command in array
echo "Task ID:"
echo ${SGE_TASK_ID}
idx=$((SGE_TASK_ID-1))
echo "Array position:"
echo ${idx}
eval "${arr[${idx}]}"
################################################################################




# #!/bin/python
#
# # Remember to specify python 3+ with "python3"
#
# # It looks like you're running in interactive mode by default, so matplotlib wants to plot everything to the screen first, which of course it can't do.
#
#
#
#
#
# # Plotting and miscellaneous imports
# import os
# import matplotlib
# # If you only want to plot figures to files (pdf, png,..), set matplotlib screen output to  'Agg' before importing pylab or pyplot.
# matplotlib.use('Cairo')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
#
# import magic
#
# %matplotlib inline
#
# # Try putting
# matplotlib.pyplot.ioff()
#
#
# # Load single-cell RNA-seq data
# # help(magic.mg.SCData.from_csv)
# scdata = magic.mg.SCData.from_csv(os.path.expanduser('/u/home/d/dpolioud/project-geschwind/RNAseq_singlecellfetal/analysis/tmp.csv'), data_type='sc-seq', normalize=False, cell_axis=1)
#
# scdata
#
# fig, ax = scdata.plot_molecules_per_cell_and_gene()
# savefig('MAGIC_1.png')
# plt.close(fig)
#
# # Minimum molecules/cell value
# CELL_MIN = 0
#
# # Maximum molecules/cell values
# CELL_MAX = 1000000
#
# # Minimum number of nonzero cells/gene
# # (None if no filtering desired)
# GENE_NONZERO = None
#
# # Minimum number of molecules/gene
# # (None if no filtering desired)
# GENE_MOLECULES = None
#
# scdata.filter_scseq_data(filter_cell_min=CELL_MIN, filter_cell_max=CELL_MAX,
#                          filter_gene_nonzero=GENE_NONZERO, filter_gene_mols=GENE_MOLECULES)
#
# scdata = scdata.normalize_scseq_data()
#
# fig, ax = scdata.plot_pca_variance_explained(n_components=40, random=True)
#
#
#
#
# # Make a square figure and axes
# figure(1, figsize=(6, 6))
# ax = axes([0.1, 0.1, 0.8, 0.8])
#
# labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
# fracs = [15, 30, 45, 10]
#
# explode = (0, 0.05, 0, 0)
# pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
# title('Raining Hogs and Dogs', bbox={'facecolor': '0.8', 'pad': 5})
#
#
