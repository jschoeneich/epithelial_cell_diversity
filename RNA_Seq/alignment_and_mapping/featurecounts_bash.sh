#!/bin/bash

### Job name
#SBATCH -J fcounts

### Set logs, remember to create logs dir
#SBATCH -e logs/log_fc_full_%j.log
#SBATCH -o logs/log_fc_full_%j.log

### Time to execute, e. g. 15 min 30 sec
#SBATCH -t 6:00:00

### Job memory needs per node, e. g. 1 GB
#SBATCH --mem=90G

### OpenMP threads
#SBATCH --cpus-per-task=10

################################################################
# PATH
if [ -r /usr/local_host/etc/bashrc ]; then
   . /usr/local_host/etc/bashrc
fi

export PATH=/usr/bin:$PATH
export PATH=/usr/local_host/bin:$PATH
################################################################
module load R/4.2.3
#source ~/miniconda3/bin/activate
#conda activate Seurat3
################################################################


date

Rscript featurecounts_hpc.R

date
