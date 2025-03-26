#!/bin/bash

### Job name
#SBATCH -J sc2d25v2

### Set logs, remember to create logs dir
#SBATCH -e logs/scrna_full_2024_both_d25_error.%j.%x.txt
#SBATCH -o logs/scrna_full_2024_both_d25_output.%j.%x.txt

### Time to execute, e. g. 15 min 30 sec
#SBATCH -t 16:00:00

### Job memory needs per node, e. g. 1 GB
#SBATCH --mem=250G

### OpenMP threads
#SBATCH --cpus-per-task=16

################################################################
# PATH
if [ -r /usr/local_host/etc/bashrc ]; then
   . /usr/local_host/etc/bashrc
fi

export PATH=/usr/bin:$PATH
export PATH=/usr/local_host/bin:$PATH
################################################################
module load scRNA
#source ~/miniconda3/bin/activate
#conda activate Seurat3
################################################################

proj_name="scrna_full_2024_both_d25_v2"
# data dir (where your results will be saved)

data_path='scrna_full_2024'

date
## 50 cores run, future memory
mkdir -p ${data_path}/${proj_name}
ln -s ${data_path}/conf/config_${proj_name}.R ${data_path}/${proj_name}
Rscript data_factory.R \
  -n 16 \
  --MaxMemMega=250000 \
  -z "lz4" \
  -c "./conf/config_${proj_name}.R" \
  -s "${data_path}/${proj_name}/save" \
  -e "${data_path}/${proj_name}/charts" \
  -a seurat_clusters \
  --allinone=TRUE \
  --nFeatureRNAfloor=600 \
  --nCountRNAfloor=1000 \
  --nCountRNAceiling=40000 \
  --pct_mitoceiling=25 \
  --pct_riboceiling=100 \
  --countmatrixformat='10X' \
  --dims4Integrate="1:30" \
  --Dims4FindNeighbors="1:50" \
  --harmony_dim="1:50" \
  -r 0.3

date
