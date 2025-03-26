#!/bin/bash

### Job name
#SBATCH -J viz_gI
#SBATCH -o logs/output_viz_gob_pan_three_inf.%j.%x.txt
#SBATCH -e logs/error_viz_gob_pan_three_inf.%j.%x.txt

### Time your job needs to execute, e. g. 15 min 30 sec 
#SBATCH -t 03:00:00
### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=50G

### OpenMP threads
#SBATCH --cpus-per-task=8

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jschoeneich@ukaachen.de

################################################################
#module load R
module load scRNA
#source ~/miniconda3/bin/activate
#conda activate Seurat3
################################################################

# could also define it here instead of taking it as an arg
proj_name="scrna_full_2024_both_d25_gob_pan_three_inf"
# data dir (where your results were saved)
# the way it is set up below -o for output will also save your report
# on that same directory path
data_path="/data/scRNA/johannes_scRNA/kinetics/scrna_seurat_pipeline/scrna_full_2024/scrna_full_2024_both_d25_gob_pan_three"
#data_path=`pwd`
#cluster="removed_clusters"
#cluster="remove_recluster"
#cluster="merged_clusters"
#cluster="annotation"
#cluster="singleton"

cluster="gob_pan_mix"

# QCC is QC for existing RObject(no rawdata)
FUNCS=(
 #Singleton
 QC
 QCC
 AmbientRNA
 DoubletDetection
 DEs
 Clusters
 Clusters_harmony
 Clusters_seurat
 EXT_MARKERS
 DEGO
 #Genesets
 progeny
 hallmark
 KEGG
 Reactome
 Genesets_stage
 progeny_stage
 DEGO_stage
 hallmark_stage
 reactome_stage
 kegg_stage
 Genesets_1v1
 DEGO_1v1
 hallmark_1v1
 reactome_1v1
 kegg_1v1
 intUMAPs
 #GET_DATA
)

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }
join_str=`join_by '", "' ${FUNCS[@]}`
json_exe_list="[\"${join_str}\"]"

date
Rscript ./viz/create_report.R \
	-a "Johannes" \
	-p ${proj_name} \
	-m "TRUE" \
	-s "${data_path}/${proj_name}/save" \
	-c "./conf/config_${proj_name}.R" \
	-o "${data_path}/${proj_name}/report_gob_pan" \
	-r "${data_path}/${proj_name}/charts_gob_pan" \
	-e ./external/Human_and_mouse_cell_markers-Markers.tsv \
	-v "./static/viz_module.ini" \
  -d "$cluster" \
  -j "${json_exe_list}" \
  -i "FALSE" \
  -f "seurat" \
  -z "lz4"
date


######Options:############################################
#	-h, --help
#		Show this help message and exit
#
#	-a CHARACTER, --author=CHARACTER
#		Author display in report [default Costalab]
#
#	-m CHARACTER, --make_element=CHARACTER
#		make_element [default FALSE]
#
#	-s CHARACTER, --savedir=CHARACTER
#		savedir [default save]
#
#	-c CHARACTER, --configfile=CHARACTER
#		configfile [default conf/config.R]
#
#	-o CHARACTER, --report_dir=CHARACTER
#		report_dir [default report]
#
#	-r CHARACTER, --charts_dir=CHARACTER
#		charts_dir [default charts]
#
#	-p CHARACTER, --project=CHARACTER
#		charts_dir [default ]
#
#	-e CHARACTER, --externalfile=CHARACTER
#		report_dir [default ./external/Human_and_mouse_cell_markers-Markers.tsv]
#
#	-d CHARACTER, --defaultclsuters=CHARACTER
#		cluster slots:
#               removed_clusters
#               remove_reclusters
#               merged_clusters
#               annotation
#               singleton [default seurat_clusters]
#
#	-j CHARACTER, --planOfreport=CHARACTER
#		plan of report [default c()]
#
#	-l CHARACTER, --singlefile=CHARACTER
#		generate single html file [default FALSE]
#
#	-i CHARACTER, --indexonly=CHARACTER
#		only generate index.html [default FALSE]

