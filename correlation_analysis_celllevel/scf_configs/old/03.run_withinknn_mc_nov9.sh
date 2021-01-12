#!/bin/bash
# Cell level correlation analysis
# mC - RNA
# - run SCF integration with downsampled cells x times (kNN by within cluster kNNs only?)
# - build pseudo-cells according to within kNN connections only

date="201108" # 
# number of nearest neighbors (both within and across modalities)
# kaList=(1000 2000 5000)
# kaList=(2 5 10 20 50 100 200 500)
# kaList=(50 100 200 500)
# kaList+=(15 25 30 35 40 45)
# kaList+=(3 4 6 7 8 9 75 300 400)
kaList=(2 5 10 20 50 100 200 500 1000)

subsample_frac=0.8 # subsample cells
subsample_times=3 # number of repeats
nprocs=${subsample_times} 

# clusters
fname_cluster="/cndd/fangming/CEMBA/data/MOp_all/results_final/miniatlas_fig4_scf_clusterings.tsv"
colname_cluster="joint_cluster_round2"

# modalities
modx='10x_cells_v3'
mody='snmcseq_gene'

# scf config template
template="config_template_${modx}_${mody}_within_cluster_knn_ka_knn.py"
echo ${kaList[@]}

for ka in ${kaList[@]}; do
	knn=$ka
	echo $ka, $knn
	# prep config file
	nameTagInConfigFile="mop_${modx}_${mody}_within_cluster_knn_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
	scf_config="config_${nameTagInConfigFile}.py"

	cp $template ${scf_config} 
	sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
	sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
	sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

	# run SCF (with 80% of cells, repeat n times, kNN within clusters)
	/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling_within_cluster_knn.py \
		-c ${scf_config} \
		-s ${subsample_frac} \
		-sn ${subsample_times} \
		-f ${fname_cluster} \
		-col ${colname_cluster}

	# # cross mod correlation analysis parallel
	# get a list of samples
	sequences=""
	for (( i=0; i<${subsample_times}; i++ )); do
		# # run corr analysis
		sequences="${sequences}${i} "
	done
	# echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} echo "this is {}" 
	echo "$sequences" | tr " " "\n"| xargs -P $nprocs -I {} \
		../correlation_analysis_celllevel_mc_nov9.py \
			-tag $nameTagInConfigFile \
			-modx $modx \
			-mody $mody \
			-ka $ka -knn $knn -dt $date -isub {}

done
