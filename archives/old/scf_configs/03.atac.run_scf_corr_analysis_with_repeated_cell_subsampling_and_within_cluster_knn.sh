#!/bin/bash


date="200826" # 
# number of nearest neighbors (both within and across modalities)
kaList=(1000 2000 5000)
# kaList=(2 5 10 20 50 100 200 500)
# kaList=(50 100 200 500)
# kaList+=(15 25 30 35 40 45)
# kaList+=(3 4 6 7 8 9 75 300 400)
echo ${kaList[@]}

subsample_frac=0.8 # subsample cells
subsample_times=3 # number of repeats
nprocs=${subsample_times} 

# clusters
fname_cluster="/cndd/fangming/CEMBA/data/MOp_all/results_final/miniatlas_fig4_scf_clusterings.tsv"
colname_cluster="joint_cluster_round2"

# modalities
modx='smarter_cells'
mody='snatac_gene'

# scf config template
template="config_template_atac_scsmart_within_cluster_knn_ka_knn.py"

for ka in ${kaList[@]}; do
	knn=$ka
		echo $ka, $knn
		# prep config file
		scf_config="config_atac_scsmart_within_cluster_knn_ka${ka}_knn${knn}_${date}.py"
		nameTagInConfigFile="mop_2mods_atacrna_within_cluster_knn_${date}_ka${ka}_knn${knn}" # need to make sure they are consistent with the config_template
		cp $template ${scf_config} 
		sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
		sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
		sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

		# run SCF (with 80% of cells, repeat n times)
		/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling_within_cluster_knn.py \
			-c ${scf_config} \
			-s ${subsample_frac} \
			-sn ${subsample_times} \
			-f ${fname_cluster} \
			-col ${colname_cluster}

		# # cross mod correlation analysis parallel
		# a list of samples
		sequences=""
		for (( i=0; i<${subsample_times}; i++ )); do
			# # run corr analysis
			sequences="${sequences}${i} "
		done
		# echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} echo "this is {}" 
		echo "$sequences" | tr " " "\n"| xargs -P $nprocs -I {} \
			../correlation_analysis_celllevel-v4-stage3.py \
				-tag $nameTagInConfigFile \
				-modx $modx \
				-mody $mody \
				-ka $ka -knn $knn -dt $date -isub {}

		# # # cross mod
		# for (( i=0; i<${subsample_times}; i++ )); do
		# 	# # run corr analysis
		# 	echo $i
		# 	# same as 01.atac.
		# 	../correlation_analysis_celllevel-v4-stage3.py \
		# 		-tag $nameTagInConfigFile \
		# 		-modx $modx \
		# 		-mody $mody \
		# 		-ka $ka -knn $knn -dt $date -isub $i
		# done

done
