#!/bin/bash

# a grid search

# knnList=(2 5 10 20 50 100 200)

# search
kaList=(2 5 10 20 50 100 200)
kaList+=(15 25 30 35 40 45 500)
kaList+=(3 4 6 7 8 9 75 300 400)
echo ${kaList[@]}


date="200811" # 
subsample_frac=0.8
subsample_times=10
fname_cluster="/cndd/fangming/CEMBA/data/MOp_all/results_final/miniatlas_fig4_scf_clusterings.tsv"
colname_cluster="joint_cluster_round2"

for ka in ${kaList[@]}; do
	knn=$ka
		echo $ka, $knn

		# # prep config file
		# output="config_mc_scsmart_within_cluster_knn_ka${ka}_knn${knn}_${date}.py"
		# cp "config_template_mc_scsmart_within_cluster_knn_ka_knn.py" $output 
		# sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" $output 
		# sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" $output 
		# sed -i "s/date = TOBEFILLED/date = ${date}/g" $output 

		# # run SCF (with 80% of cells, repeat n times)
		# /cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling_within_cluster_knn.py \
		# 	-c $output \
		# 	-s ${subsample_frac} \
		# 	-sn ${subsample_times} \
		# 	-f ${fname_cluster} \
		# 	-col ${colname_cluster}

		for (( i=0; i<${subsample_times}; i++ )); do
			# # run corr analysis
			echo $i
			../correlation_analysis_celllevel-v3-stage4.py -ka $ka -knn $knn -dt $date -isub $i
		done
done