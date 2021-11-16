#!/bin/bash

# a grid search

# search
# kaList=(2 5 10 20 50 100 200)
# kaList+=(15 25 30 35 40 45 500)
# kaList+=(3 4 6 7 8 9 75 300 400)
kaList=(2 5 10 20 50 100 200 500 1000)
echo ${kaList[@]}

date="201108"
subsample_frac=0.8
subsample_times=3

# modalities
modx='10x_cells_v3'
mody='snatac_gene'

# scf config template
template="config_template_${modx}_${mody}_ka_knn.py"

nprocs=3

for ka in ${kaList[@]}; do
	knn=$ka
	# for k in ${knnList[@]}; do
	echo $ka, $knn
	# prep config file
	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
	scf_config="config_${nameTagInConfigFile}.py"

	cp $template ${scf_config} 
	sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
	sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
	sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

	# run SCF (with 80% of cells, repeat n times)
	/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py \
		-c ${scf_config} \
		-s ${subsample_frac} \
		-sn ${subsample_times}

	# get a list of samples
	sequences=""
	for (( i=0; i<${subsample_times}; i++ )); do
		# # run corr analysis
		sequences="${sequences}${i} "
	done
	# echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} echo "this is {}" 
	echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} \
		../correlation_analysis_celllevel_atac_nov9.py \
			-tag $nameTagInConfigFile \
			-modx $modx \
			-mody $mody \
			-ka $ka -knn $knn -dt $date -isub {}

done	
