#!/bin/bash

date="201108"

# kaList=(2 5 10 20 50 100 200)
# kaList+=(15 25 30 35 40 45 500)
# kaList+=(3 4 6 7 8 9 75 300 400)
kaList=(2 5 10 20 50 100 200 500 1000)

subsample_frac=0.8
subsample_times=1
nprocs=${subsample_times}

# modalities
modx='smarter_cells'
mody='snmcseq_gene'

# scf config template
template="./configs/config_template_${modx}_${mody}_ka_knn.py"

echo ${kaList[@]}

for ka in ${kaList[@]}; do
	knn=$ka
	echo $ka, $knn
	# prep config file
	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
	scf_config="./configs/config_${nameTagInConfigFile}.py"

	# fill knn, ka_smooth, date
	cp $template ${scf_config} 
	sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
	sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
	sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

	# run SCF (with 80% of cells, repeat n times)
	/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py \
		-c ${scf_config} \
		-s ${subsample_frac} \
		-sn ${subsample_times}

	# run correlation analysis
	# get a list of samples
	sequences=""
	for (( i=0; i<${subsample_times}; i++ )); do
		# # run corr analysis
		sequences="${sequences}${i} "
	done
	# echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} echo "this is {}" 
	echo "$sequences" | tr " " "\n" | xargs -P $nprocs -I {} \
		../correlation_analysis_celllevel_mc.py \
			-tag $nameTagInConfigFile \
			-modx $modx \
			-mody $mody \
			-ka $ka -knn $knn -dt $date -isub {}

	break # for testing only
done	
