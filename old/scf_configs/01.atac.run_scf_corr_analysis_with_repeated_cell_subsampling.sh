#!/bin/bash

# search
# kaList=(2 5 10 20 50 100 200)
# kaList+=(15 25 30 35 40 45 500)
# kaList+=(3 4 6 7 8 9 75 300 400)
kaList=(1000 2000 5000)
echo ${kaList[@]}


date="200826"
subsample_frac=0.8
subsample_times=3  # change to 3/10 in the future
nprocs=${subsample_times}

modx='smarter_cells'
mody='snatac_gene'

for ka in ${kaList[@]}; do
	knn=$ka
	# for k in ${knnList[@]}; do
		echo $ka, $knn

		# prep config file
		template="config_template_atac_scsmart_ka_knn.py"
		output="config_atac_scsmart_ka${ka}_knn${knn}_${date}.py"
		nameTagInConfigFile="mop_2mods_atacrna_${date}_ka${ka}_knn${knn}" # need to make sure they are consistent with the config_template

		cp $template $output 
		sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" $output 
		sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" $output 
		sed -i "s/date = TOBEFILLED/date = ${date}/g" $output 

		# run SCF (with 80% of cells, repeat n times)
		/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py \
			-c $output \
			-s ${subsample_frac} \
			-sn ${subsample_times}

		# # cross mod parallel
		# a list of samples
		sequences=""
		for (( i=0; i<${subsample_times}; i++ )) do
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

		# # cross mod
		# for (( i=0; i<${subsample_times}; i++ )) do
		# 	# # run corr analysis
		# 	echo $i
		# 	../correlation_analysis_celllevel-v4-stage3.py \
		# 		-tag $nameTagInConfigFile \
		# 		-modx $modx \
		# 		-mody $mody \
		# 		-ka $ka -knn $knn -dt $date -isub $i
		# 	# break
		# done

done	
