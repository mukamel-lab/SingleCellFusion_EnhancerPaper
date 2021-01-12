#!/bin/bash

# a grid search

# kaList=(15 25 30 35 40 45 500)
# kList=(2 5 10 20 50 100 200)

# search
kaList=(3 4 6 7 8 9 75 300 400)

for ka in ${kaList[@]}; do
	k=$ka
	# for k in ${kList[@]}; do
		echo $ka, $k
		output="config_mc_scsmart_ka${ka}_k${k}.py"
		cp config_template_mc_scsmart_ka_k.py $output 
		sed -i "s/k = TOBEFILLED/k = ${k}/g" $output 
		sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" $output 

		# for a given k run this
		# run SCF
		/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main.py -c $output > ../logs/${output}.log 2>&1
		# run corr analysis
		../correlation_analysis_celllevel-v2-stage3.py -ka $ka -k $k -dt "200803"  # > "./logs/${name}.log" 2>&1 
		# break
	# done
	# break
done