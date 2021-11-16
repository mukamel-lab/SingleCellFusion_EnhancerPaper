#!/bin/bash


kaList=(2 5 10 20 50 100 200)
kList=(2 5 10 20 50 100 200)

for ka in ${kaList[@]}; do
	for k in ${kList[@]}; do
		name="mc_scsmart_ka${ka}_k${k}"
		echo $ka, $k, $name

		# output="config_mc_scsmart_ka${ka}_k${k}.py"
		# cp config_template_mc_scsmart_ka_k.py $output 
		# sed -i "s/k = TOBEFILLED/k = ${k}/g" $output 
		# sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" $output 

		# # run SCF
		# /cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main.py -c $output > ../logs/${output}.log 2>&1

		# run corr analysis
		./correlation_analysis_celllevel-v2-stage3.py -ka $ka -k $k -dt "200803"  # > "./logs/${name}.log" 2>&1 


		# break
	done
	# break
done