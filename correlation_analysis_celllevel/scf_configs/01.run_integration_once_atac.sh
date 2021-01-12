#!/bin/bash

date="201120"
# modalities
modx='10x_cells_v3'
mody='snatac_gene'
ka=30 #smoothing
# knn=30
subsample_frac=1
subsample_times=1
knns=(50 100)
for knn in ${knns[@]}; do
	echo $modx, $mody, $ka, $knn
	# scf config template
	template="./configs/config_template_${modx}_${mody}_ka_knn.py"

	# prep config file
	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
	scf_config="./configs/config_${nameTagInConfigFile}.py"

	# fill knn, ka_smooth, date
	cp $template ${scf_config} 
	sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
	sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
	sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

	# # run SCF
	/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling.py \
		-c ${scf_config} \
		-s ${subsample_frac} \
		-sn ${subsample_times}

done



