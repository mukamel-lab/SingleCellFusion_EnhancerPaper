#!/bin/bash

# 1. run_scf_repeatedly_ with different knn and downsampling
# output: knn_xy matrices, and cells (i, knn)

# 2. within-modal clustering ...
# for each cell set from (i, knn) -> run Leiden clustering with different resolution r
# output: clustering result.tsv (i, knn, r): (i, knn) separate tables, r columns per table

# 3. correlation analysis ...
# for each combination of (i, knn, r), run a correlation analysis
# output: corrs.pkl (i, knn, r) 

date="211115" # still use the (finished) 1130 results

num_metacell_limit=1001
# generate (i, knn) knn_xy matrices
# modalities
modx='rna'
mody='atac'
ka=30
knns=(30)
corr_type='spearmanr'
subsample_frac=0.1
subsample_times=1
resolutions=(1) # Leiden clustering resolutions used to generate metacells -- just 1 for demo

# 1.
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
	./scf/SCF_main_repeat_subsampling.py \
		-c ${scf_config} \
		-s ${subsample_frac} \
		-sn ${subsample_times}
done

# 2.
# run leiden clustering for each (i, knn) 
# get a list of samples
inputNameTag="mop_${modx}_${mody}_ka${ka}_knn{}_${date}"
./generate_metacells_rna.py --mod $modx --knns ${knns[@]} -sn ${subsample_times} -tag $inputNameTag -r ${resolutions}

# 3.
# correlation analysis (i, knn, r)
for knn in ${knns[@]}; do
	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # 	
	for (( i=0; i<${subsample_times}; i++ )); do
	# 	# # run corr analysis
	./correlate_metacells_atac_rna.py \
		-modx $modx \
		-mody $mody \
		-tag $nameTagInConfigFile \
		-isub $i \
		--corr_type ${corr_type} \
		-n ${num_metacell_limit} \
		-f
	done
done