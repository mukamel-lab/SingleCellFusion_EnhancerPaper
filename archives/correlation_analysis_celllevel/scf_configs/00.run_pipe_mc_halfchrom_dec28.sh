#!/bin/bash
# thoughts 
# 1. run_scf_repeatedly_ with different knn and downsampling
# output: knn_xy matrices, and cells (i, knn)

# 2. within-modal clustering ...
# for each cell set from (i, knn) -> run Leiden clustering with different resolution r
# output: clustering result.tsv (i, knn, r): (i, knn) separate tables, r columns per table

# 3. correlation analysis ...
# for each combination of (i, knn, r), run a correlation analysis
# output: corrs.pkl (i, knn, r) 

date="210101" # 201130, 201206, 201228 

num_metacell_limit=10001
# generate (i, knn) knn_xy matrices
# modalities
modx='10x_cells_v3'
mody='snmcseq_gene'
ka=30
knns=(30 5 10 50 100)
corr_type='spearmanr'
# bootstrapping 80% of cells 5 times repeatedly
subsample_frac=1
subsample_times=1
feature_subsample_frac=0.5
feature_subsample_times=3
gene_annot="/cndd2/fangming/references/gencode.vM16.annotation_genes_19x.bed"

# for knn in ${knns[@]}; do
# 	echo $modx, $mody, $ka, $knn

# 	# scf config template
# 	template="./configs/config_template_${modx}_${mody}_ka_knn.py"

# 	# prep config file
# 	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # need to make sure they are consistent with the config_template
# 	scf_config="./configs/config_${nameTagInConfigFile}.py"

# 	# fill knn, ka_smooth, date
# 	cp $template ${scf_config} 
# 	sed -i "s/knn = TOBEFILLED/knn = ${knn}/g" ${scf_config} 
# 	sed -i "s/ka_smooth = TOBEFILLED/ka_smooth = ${ka}/g" ${scf_config} 
# 	sed -i "s/date = TOBEFILLED/date = ${date}/g" ${scf_config} 

# 	# # run SCF
# 	/cndd2/fangming/projects/SCF_package/SingleCellFusion/scripts/SCF_main_repeat_subsampling_subset_features.py \
# 		-c ${scf_config} \
# 		-s ${subsample_frac} \
# 		-sn ${subsample_times} \
# 		-fs ${feature_subsample_frac} \
# 		-fsn ${feature_subsample_times} \
# 		-ga ${gene_annot}
# done

# # run leiden clustering for each (i, knn) 
# # get a list of samples
# inputNameTag="mop_${modx}_${mody}_ka${ka}_knn{}_${date}"
# ../generate_rna_dataset_clustering.py --mod $modx --knns ${knns[@]} -sn ${subsample_times} -tag $inputNameTag

# correlation analysis (i, j, knn, r)
# most time consuming and memory consuming
# limit r to 10^3

for knn in ${knns[@]}; do
	nameTagInConfigFile="mop_${modx}_${mody}_ka${ka}_knn${knn}_${date}" # 	
	for (( i=0; i<${subsample_times}; i++ )); do
		for (( j=0; j<${feature_subsample_times}; j++ )); do
		# 	# # run corr analysis (for different r (cluster resolution))
		../correlation_analysis_celllevel_mcrna_halfchrom_Dec28.py \
			-modx $modx \
			-mody $mody \
			-tag $nameTagInConfigFile \
			-isub $i \
			-ifsub $j \
			--corr_type ${corr_type} \
			-n ${num_metacell_limit}
		# break
		done
	# break
	done
# break
done

