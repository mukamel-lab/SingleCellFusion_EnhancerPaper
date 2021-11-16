#!/usr/bin/env python3
"""An example configuration file
"""
import sys
import os


# Assuming the cell order in the metadata tables are the same as those in the gene level matrices
# The output knn matrices follow such order as well 

ka_smooth = 15 
knn = 15 
date = 200826


# # Configs  
name = 'mop_2mods_atacrna_{}_ka{}_knn{}'.format(date, ka_smooth, knn)
outdir = '/cndd2/fangming/projects/miniatlas/results'
output_pcX_all = outdir + '/pcX_all_{}.npy'.format(name)
output_cells_all = outdir + '/cells_all_{}.npy'.format(name) 
output_imputed_data_format = outdir + '/imputed_data_{}_{{}}.npy'.format(name)
output_clst_and_umap = outdir + '/intg_summary_{}.tsv'.format(name)
output_figures = outdir + '/figures/{}_{{}}.{{}}'.format(name)
output_cluster_centroids = outdir + '/centroids_{}.pkl'.format(name)

save_knn = True # new required arguments (7/27/2020) 
output_knn_within = outdir + "/knn_within_{}_{{}}.npz".format(name)
output_knn_across = outdir + "/knn_across_{}_{{}}_{{}}.npz".format(name)
# end of new required arguments (7/27/2020)

# required for downsamp (8/7/2020)
output_cells = outdir + "/cells_{{}}_{}.npy".format(name)



DATA_DIR = '/cndd/fangming/CEMBA/data/MOp_all/data_freeze_neurons'
# fixed dataset configs
sys.path.insert(0, DATA_DIR)
from __init__datasets import *

meta_f = os.path.join(DATA_DIR, '{0}_metadata.tsv')
hvftrs_f = os.path.join(DATA_DIR, '{0}_hvfeatures.{1}')
hvftrs_gene = os.path.join(DATA_DIR, '{0}_hvfeatures.gene')
hvftrs_cell = os.path.join(DATA_DIR, '{0}_hvfeatures.cell')

mods_selected = [
    # 'snmcseq_gene',
    'snatac_gene',
    'smarter_cells',
    # 'smarter_nuclei',
    # '10x_cells_v2', 
    # '10x_cells_v3',
    # '10x_nuclei_v3',
    # '10x_nuclei_v3_macosko',
    ]

# features_selected = ['smarter_cells']
# features_selected = ['snmcseq_gene']
features_selected = ['snatac_gene']

# check features
for features_modality in features_selected:
    assert (features_modality in mods_selected)

# within modality
ps = {'mc': 0.9,
      'atac': 0.1,
      'rna': 0.7,
     }
drop_npcs = {
      'mc': 0,
      'atac': 0,
      'rna': 0,
     }
ka_smooth = ka_smooth # default: 5

# across modality
cross_mod_distance_measure = 'correlation' # cca
knn = knn 
relaxation = 3
n_cca = 30

# PCA
npc = 50

# clustering
k = 30 # default: 30
resolutions = [0.1, 1, 2, 4]
# umap
umap_neighbors = 60
min_dist = 0.5
