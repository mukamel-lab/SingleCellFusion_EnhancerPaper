#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.insert(0, "/cndd/fangming/CEMBA/snmcseq_dev")
from multiprocessing import Pool,cpu_count
from functools import partial
from scipy import sparse
from scipy import stats
import pickle
import datetime
import argparse
import logging

import snmcseq_utils
from __init__ import *
from __init__jupyterlab import *

import CEMBA_clst_utils
import fbpca
import itertools
import argparse

def pipe_singlemod_clustering(
        f_gene, f_cell, f_mat, 
        flist_selected_cells, flist_res,
        resolutions,
        npc=50,
        kforleiden=30,
    ):
    """
    """
    # read input
    gc_mat = snmcseq_utils.load_gc_matrix(f_gene, f_cell, f_mat)

    for f_selected_cells, f_res in zip(flist_selected_cells, flist_res):
        print("processing {}".format(f_selected_cells))
        selected_cells = np.load(f_selected_cells, allow_pickle=True)
        # trim the matrix
        cg_mat_dense = gc_mat.data.T.todense()
        cg_mat_dense = cg_mat_dense[snmcseq_utils.get_index_from_array(gc_mat.cell, selected_cells)]

        # Louvain clustering for different resolutions
        # X should be selected from highly variable genes and normalized
        U, s, Vt = fbpca.pca(cg_mat_dense, k=npc)
        pcX = U.dot(np.diag(s))

        res_clsts = CEMBA_clst_utils.clustering_routine_multiple_resolutions(
                            pcX, selected_cells, kforleiden, 
                            seed=1, verbose=True,
                            resolutions=resolutions, metric='euclidean', option='plain', 
                            n_trees=10, search_k=-1, num_starts=None
                            )

        # organize and save results
        # res_clsts = res_clsts.reindex(gc_mat.cell)
        print(f_res)
        res_clsts.to_csv(f_res, sep='\t', na_rep='NA', header=True, index=True)
    return 

def wrapper_singlemod_clustering(
        mod, knns, subsample_times, 
        input_name_tag,
    ):
    """
    """
    # input data
    f_cell = '/cndd/fangming/CEMBA/data/MOp_all/data_freeze_neurons/{}_hvfeatures.cell'.format(mod)
    f_gene = '/cndd/fangming/CEMBA/data/MOp_all/data_freeze_neurons/{}_hvfeatures.gene'.format(mod)
    f_mat = '/cndd/fangming/CEMBA/data/MOp_all/data_freeze_neurons/{}_hvfeatures.npz'.format(mod)

    # input cell lists
    # input_name_tag = "mop_10x_cells_v3_snmcseq_gene_ka30_knn{{}}_201130"
    flist_selected_cells = [
        ('/cndd2/fangming/projects/miniatlas/results/cells_{}_{}.npy.{}.npy'
         .format(mod, input_name_tag.format(knn), i_sub))
            for knn, i_sub in itertools.product(knns, np.arange(subsample_times))
        ]
    # output files
    flist_res = [
        ('/cndd2/fangming/projects/miniatlas/results/clusterings_{}_{}_sub{}.tsv.gz'
         .format(mod, input_name_tag.format(knn), i_sub))
            for knn, i_sub in itertools.product(knns, np.arange(subsample_times))
        ]

    # resolution
    start, end = 0, 4
    resolutions = np.logspace(start, end, 10*(end-start)+1)
    print(resolutions)

    npc = 50
    kforleiden = 30

    pipe_singlemod_clustering(
        f_gene, f_cell, f_mat, 
        flist_selected_cells, flist_res,
        resolutions,
        npc=npc,
        kforleiden=kforleiden,
    )
    return 

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mod", help="modality", required=True)
    parser.add_argument("-ks", "--knns", help="a list of knns", nargs="+", required=True)
    parser.add_argument("-sn", "--subsample_times", help=">1", type=int, required=True)
    parser.add_argument("-tag", "--input_name_tag", help="input_name_tag", required=True)
    return parser

if __name__ == "__main__":
    # 
    parser = create_parser()
    args = parser.parse_args()

    # output setting
    # run this with each combination of (i_sub, knn)
    mod = args.mod
    knns = np.array(args.knns).astype(int)
    subsample_times = args.subsample_times
    input_name_tag = args.input_name_tag

    wrapper_singlemod_clustering(
        mod, knns, subsample_times, input_name_tag 
    )
