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

sys.path.insert(0, "/cndd2/fangming/projects/scf_enhancers/scripts/scf_enhancer_paper")
import enhancer_gene_utils

logger = snmcseq_utils.create_logger()

def pipe_corr_analysis_mc(
        common_rna_cells, common_mc_cells,
        cell_cell_knn_xaxis, cell_cell_knn_yaxis,
        common_genes,
        X, Y_cg, Y_mcg, 
        modx_clsts, knn_xy, 
        enhancer_gene_to_eval,
        output_corrs,
    ):
    """
    """
    # new cells  
    common_rna_cells_updated = np.intersect1d(common_rna_cells, cell_cell_knn_xaxis)
    common_mc_cells_updated = np.intersect1d(common_mc_cells, cell_cell_knn_yaxis)

    # make sure the original matrices have the correct index
    x_idx = snmcseq_utils.get_index_from_array(common_rna_cells, common_rna_cells_updated)
    y_idx = snmcseq_utils.get_index_from_array(common_mc_cells, common_mc_cells_updated)
    X = X.tocsc()[:, x_idx] 
    Y_cg = Y_cg.tocsc()[:, y_idx]
    Y_mcg = Y_mcg.tocsc()[:, y_idx] 

    # make sure knn_xy, knn_xx have the right cell index
    cell_idx_xaxis = snmcseq_utils.get_index_from_array(cell_cell_knn_xaxis, common_rna_cells_updated)
    cell_idx_yaxis = snmcseq_utils.get_index_from_array(cell_cell_knn_yaxis, common_mc_cells_updated)
    knn_xy = knn_xy.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_yaxis] # x-by-y
    modx_clsts = modx_clsts.reindex(common_rna_cells_updated)

    logging.info("{}_{}_{}_{}_{}".format(knn_xy.shape, modx_clsts.shape, X.shape, Y_cg.shape, Y_mcg.shape))

    for clst_col in modx_clsts.columns: 
        print(clst_col)

        # choose one clustering to proceed
        uniq_labels = np.sort(modx_clsts[clst_col].unique()) 
        print("Number of metacells: {}".format(len(uniq_labels)))

        knn_xz = enhancer_gene_utils.turn_cluster_labels_to_knn(modx_clsts[clst_col].values, 
                                            uniq_labels,
                                           )
        # normalization - such that metacells made of more cells still sums to 1
        knn_xz = knn_xz.dot(sparse.diags(np.ravel(1.0/knn_xz.sum(axis=0))))

        # gene by metacell
        gc_rna = X.dot(knn_xz).todense() 

        # enhancer by metacell
        knn_yz = knn_xy.T.dot(knn_xz)
        ec_cg = Y_cg.dot(knn_yz).todense() 
        ec_mcg = Y_mcg.dot(knn_yz).todense()  
        print(gc_rna.shape, ec_cg.shape, ec_mcg.shape)

        # mC
        ec_mccg = snmcseq_utils.get_mcc_lite_v4(
                                       pd.DataFrame(ec_cg).astype(np.float32), 
                                       pd.DataFrame(ec_mcg).astype(np.float32), 
                                       base_call_cutoff=5, sufficient_coverage_fraction=0.8, fillna=True)
        print(ec_mccg.shape)

        output_corr = output_corrs.format(clst_col)
        # corr analysis
        (to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells) = enhancer_gene_utils.compute_enh_gene_corrs(
            gc_rna, ec_mccg, 
            common_genes, ec_mccg.index.values,
            enhancer_gene_to_eval['gene'].values, 
            enhancer_gene_to_eval['ens'].values, 
            output_file=output_corr, chunksize=100000, verbose_level=0,
            )
    return

def wrap_corr_analysis_mc(
        mod_x, mod_y, 
        input_nme_tag, i_sub,
    ):
    """
    """
    # (i, k, --r)

    # input enh-gene tables, gene-by-cell, enhancer-by-cell matrices
    input_enh_gene_table = '/cndd2/fangming/projects/scf_enhancers/results/200521_to_evals.tsv' 
    input_bundle_dirc = '/cndd2/fangming/projects/scf_enhancers/data/organized_cell_level/version_nov9'
    bundle_fnames = (
        'cell_10x_cells_v3.txt',
        'cell_snmcseq_gene.txt',

        'gene_10x_cells_v3.txt',
        'enh.tsv',

        'mat_10x_cells_v3.npz',
        'mat_mcg_snmcseq_gene.npz',
        'mat_cg_snmcseq_gene.npz',
    )

    # for cross modal knn
    # name = "corr_analysis_{}_{}_{}".format(mod_x, mod_y, input_name_tag) 
    output_corrs = '/cndd2/fangming/projects/scf_enhancers/results/{}_{}_{{}}_corrs.pkl'.format(input_name_tag, i_sub)

    # for knn_xx
    input_knn_dirc = '/cndd2/fangming/projects/miniatlas/results'
    input_modx_clsts = [
        'clusterings_{}_{}_sub{}.tsv.gz'.format(mod_x, input_name_tag, i_sub),
    ]

    # for knn_xy
    input_knn_xy = 'knn_across_{}_{}_{}.npz.{}.npz'.format(input_name_tag, mod_x, mod_y, i_sub) 
    input_knn_cells_xaxis = 'cells_{}_{}.npy.{}.npy'.format(mod_x, input_name_tag, i_sub)
    input_knn_cells_yaxis = 'cells_{}_{}.npy.{}.npy'.format(mod_y, input_name_tag, i_sub)

    # # Load data 
    # input_bundle
    with snmcseq_utils.cd(input_bundle_dirc):
        bundle = []
        for fname in bundle_fnames:
            #  save all as pickle file
            with open(fname, "rb") as fh:
                item = pickle.load(fh)
            bundle.append(item)
            logging.info("{}_{}_{}".format(type(item), item.shape, fname))

    (common_rna_cells, common_mc_cells, 
     common_genes, common_enhancer_regions,
     X, Y_mcg, Y_cg, 
    #  knn_xy, knn_xx,
    ) = bundle

    # input knn networks 
    with snmcseq_utils.cd(input_knn_dirc):
        # for knn_xx 
        # modx_clsts = pd.read_csv(input_modx_clsts, sep='\t',index_col=0)
        modx_clsts = pd.concat([
            pd.read_csv(fname, sep='\t',index_col=0)
            for fname in input_modx_clsts
        ], axis=1)
        # for knn_xy 
        knn_xy = sparse.load_npz(input_knn_xy)  
        cell_cell_knn_xaxis = np.load(input_knn_cells_xaxis, allow_pickle=True)
        cell_cell_knn_yaxis = np.load(input_knn_cells_yaxis, allow_pickle=True)

        print(modx_clsts.shape, knn_xy.shape, 
              cell_cell_knn_xaxis.shape, 
              cell_cell_knn_yaxis.shape,
             )

    # enhancer-gene linkage
    enhancer_gene_to_eval = pd.read_csv(input_enh_gene_table, sep='\t')

    pipe_corr_analysis_mc(
        common_rna_cells, common_mc_cells,
        cell_cell_knn_xaxis, cell_cell_knn_yaxis,
        common_genes,
        X, Y_cg, Y_mcg, 
        modx_clsts, knn_xy, 
        enhancer_gene_to_eval,
        output_corrs,
    )
    return 

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-modx', '--mod_x', required=True, type=str)
    parser.add_argument('-mody', '--mod_y', required=True, type=str)
    parser.add_argument('-tag', '--input_name_tag', required=True, type=str)
    # parser.add_argument('-knn', '--knn_across', required=True, type=int)
    parser.add_argument('-isub', '--i_sub', type=str, help="[0-9]")
    return parser


if __name__ == "__main__":
    # 
    parser = create_parser()
    args = parser.parse_args()

    # output setting
    # run this with each combination of (i_sub, knn)
    mod_x = args.mod_x
    mod_y = args.mod_y
    input_name_tag = args.input_name_tag
    i_sub = args.i_sub

    print(mod_x, mod_y, input_name_tag, i_sub)

    # run
    wrap_corr_analysis_mc(
        mod_x, mod_y, 
        input_name_tag, i_sub,
    )

