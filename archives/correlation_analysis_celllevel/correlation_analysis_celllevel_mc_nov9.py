#!/usr/bin/env python3

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


def row_dot_product_norm_by_numcol(X_zscore, Y_zscore, x_idx, y_idx, 
                                   chunksize=10000, verbose_level=100000):
    """compute (X_zscore[x_idx]*Y_zscore[y_idx]).mean(axis=1)
    correlation values given matched x_idx and y_idx...
    """
    assert len(x_idx) == len(y_idx)
    num_pairs = len(x_idx)
    corrs = []
    for pair_idx in snmcseq_utils.chunks(np.arange(num_pairs), chunksize):
        if pair_idx[0] % verbose_level == 0:
            logging.info(pair_idx[0])

        _res = (X_zscore[x_idx[pair_idx]]*Y_zscore[y_idx[pair_idx]]).mean(axis=1)
        corrs.append(_res)
    corrs = np.hstack(corrs) 
    return corrs 

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-tag', '--input_name_tag', required=True, type=str)
    parser.add_argument('-modx', '--mod_x', required=True, type=str)
    parser.add_argument('-mody', '--mod_y', required=True, type=str)

    parser.add_argument('-knn', '--knn_across', required=True, type=int)
    parser.add_argument('-ka', '--knn_within', required=True, type=int)
    parser.add_argument('-dt', '--result_date', type=str, help="eg: 200803")
    parser.add_argument('-isub', '--i_sub', type=str, help="[0-9]")
    return parser


if __name__ == '__main__':
    ti = time.time()

    parser = create_parser()
    args = parser.parse_args()
    logger = snmcseq_utils.create_logger()

    input_name_tag = args.input_name_tag
    mod_x = args.mod_x
    mod_y = args.mod_y

    knn = args.knn_across 
    ka = args.knn_within 
    result_date = args.result_date
    i_sub = args.i_sub

    name = "corr_analysis_{}_{}_{}".format(mod_x, mod_y, input_name_tag) 
    logging.info(name)
    output_corr = '/cndd2/fangming/projects/scf_enhancers/results/{}_{}_corrs.pkl'.format(name, i_sub)
    output_to_correlate = '/cndd2/fangming/projects/scf_enhancers/results/{}_{}_corrs_idx.pkl'.format(name, i_sub)
    output_figure1 = '/cndd2/fangming/projects/scf_enhancers/results/figures/{}_{}_plot1.pdf'.format(name, i_sub)
    output_figure2 = '/cndd2/fangming/projects/scf_enhancers/results/figures/{}_{}_plot1_v2.pdf'.format(name, i_sub)

    input_knn_xy = '/cndd2/fangming/projects/miniatlas/results/knn_across_{}_{}_{}.npz.{}.npz'.format(input_name_tag, mod_x, mod_y, i_sub)
    input_knn_xx = '/cndd2/fangming/projects/miniatlas/results/knn_within_{}_{}.npz.{}.npz'.format(input_name_tag, mod_x, i_sub) 
    input_knn_cells_xaxis = '/cndd2/fangming/projects/miniatlas/results/cells_{}_{}.npy.{}.npy'.format(mod_x, input_name_tag, i_sub)
    input_knn_cells_yaxis = '/cndd2/fangming/projects/miniatlas/results/cells_{}_{}.npy.{}.npy'.format(mod_y, input_name_tag, i_sub)

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

    # knn_xy, knn_xx
    knn_xy = sparse.load_npz(input_knn_xy)  
    knn_xx = sparse.load_npz(input_knn_xx) 

    # import their axes
    cell_cell_knn_xaxis = np.load(input_knn_cells_xaxis, allow_pickle=True)
    cell_cell_knn_yaxis = np.load(input_knn_cells_yaxis, allow_pickle=True)

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
    knn_xx = knn_xx.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_xaxis] # x-by-x

    logging.info("{}_{}_{}_{}_{}".format(knn_xy.shape, knn_xx.shape, X.shape, Y_cg.shape, Y_mcg.shape))

    # enhancer-gene linkage
    enhancer_gene_to_eval = pd.read_csv(input_enh_gene_table, sep='\t')

    # # Compute metacell level signal 
    # random sampling
    if n_metacell_samples > 0:
        # sample n indices from a list ...
        np.arange(len(common_rna_cells))


    # gene by metacell
    gc_rna = X.dot(knn_xx.T).todense() 
    # enhancer by metacell
    ec_cg = Y_cg.dot(knn_xy.T).todense() 
    ec_mcg = Y_mcg.dot(knn_xy.T).todense()  

    logging.info("{}_{}_{}".format(gc_rna.shape, ec_cg.shape, ec_mcg.shape))

    # get mcc
    ec_mccg =  snmcseq_utils.get_mcc_lite_v4(
                                    pd.DataFrame(ec_cg).astype(np.float32), 
                                    pd.DataFrame(ec_mcg).astype(np.float32), 
                                    base_call_cutoff=5, sufficient_coverage_fraction=0.8, fillna=True)
    logging.info(ec_mccg.shape)

    # ## Correlate enhancer and gene (a new notebook)
    gc_rna_zscore = stats.zscore(np.array(gc_rna), axis=1, ddof=0,)
    ec_mccg_zscore = stats.zscore(ec_mccg.values, axis=1, ddof=0,)

    # use ec_mccg and gc_rna
    # correlate e-g according to a e-g table
    gene_idx = snmcseq_utils.get_index_from_array(common_genes, enhancer_gene_to_eval['gene'])
    enh_idx = snmcseq_utils.get_index_from_array(ec_mccg.index.values, enhancer_gene_to_eval['ens']) # be careful here!
    to_correlate = ~np.logical_or(gene_idx==-1, enh_idx==-1)
    with open(output_to_correlate, "wb") as fh:
        pickle.dump(to_correlate, fh)
    
    gene_idx = gene_idx[to_correlate]
    enh_idx = enh_idx[to_correlate]

    # corr
    corrs = row_dot_product_norm_by_numcol(gc_rna_zscore, ec_mccg_zscore, gene_idx, enh_idx)

    # corr shuffled cells
    corrs_shuffled_cells = row_dot_product_norm_by_numcol(
        gc_rna_zscore[:,np.random.permutation(gc_rna_zscore.shape[1])], 
        ec_mccg_zscore, 
        gene_idx, enh_idx)

    # corr shuffled genes (break up the pairs)
    gene_idx_uniq = np.unique(gene_idx)
    shuff_genes = {
        gene: gene_shuff for gene, gene_shuff in 
            zip(gene_idx_uniq, gene_idx_uniq[np.random.permutation(len(gene_idx_uniq))])
        }
    gene_idx_shuff = np.array([shuff_genes[gene] for gene in gene_idx])

    corrs_shuffled = row_dot_product_norm_by_numcol(
        gc_rna_zscore, 
        ec_mccg_zscore, 
        gene_idx_shuff, enh_idx)

    # save ec_mccg and gc_rna
    with open(output_corr, 'wb') as fh:
        pickle.dump((corrs, corrs_shuffled, corrs_shuffled_cells), fh)


    # # plotting results
    with open(output_corr, 'rb') as fh:
        corrs, corrs_shuffled, corrs_shuffled_cells = pickle.load(fh)
    logging.info("{}_{}_{}".format(corrs.shape, corrs_shuffled.shape, corrs_shuffled_cells.shape))

    dists = enhancer_gene_to_eval.loc[to_correlate, 'dist'].values
    logging.info("{}_{}".format(np.min(dists), np.max(dists)))

    config = {
        'kde': False,
        "hist_kws": {
                    'histtype': 'step', 
    #                 'edgecolor': 'none',
                    'alpha': 1, 
                    'density': True, 
                    },
    }

    colors = snmcseq_utils.get_grad_colors(5, cmap='Blues_r')
    tracks = {
        'pairs (<100kb)': corrs[dists<1e5], 
        'pairs (<500kb)': corrs[dists<5e5], 
        'pairs (<1Mb)': corrs, 
        'shuffled pairs': corrs_shuffled, 
    #     'shuffled cells': corrs_shuffled_cells, 
        }
    track_colors = {
        'pairs (<100kb)': colors[0], 
        'pairs (<500kb)': colors[1], 
        'pairs (<1Mb)': colors[2], 
        'shuffled pairs': 'gray', 
    #     'shuffled cells': 'gray', 
        }

    # figure 1
    num_bins = 200
    bins = np.linspace(-0.3, 0.3, num_bins)

    tracks_hist_ratios = {}
    track_cdfs = {}
    fdr_cdfs = {}

    hist_shuff, _ = np.histogram(corrs_shuffled, bins=bins, density=True)
    cdf_shuff = np.cumsum(hist_shuff)
    for label, track in tracks.items():
        hist, _ = np.histogram(track, bins=bins, density=True)
        cdf = np.cumsum(hist)
        fdr = cdf_shuff/cdf
        tracks_hist_ratios[label] = hist/hist_shuff
        track_cdfs[label] = cdf 
        fdr_cdfs[label] = fdr

    fig, axs = plt.subplots(3, 1, figsize=(4*2,4*3), sharex=True)

    ax = axs[0]
    for label, track in tracks.items():
        sns.distplot(track, bins=bins, ax=ax, label=label, color=track_colors[label], **config)
    ax.legend(bbox_to_anchor=(1,1), loc='upper left')
    ax.set_ylabel('fraction of pairs')

    ax = axs[1]
    for label, track in tracks_hist_ratios.items():
        ax.plot(bins[1:], track, label="data/shuffled", color=track_colors[label])
    ax.axhline(1, linestyle='--', color='gray')
    ax.set_ylabel('ratio (data/shuffled)')

    ax = axs[2]
    for label, track in fdr_cdfs.items():
        ax.plot(bins[1:], track, color=track_colors[label])
    ax.set_ylim([0, 1.1])
    ax.set_ylabel('FDR')
    ax.set_xlabel('Pearson r')

    snmcseq_utils.savefig(fig, output_figure1)
    plt.show()

    # figure 2
    num_bins = 200
    bins = np.linspace(-1, 0.3, num_bins)

    tracks_hist_ratios = {}
    track_cdfs = {}
    fdr_cdfs = {}

    hist_shuff, _ = np.histogram(corrs_shuffled, bins=bins, density=True)
    cdf_shuff = np.cumsum(hist_shuff)
    for label, track in tracks.items():
        hist, _ = np.histogram(track, bins=bins, density=True)
        cdf = np.cumsum(hist)
        fdr = cdf_shuff/cdf
        tracks_hist_ratios[label] = hist/hist_shuff
        track_cdfs[label] = cdf 
        fdr_cdfs[label] = fdr
    fig, axs = plt.subplots(3, 1, figsize=(4*2,4*3), sharex=True)

    ax = axs[0]
    for label, track in tracks.items():
        sns.distplot(track, bins=bins, ax=ax, label=label, color=track_colors[label], **config)
    ax.legend(bbox_to_anchor=(1,1), loc='upper left')
    ax.set_ylabel('fraction of pairs')

    ax = axs[1]
    for label, track in tracks_hist_ratios.items():
        ax.plot(bins[1:], track, label="data/shuffled", color=track_colors[label])
    ax.set_yscale('log')
    ax.axhline(1, linestyle='--', color='gray')
    ax.set_ylabel('ratio (data/shuffled)')

    ax = axs[2]
    for label, track in fdr_cdfs.items():
        ax.plot(bins[1:], track, color=track_colors[label])
    ax.set_ylim([1e-5, 2])
    ax.set_yscale('log')
    ax.set_ylabel('FDR')
    ax.set_xlabel('Pearson r')

    snmcseq_utils.savefig(fig, output_figure2)
    plt.show()

