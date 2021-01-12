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


# In[7]:


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

    ka = args.knn_within 
    knn = args.knn_across 
    result_date = args.result_date
    i_sub = args.i_sub

    today = datetime.date.today()
    name = "corr_analysis_withinmc_ka{}_knn{}_{}_{}".format(ka, knn, result_date, today) 
    logging.info(name)

    input_knn_xx = '/cndd2/fangming/projects/miniatlas/results/knn_within_mop_2mods_{}_ka{}_knn{}_snmcseq_gene.npz.{}.npz'.format(result_date, ka, knn, i_sub) 
    input_knn_cells_xaxis = '/cndd2/fangming/projects/miniatlas/results/cells_{}_mop_2mods_{}_ka{}_knn{}.npy.{}.npy'.format('snmcseq_gene', result_date, ka, knn, i_sub)

    input_enh_gene_table = '/cndd2/fangming/projects/scf_enhancers/results/200521_to_evals.tsv' 

    input_bundle_dirc = '/cndd2/fangming/projects/scf_enhancers/data/organized_cell_level/version_mc_only_aug9'
    bundle_fnames = (
        'cell_snmcseq_gene.txt',
        'gene_snmcseq_gene.txt',
        'enh_snmcseq_gene.tsv',
        
        'mat_genebody_mch_snmcseq_gene.npz',
        'mat_genebody_ch_snmcseq_gene.npz',
        'mat_enhancer_mcg_snmcseq_gene.npz',
        'mat_enhancer_cg_snmcseq_gene.npz',
    )


    output_corr = '/cndd2/fangming/projects/scf_enhancers/results/{}_{}_corrs.pkl'.format(name, i_sub)
    output_to_correlate = '/cndd2/fangming/projects/scf_enhancers/results/{}_{}_corrs_idx.pkl'.format(name, i_sub)
    output_figure1 = '/cndd2/fangming/projects/scf_enhancers/results/figures/{}_{}_plot1.pdf'.format(name, i_sub)
    output_figure2 = '/cndd2/fangming/projects/scf_enhancers/results/figures/{}_{}_plot1_v2.pdf'.format(name, i_sub)
    logging.info("Output: {}\n{}\n{}\n{}".format(output_corr, output_to_correlate, output_figure1, output_figure2))

    # # Load data 

    # In[4]:


    # input_bundle
    with snmcseq_utils.cd(input_bundle_dirc):
        bundle = []
        for fname in bundle_fnames:
            #  save all as pickle file
            with open(fname, "rb") as fh:
                item = pickle.load(fh)
            bundle.append(item)
            logging.info("{}_{}_{}".format(type(item), item.shape, fname))

    (common_cells,  
     common_genes, common_enhancer_regions,
     X_mch, X_ch, Y_mcg, Y_cg, 
    ) = bundle


    # In[7]:
    # knn_xy, knn_xx
    # knn_xy = sparse.load_npz(input_knn_xy)  
    knn_xx = sparse.load_npz(input_knn_xx) 

    # import their axes
    cell_cell_knn_xaxis = np.load(input_knn_cells_xaxis, allow_pickle=True)
    # cell_cell_knn_yaxis = np.load(input_knn_cells_yaxis)

    # new cells  
    common_cells_updated = np.intersect1d(common_cells, cell_cell_knn_xaxis)

    # make sure the original matrices have the correct index
    x_idx = snmcseq_utils.get_index_from_array(common_cells, common_cells_updated)
    X_ch = X_ch.tocsc()[:, x_idx] 
    X_mch = X_mch.tocsc()[:, x_idx] 
    Y_cg = Y_cg.tocsc()[:, x_idx]
    Y_mcg = Y_mcg.tocsc()[:, x_idx] 

    # make sure knn_xy, knn_xx have the right cell index
    cell_idx_xaxis = snmcseq_utils.get_index_from_array(cell_cell_knn_xaxis, common_cells_updated)
    knn_xx = knn_xx.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_xaxis] # x-by-x

    logging.info("{}_{}_{}_{}_{}".format(knn_xx.shape, X_ch.shape, X_mch.shape, Y_cg.shape, Y_mcg.shape))




    # In[ ]:


    # enhancer-gene linkage
    enhancer_gene_to_eval = pd.read_csv(input_enh_gene_table, sep='\t')

    # # Compute metacell level signal 

    # In[5]:

    # gene by metacell
    knn_xx.data = [1]*len(knn_xx.data) # remove all weights
    gc_ch = X_ch.dot(knn_xx.T).todense() 
    gc_mch = X_mch.dot(knn_xx.T).todense() 


    # enhancer by metacell
    ec_cg = Y_cg.dot(knn_xx.T).todense() 
    ec_mcg = Y_mcg.dot(knn_xx.T).todense()  

    logging.info("{}_{}_{}".format(gc_ch.shape, gc_mch.shape, ec_cg.shape, ec_mcg.shape))


    # In[6]:


    # get mcc
    ec_mccg =  snmcseq_utils.get_mcc_lite_v4(
                                    pd.DataFrame(ec_cg).astype(np.float32), 
                                    pd.DataFrame(ec_mcg).astype(np.float32), 
                                    base_call_cutoff=5, sufficient_coverage_fraction=0.8, fillna=True)
    gc_mcch =  snmcseq_utils.get_mcc_lite_v4(
                                    pd.DataFrame(gc_ch).astype(np.float32), 
                                    pd.DataFrame(gc_mch).astype(np.float32), 
                                    base_call_cutoff=100, sufficient_coverage_fraction=0.8, fillna=True)
    logging.info("{} {}".format(ec_mccg.shape, gc_mcch.shape))

    # ## Correlate enhancer and gene (a new notebook)

    # In[8]:


    gc_mcch_zscore = stats.zscore(gc_mcch.values, axis=1, ddof=0,)
    ec_mccg_zscore = stats.zscore(ec_mccg.values, axis=1, ddof=0,)


    # In[9]:


    # use ec_mccg and gc_rna
    # correlate e-g according to a e-g table
    gene_idx = snmcseq_utils.get_index_from_array(common_genes[gc_mcch.index.values], enhancer_gene_to_eval['gene'])
    enh_idx = snmcseq_utils.get_index_from_array(ec_mccg.index.values, enhancer_gene_to_eval['ens']) # be careful here!
    to_correlate = ~np.logical_or(gene_idx==-1, enh_idx==-1)
    with open(output_to_correlate, "wb") as fh:
        pickle.dump(to_correlate, fh)
    
    gene_idx = gene_idx[to_correlate]
    enh_idx = enh_idx[to_correlate]

    # corr
    corrs = row_dot_product_norm_by_numcol(gc_mcch_zscore, ec_mccg_zscore, gene_idx, enh_idx)

    # corr shuffled cells
    corrs_shuffled_cells = row_dot_product_norm_by_numcol(
        gc_mcch_zscore[:,np.random.permutation(gc_mcch_zscore.shape[1])], 
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
        gc_mcch_zscore, 
        ec_mccg_zscore, 
        gene_idx_shuff, enh_idx)

    # In[10]:


    # save ec_mccg and gc_rna
    with open(output_corr, 'wb') as fh:
        pickle.dump((corrs, corrs_shuffled, corrs_shuffled_cells), fh)


    # # plotting results

    # In[11]:


    with open(output_corr, 'rb') as fh:
        corrs, corrs_shuffled, corrs_shuffled_cells = pickle.load(fh)
    logging.info("{}_{}_{}".format(corrs.shape, corrs_shuffled.shape, corrs_shuffled_cells.shape))


    # In[12]:


    dists = enhancer_gene_to_eval.loc[to_correlate, 'dist'].values
    logging.info("{}_{}".format(np.min(dists), np.max(dists)))


    # In[13]:


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


    # In[14]:


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

