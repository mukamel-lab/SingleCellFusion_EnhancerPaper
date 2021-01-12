"""HiC related notebook functions"""
import numpy as np
import matplotlib.pyplot as plt
# import multiprocessing as mp
# import h5py
# import pickle
import seaborn as sns
import tqdm

from scipy import stats

from statsmodels.stats.multitest import multipletests

import sys
# import itertools
# import time

sys.path.insert(0, '/cndd/fangming/CEMBA/snmcseq_dev')
from __init__ import *
from __init__jupyterlab import * 
import snmcseq_utils

from enhancer_gene_hic_validation_utils import *
from hic_plot_utils import *


def ttest_routine(contacts_mean, contacts_std, contacts_n, 
                  contacts_sig_pairs_mean,
                  contacts_sig_pairs_std,
                  contacts_sig_pairs_n,
                  p_th1=0.05, p_th2=0.001):
    """
    """
    # vs background
    mean_fcs_bck = collections.OrderedDict()
    padjs_bck = collections.OrderedDict()
    for key, item in contacts_sig_pairs_mean.items():
        # mean FC
        mean_fcs_bck[key] = item/contacts_mean

        # significance t-test
        t, p = stats.ttest_ind_from_stats(contacts_mean, contacts_std, contacts_n, 
                                          contacts_sig_pairs_mean[key],
                                          contacts_sig_pairs_std[key],
                                          contacts_sig_pairs_n[key],
                                          equal_var=True)
        # multiple comparison 
        _, padj, _, _ = multipletests(p, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)

        # 
        padjs_bck[key] = padj

    # vs correlated
    mean_fcs_correlated = collections.OrderedDict()
    padjs_correlated = collections.OrderedDict()
    for key, item in contacts_sig_pairs_mean.items():
        if key.startswith('linked_'):
            # mean FC
            key_matched = key.replace('linked_', 'correlated_') 
            mean_fcs_correlated[key] = item/contacts_sig_pairs_mean[key_matched]

            # significance t-test
            t, p = stats.ttest_ind_from_stats( 
                                              contacts_sig_pairs_mean[key_matched],
                                              contacts_sig_pairs_std[key_matched],
                                              contacts_sig_pairs_n[key_matched],
                                              contacts_sig_pairs_mean[key],
                                              contacts_sig_pairs_std[key],
                                              contacts_sig_pairs_n[key],
                                              equal_var=True)
            # multiple comparison 
            _, padj, _, _ = multipletests(p, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)

            # 
            padjs_correlated[key] = padj
            
    return mean_fcs_bck, padjs_bck, mean_fcs_correlated, padjs_correlated
    
def aggregate_mean_std(mu, sigma, n):
    """Given a list of mu, sigma, n 
    return the aggregated mu, sigma, n    
    """
    mu = np.array(mu)
    sigma = np.array(sigma)
    n = np.array(n)
    
    # add robustness to n=0, mu and sigma=np.nan
    cond = (n!=0)
    n = n[cond].copy()
    mu = mu[cond].copy()
    sigma = sigma[cond].copy()
    
    n_agg = np.sum(n)
    
    ratio = n/n_agg
    mu_agg = np.dot(ratio, mu)
    sigma_agg = np.sqrt(np.dot(ratio, np.power(sigma, 2)) 
                      + np.dot(ratio, np.power((mu-mu_agg), 2))
                       )
    return mu_agg, sigma_agg, n_agg

def aggregate_mean_std_matrix(mu, sigma, n):
    """Given a matrix of mu, sigma, n  
    return the aggregated array mu, sigma, n (remove the first matrix axis) 
    """
    mu = np.array(mu)
    sigma = np.array(sigma)
    n = np.array(n)
    nrow, ncol = mu.shape
    
    mu_agg_all = []
    sigma_agg_all = []
    n_agg_all = []
    for col_idx in np.arange(ncol):
        mu_agg, sigma_agg, n_agg = aggregate_mean_std(mu[:,col_idx], sigma[:,col_idx], n[:,col_idx])
        mu_agg_all.append(mu_agg)
        sigma_agg_all.append(sigma_agg)
        n_agg_all.append(n_agg)
        
    return np.array(mu_agg_all), np.array(sigma_agg_all), np.array(n_agg_all)

def plot2(distances, track_names, mean_fc_all, celltype_palette, output_fig, 
          ylim=[-0.7, 0.7], yticks=[-0.5, 0, 0.5],
         ):
    """
    """
    fig, axs = plt.subplots(2, 3, figsize=(5*3,5*2), sharex=True,)
    for idx, row in mean_fc_all.iterrows():
        celltype = row['celltype']
        for ax, track_name, in zip(axs.flat, track_names,):
            ax.plot(distances, np.clip(np.log2(row['mean_fc'][track_name]), -1, 1), 
                    label=celltype, color=celltype_palette[celltype],
                   )
            
    for ax, track_name, in zip(axs.flat, track_names,):
        ax.xaxis.set_major_formatter(mtick.EngFormatter())
        ax.set_xlabel('Genomic distance')
        ax.set_xlim([0, 1e5])
        ax.set_ylim(ylim)
        ax.set_yticks(yticks)

        ax.set_title(track_name)
        ax.axhline(0, linestyle='--', color='gray')

    axs.flat[0].set_ylabel('log2 FC')
    ax.legend(bbox_to_anchor=(1,1), loc='upper left')

    # fig.suptitle("", fontsize=20)
    snmcseq_utils.savefig(fig, output_fig)
    plt.show() 
    return 

def plot3(distances, mean_fc_all, 
          track_names, track_palette, 
          p_adjs,
          output_fig, 
          ylim=[-0.3, 0.6], 
          yticks=[-0.3,0,0.3,0.6],
          p_th1=0.05,
          p_th2=0.001,
         ):
    """
    """
    
    ys = {track_name: [] for track_name in track_names}
    for idx, row in mean_fc_all.iterrows():
        celltype = row['celltype']
        for track_name in track_names:
            ys[track_name].append(
                np.clip(np.log2(row['mean_fc'][track_name]), -1, 1)
                )

    for track_name in track_names:
        ys[track_name] = np.array(ys[track_name])

    ys_mean = {track_name: np.nanmean(ys[track_name], axis=0) for track_name in track_names}
    ys_std = {track_name: np.nanstd(ys[track_name], axis=0) for track_name in track_names}
    ys_err = {track_name: ys_std[track_name]*1.96/np.sqrt(8) for track_name in track_names}

    p_adjs = {'mc': [], 'atac': [], 'both': [],}

    # t test compare linked vs correlated
    num_celltypes, num_dists = ys['linked_mc'].shape
    for catg in ['mc', 'atac', 'both']:
        for i in np.arange(num_dists):
            _a = ys['linked_{}'.format(catg)][:, i]
            _b = ys['correlated_{}'.format(catg)][:, i]
            _t, _pval = stats.ttest_rel(_a, _b)

            # multiple comparison 
            _, p_adj, _, _ = multipletests(_pval, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)

            # 
            p_adjs[catg].append(p_adj)


    fig, axs = plt.subplots(1, 3, figsize=(5*3,5), sharey=True, sharex=True)
    for i, (ax, _type) in enumerate(zip(axs, ['mc', 'atac', 'both'])):
        for track_name in track_names:
            color = track_palette[track_name]
            if track_name.endswith(_type):
                pvals = p_adjs[_type]
                for idx, dist in enumerate(distances):
                    if pvals[idx] < p_th2:
                        ax.text(dist, -0.3, '*\n*\n*', fontsize=15, linespacing=0.3, ha='center')
                    elif pvals[idx] < p_th1: 
                        ax.text(dist, -0.3, '*', fontsize=15, ha='center')

                ax.fill_between(distances, 
                        ys_mean[track_name]-ys_err[track_name], 
                        ys_mean[track_name]+ys_err[track_name], 
                        color=color,
                        alpha=0.2
                       )
                ax.plot(distances, ys_mean[track_name], 
                        label=track_name, color=color,
                        linewidth=5,
                       )
        ax.axhline(0, linestyle='--', color='gray')
        
        ax.set_title(_type)
        ax.set_xlabel('Genomic distance')
        if i == 0:
            ax.set_ylabel('log2(FC)\n(+/- 95% CI; n=8 cell types)')
            ax.xaxis.set_major_formatter(mtick.EngFormatter())
            ax.set_xlim([-5e3, 1.05*1e5])
            ax.set_ylim(ylim)
            ax.set_yticks(yticks)

    axs[2].annotate("*: FDR < 0.05\n***: FDR<0.001", (1.05, 0.1), xycoords='axes fraction', fontsize=15)
    handles, labels = snmcseq_utils.combine_legends(axs.flat)
    # handles, labels = snmcseq_utils.dedup_legends(handles, labels)
    ax.legend(handles, labels, bbox_to_anchor=(1,1), loc='upper left')

    snmcseq_utils.savefig(fig, output_fig)
    plt.show() 


def plot1_v3(distances, celltype, resolution, 
          track_names, track_palette, 
          contacts_mean, contacts_std, contacts_n, 
          contacts_sig_pairs_mean, contacts_sig_pairs_std, contacts_sig_pairs_n,
          padjs_correlated,
          output_fig,
          pval_y=0.003,
          ylim=[3e-3, 1e-1],
          yticks=[3e-3, 1e-2, 1e-1],
          p_th1=0.05,
          p_th2=0.001,
         ):
    """
    """
    # plot 
    fig, axs = plt.subplots(1, 3, figsize=(5*3,5*1), sharey=True, sharex=True)
    for i, (_name) in enumerate([_name for _name in track_names if _name.startswith('linked')]):
        if _name.endswith('_mc'):
            ax = axs[0]
        elif _name.endswith('_atac'):
            ax = axs[1]
        elif _name.endswith('_both'):
            ax = axs[2]
        _name_matched = _name.replace('linked', 'correlated')
            
        # define tracks
        track_info = [
            (distances, contacts_mean, contacts_std*1.96/np.sqrt(contacts_n), 
             'gray', 'all_bins', '--', 3),
            
            (distances, contacts_sig_pairs_mean[_name_matched], 
             contacts_sig_pairs_std[_name_matched]*1.96/np.sqrt(contacts_sig_pairs_n[_name_matched]), 
             'black', _name_matched, '-', 5),
            
            (distances, contacts_sig_pairs_mean[_name], 
             contacts_sig_pairs_std[_name]*1.96/np.sqrt(contacts_sig_pairs_n[_name]), 
             track_palette[_name], _name, '-', 5),
        ]
        for (_x, _y, _yerr, _color, _label, _linestyle, _linewidth) in track_info:
            ax.fill_between(_x, _y-_yerr, _y+_yerr, color=_color, alpha=0.2)
            ax.plot(_x, _y, _linestyle, color=_color, linewidth=_linewidth, label=_label, alpha=1)
        
        # significance
        pvals = padjs_correlated[_name]
        for idx, dist in enumerate(distances):
            if pvals[idx] < p_th2:
                ax.text(dist, pval_y, '*\n*\n*', fontsize=15, linespacing=0.3, ha='center')
            elif pvals[idx] < p_th1: 
                ax.text(dist, pval_y, '*', fontsize=15, ha='center')
        if i == 0:
            ax.set_ylabel('Contact frequency\n(mean +/- 95% CI)')
        ax.set_xlabel('Genomic distance')
        
        ax.set_title(_name.split('_')[-1].upper())
        
    ax.set_yscale('log')
    ax.set_xlim([0, 1e5])
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)
    ax.xaxis.set_major_formatter(mtick.EngFormatter())
    ax.set_xlabel('Genomic distance')

    # legends
    handles, labels = snmcseq_utils.combine_legends(axs) #.flat)
    handles, labels = snmcseq_utils.dedup_legends(handles, labels)
    ax.legend(handles, labels, bbox_to_anchor=(1,1), loc='upper left')
    ax.annotate("linked vs correlated pairs\n*: FDR < 0.05\n***: FDR<0.001", (1.05, 0.1), xycoords='axes fraction', fontsize=15)

    for ax in axs:
        sns.despine(ax=ax)
    
    fig.suptitle("{} {}, {} HiC resolution".format('All chroms', celltype, resolution), y=1.02)

    snmcseq_utils.savefig(fig, output_fig)
    plt.show()
    return


def plot1_v3_cov(distances, celltype, resolution, 
          track_names, track_palette, 
          contacts_mean, contacts_std, contacts_n, 
          contacts_sig_pairs_mean, contacts_sig_pairs_std, contacts_sig_pairs_n,
          output_fig,
         ):
    """
    """
    # plot 
    fig, axs = plt.subplots(1, 3, figsize=(5*3,5*1), sharey=True, sharex=True)
    for i, (_name) in enumerate([_name for _name in track_names if _name.startswith('linked')]):
        if _name.endswith('_mc'):
            ax = axs[0]
        elif _name.endswith('_atac'):
            ax = axs[1]
        elif _name.endswith('_both'):
            ax = axs[2]
        _name_matched = _name.replace('linked', 'correlated')
            
        # define tracks
        track_info = [
            (distances, contacts_n, 
             'gray', 'all_bins', '--', 3),
            
            (distances, contacts_sig_pairs_n[_name_matched], 
             'black', _name_matched, '-', 5),
            
            (distances, contacts_sig_pairs_n[_name], 
             track_palette[_name], _name, '-', 5),
        ]
        for (_x, _y, _color, _label, _linestyle, _linewidth) in track_info:
            ax.plot(_x, _y, _linestyle, color=_color, linewidth=_linewidth, label=_label, alpha=1)
        
        if i == 0:
            ax.set_ylabel('Number of pairs fall in the distance bin\n(10k resolution)')
        ax.set_xlabel('Genomic distance')
        
        ax.set_title(_name.split('_')[-1].upper())
        
    ax.set_yscale('log')
    ax.set_xlim([0, 1e5])
    ax.xaxis.set_major_formatter(mtick.EngFormatter())
    ax.yaxis.set_major_formatter(mtick.EngFormatter())
    ax.set_xlabel('Genomic distance')

    # legends
    handles, labels = snmcseq_utils.combine_legends(axs) #.flat)
    handles, labels = snmcseq_utils.dedup_legends(handles, labels)
    ax.legend(handles, labels, bbox_to_anchor=(1,1), loc='upper left')
    
    for ax in axs:
        sns.despine(ax=ax)

    fig.suptitle("{} {}, {} HiC resolution".format('All chroms', celltype, resolution), y=1.02)

    snmcseq_utils.savefig(fig, output_fig)
    plt.show()
    return


def plot1_v4(distances, celltype, resolution, 
          track_names, track_palette,
          contacts_mean, contacts_std, contacts_n, 
          contacts_sig_pairs_mean, contacts_sig_pairs_std, contacts_sig_pairs_n,
          padjs_correlated,
          output_fig,
          ylim=[-0.3, 0.6],
          yticks=[-0.3,0,0.3,0.6],
          pval_y=-0.4,
          p_th1=0.05,
          p_th2=0.001,
         ):
    """
    """
    # plot 
    fig, axs = plt.subplots(1, 3, figsize=(5*3,5*1), sharey=True, sharex=True)
    for i, (_name) in enumerate([_name for _name in track_names if _name.startswith('linked')]):
        if _name.endswith('_mc'):
            ax = axs[0]
        elif _name.endswith('_atac'):
            ax = axs[1]
        elif _name.endswith('_both'):
            ax = axs[2]
        _name_matched = _name.replace('linked', 'correlated')
            
        # define tracks
        track_info = [
#             (distances, contacts_mean, contacts_std*1.96/np.sqrt(contacts_n), 
#              'gray', 'all_bins', '--', 3),
            
            (distances, contacts_sig_pairs_mean[_name_matched], 
             contacts_sig_pairs_std[_name_matched]*1.96/np.sqrt(contacts_sig_pairs_n[_name_matched]), 
             'black', _name_matched, '-', 5),
            
            (distances, contacts_sig_pairs_mean[_name], 
             contacts_sig_pairs_std[_name]*1.96/np.sqrt(contacts_sig_pairs_n[_name]), 
             track_palette[_name], _name, '-', 5),
        ]
        for (_x, _y, _yerr, _color, _label, _linestyle, _linewidth) in track_info:
            ax.fill_between(_x, np.log2((_y-_yerr)/contacts_mean), np.log2((_y+_yerr)/contacts_mean), color=_color, alpha=0.2)
            ax.plot(_x, np.log2(_y/contacts_mean), _linestyle, color=_color, linewidth=_linewidth, label=_label, alpha=1)
        
        # significance
        pvals = padjs_correlated[_name]
        for idx, dist in enumerate(distances):
            if pvals[idx] < p_th2:
                ax.text(dist, pval_y, '*\n*\n*', fontsize=15, linespacing=0.3, ha='center')
            elif pvals[idx] < p_th1: 
                ax.text(dist, pval_y, '*', fontsize=15, ha='center')
                
        # labels
        if i == 0:
            ax.set_ylabel('log2(Fold change in contact frequency)\n(mean +/- 95% CI)')
        ax.axhline(0, linestyle='--', color='gray')
        ax.set_xlabel('Genomic distance')
        ax.set_title(_name.split('_')[-1].upper())
    
    ax.set_xlim([-5e3, 1.05*1e5])
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)

    ax.xaxis.set_major_formatter(mtick.EngFormatter())
    ax.set_xlabel('Genomic distance')

    # legends
    handles, labels = snmcseq_utils.combine_legends(axs) #.flat)
    handles, labels = snmcseq_utils.dedup_legends(handles, labels)
    ax.legend(handles, labels, bbox_to_anchor=(1,1), loc='upper left')
    ax.annotate("linked vs correlated pairs\n*: FDR < 0.05\n***: FDR<0.001", (1.05, 0.1), xycoords='axes fraction', fontsize=15)

    for ax in axs:
        sns.despine(ax=ax)
    
    fig.suptitle("{} {}, {} HiC resolution".format('All chroms', celltype, resolution), y=1.02)

    snmcseq_utils.savefig(fig, output_fig)
    plt.show()
    return
