import numpy as np
import pandas as pd
import collections
import tqdm

def kth_diag_indices(n, k=0):
    """
    """
    assert k >= 0
    
    row_idx = np.arange(n)
    col_idx = np.arange(n)
    
    col_idx = row_idx + k
    
    if k > 0:
        return row_idx[:-k], col_idx[:-k]
    if k == 0:
        return row_idx, col_idx

def mat_values_kth_diag_indices(mat, k=0):
    """
    """
    return mat[kth_diag_indices(len(mat), k=k)]

def enh_gene_id_to_binidx(_table, tsss, enhs, resolution):
    """
    Given a dataframe of enh-gene pairs (by IDs)
    return which bins in HiC data they correspond to
    
    Return a dataframe with dist_bin as index, 'gene_bin', 'enh_bin' as columns
    """
    # merge gid, enh_id with genome locations
    _tmp = pd.merge(_table, tsss.rename(columns={'chr': 'gene_chr',
                                                 'start': 'gene_start',
                                                 'end': 'gene_end',
                                                 'gid': 'gene',
                                                }), on='gene')

    _tmp = pd.merge(_tmp, enhs.rename(columns={'chr': 'enh_chr',
                                               'start': 'enh_start',
                                               'end': 'enh_end',
                                               'center': 'enh_center',
                                               'enh_id': 'enh',
                                              }), on='enh')
    # find the correct bins
    _tmp['enh_bin'] = (_tmp['enh_center'].values/resolution).astype(int)
    _tmp['gene_bin'] = (_tmp['gene_start'].values/resolution).astype(int)
    _tmp['dist_bin'] = np.abs(_tmp['enh_bin'] - _tmp['gene_bin'])
    _tmp = _tmp.set_index('dist_bin')
    
    return _tmp

def get_contact_stats(mat, paired_bin_tables, resolution, distance_cutoff=1e6):
    """
    """
    distance_idx = np.arange(1, int(distance_cutoff/resolution))
    distances = distance_idx * resolution

    contacts_mean = []
    contacts_std = []
    contacts_sig_pairs_mean = collections.OrderedDict({_label: [] for _label, _table in paired_bin_tables.items()})
    contacts_sig_pairs_std = collections.OrderedDict({_label: [] for _label, _table in paired_bin_tables.items()})

    for i in tqdm.tqdm(distance_idx):
        contacts = mat_values_kth_diag_indices(mat, k=i)
        contacts_mean.append(np.nanmean(contacts)) 
        contacts_std.append(np.nanstd(contacts)) 

        for _label, _table in paired_bin_tables.items():
            try:
                pairs = _table.loc[i]
                contacts_sig_pairs = mat[(pairs['enh_bin'].values, pairs['gene_bin'].values)]
            except:
                contacts_sig_pairs = []

            contacts_sig_pairs_mean[_label].append(np.nanmean(contacts_sig_pairs)) 
            contacts_sig_pairs_std[_label].append(np.nanstd(contacts_sig_pairs)) 

    contacts_mean = np.array(contacts_mean)
    contacts_std = np.array(contacts_std)

    contacts_sig_pairs_mean = collections.OrderedDict({
                    _label: np.array(_table) for _label, _table in contacts_sig_pairs_mean.items()})
    contacts_sig_pairs_std = collections.OrderedDict({
                    _label: np.array(_table) for _label, _table in contacts_sig_pairs_std.items()})

    return (contacts_mean, contacts_std, 
            contacts_sig_pairs_mean, contacts_sig_pairs_std,
           )

def get_contacts(mat, paired_bin_tables, resolution, distance_idx):
    """
    given distances/distance idx; 
    - contacts: a list of numpy arrays
    - contacts_sig_pairs: a dict of list of numpy arrays 
    """

    contacts = []
    contacts_sig_pairs = collections.OrderedDict({_label: [] for _label, _table in paired_bin_tables.items()})

    for i in tqdm.tqdm(distance_idx):
        _contacts = mat_values_kth_diag_indices(mat, k=i) # get all contacts from only that distance bin
        contacts.append(_contacts) 

        for _label, _table in paired_bin_tables.items():
            try:
                pairs = _table.loc[i] # get pairs from only that distance bin
                _contacts_sig_pairs = mat[(pairs['enh_bin'].values, pairs['gene_bin'].values)]
            except:
                _contacts_sig_pairs = []
            contacts_sig_pairs[_label].append(_contacts_sig_pairs) 

    return (contacts, contacts_sig_pairs)

