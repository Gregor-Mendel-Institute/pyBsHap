# Summary statistics
import numpy as np
from scipy import stats
import logging
import h5py as h5
import pandas as pd


log = logging.getLogger(__name__)

def exact_mc_perm_test(xm, ym, nmc):
    ### Exact permutation test by randomizing
    ## xm and ym are numpy arrays. nmc is exact number of tests
    n, k = len(xm), 0
    diff = np.mean(np.abs(xm - ym))
    zs = np.concatenate([xm, ym])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.mean(np.abs(zs[:n] - zs[n:]))
    return(float(k) / nmc)

def get_random_permeth_diffs(meths_x, meths_y, req_inds, samplenos, nmc=1000):
    ## Cannot take very huge number of permutations.
    ## Here meths_x and meths_y have filter_pos_ix
    rand_context_inds = np.sort( np.random.choice(req_inds, nmc * samplenos) )
    rand_diffs_context_inds = np.abs(meths_x.get_permeths(meths_x.filter_pos_ix[rand_context_inds]) - meths_y.get_permeths(meths_y.filter_pos_ix[rand_context_inds]))
    np.random.shuffle(rand_diffs_context_inds)
    return(np.mean( rand_diffs_context_inds.reshape((nmc, samplenos)), axis = 1 ))

def differential_methylation(meths_x, meths_y, req_inds):
    ## Here meths_x and meths_y have filter_pos_ix
    import statsmodels.stats.api as sms
    permeths_x = meths_x.get_permeths(meths_x.filter_pos_ix[req_inds])
    permeths_y = meths_y.get_permeths(meths_y.filter_pos_ix[req_inds])
    cm = sms.CompareMeans(sms.DescrStatsW(permeths_x), sms.DescrStatsW(permeths_y))
    return(cm.ttest_ind())


def pivot_mCs_required_bed(mc_bed, req_start_df, min_bp = 5000):
    """
    Function to take mean at each position based on its distance from another start file
    """
    assert type(mc_bed) is pd.DataFrame, "provide a dataframe"
    assert type(req_start_df) is pd.DataFrame, "provide a dataframe"
    assert mc_bed.shape[1] == 2, "provide a dataframe with index and mC/phenotype"
    mc_start_mat = np.tile( mc_bed.iloc[:,0].values, (req_start_df.shape[0], 1) ).T
    mc_pheno_mat = np.tile( mc_bed.iloc[:,1].values, (req_start_df.shape[0], 1) ).T
    req_bed_mat = np.tile( req_start_df.iloc[:,0].values, (mc_bed.shape[0], 1) )
    start_diff_mat = mc_start_mat - req_bed_mat
    req_positions_ix = np.where(np.abs(start_diff_mat) < min_bp)

    dist_positions = start_diff_mat[req_positions_ix]
    pheno_positions = mc_pheno_mat[req_positions_ix]

    sorted_by_dist_ix = np.argsort(dist_positions)
    output_df = pd.DataFrame( np.column_stack([dist_positions[sorted_by_dist_ix], pheno_positions[sorted_by_dist_ix]]), columns = ['dist', 'pheno'] )
    output_df = output_df.groupby('dist', sort = False).mean()
    return( output_df )