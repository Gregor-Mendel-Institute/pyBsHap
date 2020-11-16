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
