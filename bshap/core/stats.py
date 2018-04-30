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

def get_random_permeth_diffs(mx, my, req_inds, samplenos, nmc=1000):
    ## Cannot take very huge number of permutations.
    rand_context_inds = np.sort( np.random.choice(req_inds, nmc * samplenos) )
    rand_diffs_context_inds = np.abs(mx.get_permeths(mx.filter_pos_ix[rand_context_inds]) - my.get_permeths(my.filter_pos_ix[rand_context_inds]))
    np.random.shuffle(rand_diffs_context_inds)
    return(np.mean( rand_diffs_context_inds.reshape((nmc, samplenos)), axis = 1 ))
