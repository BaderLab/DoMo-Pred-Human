from scipy.stats import pearsonr
import os
import math
import numpy as np

matrix = {}

def load_expression(path):
    '''
    '''
    if os.path.exists(path + '/profile_binary.npz'):
        return np.load(path + '/profile_binary.npz')

    return False


def _fisher_transform(lst):
    '''
    '''
    for i in range(len(lst)):
        if lst[i] == 1.0:
            lst[i] = lst[i] - 0.0000000000000001
        lst[i] = 0.5*(math.log(1 + lst[i]) - math.log(1 - lst[i]))
    return lst



def _fisher_mean(lst):
    '''
    '''
    return sum(lst) / len(lst)


def _inversion(fisher_mean):
    '''
    '''
    return (math.exp(2 * fisher_mean) - 1) / (math.exp(2 * fisher_mean) + 1)


def _calculate_corr(matA, matB):
    '''
    number of datasets is 140
    '''
    pc = []
    for i in range(140):
        if not np.isnan(matA[i][0]) and not np.isnan(matB[i][0]):
            p = pearsonr(matA[i], matB[i])[0]
            if not np.isnan(p):
                pc.append(p)
    if len(pc) != 0:
        return _inversion(_fisher_mean(_fisher_transform(pc)))
    return np.nan


def expression(a, b, profile):

    if a == b:
        return np.nan
    else:
        if a in profile and b in profile:
            apc = _calculate_corr(profile[a], profile[b])
            if not np.isnan(apc):
                apc = 1 / (1 + math.e**(-(5 * apc)))
                # apc = apc
            return apc
        else:
            return np.nan
