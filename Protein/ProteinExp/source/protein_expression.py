import numpy as np
from scipy.stats.mstats import pearsonr
import math
from cPickle import dump, load


def save_protein_expression(path, out_path):
    '''
    '''

    fd = open(path)
    fd.readline()
    data = {}
    for line in fd:
        line = line.strip().split(",")
        pnts = map(float, line[1:])
        pid = line[0]
        if pnts.count(0.0) < 27:
            pnts = np.array(pnts)
            # pnts = np.ma.masked_where(pnts == 0.0, pnts)
            if pid not in data:
                data[pid] = np.array([pnts])
            else:
                val = np.array([pnts])
                data[pid] = np.row_stack((data[pid], val))

    with open(out_path + '/protein_exp.pck', 'wb') as out_strm:
        dump(data, out_strm)


def load_proteinx(path):
    '''
    '''

    with open(path + '/protein_exp.pck', 'rb') as in_strm:
        data = load(in_strm)

    return data


def pexpression(a, b, data):
    '''
    '''
    if a == b:
        return np.nan

    if a in data and b in data:
        tmp = []
        for row1 in data[a]:
            for row2 in data[b]:
                tmp.append(pearsonr(row1, row2)[0])
            if max(tmp) is not np.ma.masked:
                pc = 1 / (1 + math.e**(-(5 * max(tmp))))
                # pc = 1 - pc
                return pc

    return np.nan
