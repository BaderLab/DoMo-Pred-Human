''' Code below is for running standards and generating histogram plots'''

import os
import matplotlib.pyplot as plt
from numpy import loadtxt, column_stack, chararray, random, savetxt, arange

os.environ["WORK_DIR_PEP"] = os.getcwd() + '/Peptide/'

from Peptide import run_peptide


def filter(array_a, array_b):
    '''
    '''

    array_x = array_a[array_a > None]
    array_y = array_b[array_b > None]

    len_x = len(array_x)
    len_y = len(array_y)

    random.shuffle(array_x)
    random.shuffle(array_y)

    if len_x > len_y:
        return array_x[:len_y], array_y
    else:
        return array_x, array_y[:len_x]


def save_text(data, int_set, file_path):
    '''
    '''
    data = column_stack((int_set, data))
    savetxt(file_path, data, delimiter="\t", fmt="%s")


def transform_data(pwm, pep, dom, pos, pseq, dseq):
    '''
    '''
    result = []
    for i in range(len(pwm)):
        pwm_id = pwm[i].split('--')[0].replace('_', '/')
        st, ed = pos[i].split(':')
        result.append([dom[i], pep[i], dseq[pwm_id], pseq[i], st, ed, range(int(ed) - int(st))])
    return result


def run_standard(file_path, data, domain):
    '''
    '''
    data_set = loadtxt(file_path, dtype=str)
    pept_set = transform_data(data_set[:, -7], data_set[:, -4], data_set[:, -3],
                              data_set[:, -1], data_set[:, 2], domain)

    peptide, pept_set = run_peptide.run_features(pept_set, data[0], data[1],
                                                data[2], data[3])

    assert data_set.shape[0] == peptide.shape[0], 'input/output mismatch'
    return peptide, data_set


if __name__ == '__main__':
    data = run_peptide.setup_peptide()
    domain = dict(loadtxt('domain.txt', dtype='str'))
    print 'Working on data...'

    path = '/Users/shobhit/Project/Mohamed/'
    alph = 'abcdefghijklmnopqrstuvwxyz'
    for char1 in alph:
        for char2 in alph:
            fname = 'x' + char1 + char2
            print fname
            if fname + '.pep' in os.listdir(path):
                print fname,  '  done!'
                continue
            file_data, pos_set = run_standard(path + fname, data, domain)
            save_text(file_data, pos_set, path + fname + '.pep')

