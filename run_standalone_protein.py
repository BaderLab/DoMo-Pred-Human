''' Code below is for running standards and generating histogram plots'''

import os, sys
import matplotlib.pyplot as plt
from numpy import loadtxt, column_stack, chararray, random, savetxt, arange

os.environ["WORK_DIR_PRO"] = os.getcwd() + '/Protein/'

from Protein import run_protein


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


def run_standard(file_path, data):
    '''
    '''
    data_set = loadtxt(file_path, dtype=str)
    # prot_set = data_set[:, -4: -2]
    prot_set = data_set
    empty = chararray((len(prot_set), 3))
    empty[:] = ''
    prot_set = column_stack((prot_set, empty))
    print(prot_set)

    protein, pro_set = run_protein.run_features(prot_set, data[0]["C"],
                                                data[0]["P"], data[0]["F"],
                                                data[1], data[2], data[3],
                                                data[4])

    assert data_set.shape[0] == protein.shape[0], 'input/output mismatch'
    return protein, data_set


if __name__ == '__main__':
    '''
    '''
    data = run_protein.setup_protein()
    print 'Working on data...'
    file_data, pos_set = run_standard('test.txt', data)
    save_text(file_data, pos_set, 'test' + '.pot')

