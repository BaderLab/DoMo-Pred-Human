import sys
sys.path.append("Classifier/")
from semi_nb import *
import numpy as np
import cPickle


def load_model(filename):
    '''
    '''
    with open(filename, 'rb') as input:
        model = cPickle.load(input)
    return model


if __name__ == '__main__':
    '''
    '''
    # choose the trained model
    classifier = load_model('saved_models/protein_lbl_cls_2018-10-31.pck')
    data = np.loadtxt('../test.pot', usecols=(2, 3, 4, 5, 6, 7, 8))
    names = np.loadtxt('../test.pot', usecols=(0, 1), dtype='str')

    # label for positive class is 0
    probabilities = []
    for case in data:
        clas, prob = classifier.posterior(case)
        if clas == 0:
            probabilities.append(prob)
        else:
            probabilities.append(1 - prob)

    results = np.column_stack([names, probabilities])

    # save the results
    np.savetxt('pdz_results.csv', results, fmt='%s', delimiter=',')



