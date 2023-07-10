import os, sys
from numpy import loadtxt, nan, column_stack, row_stack, array, savetxt

os.environ["WORK_DIR_CLS"] = os.getcwd() + '/Classifier/'

from Classifier import semi_nb


def load_classifiers(pep_path, pot_path):
    '''
    '''
    pep_modl = semi_nb.load_model(pep_path)
    pot_modl = semi_nb.load_model(pot_path)

    return pep_modl, pot_modl


def calculate_posterior(data, model):
    '''
    '''
    result = []
    for case in data:
        label, prob = model.posterior(case)
        result.append({label: prob, 1 - label: 1 - prob})

    return result


def load_dataset(path, col):
    '''
    '''
    data = loadtxt(path, dtype=str)
    if len(data.shape) == 1:
        data = array([data])
    data[data == 'None'] = nan
    return data[:, :col], data[:, col:].astype(float)


def bayes(prob1, prob2):
    '''
    '''

    tmp = {}
    for class_label in prob1:
        tmp[class_label] = (prob1[class_label] * prob2[class_label]) / 0.5

    # print prob1, prob2
    den = sum(tmp.values())
    for class_label in tmp:
        tmp[class_label] /= den

    return tmp


def save_text(file_path, data):
    '''
    '''
    savetxt(file_path, data, delimiter="\t", fmt="%s")


if __name__ == '__main__':
    '''
    '''
    MODL = 'Classifier/'
    PATH = '/Users/shobhit/Project/Mohamed/'
    pep_modl, pot_modl = load_classifiers(MODL + 'peptide.pck',
                                          MODL + 'prot_negatome.pck')

    alph = 'abcdefghijklmnopqrstuvwxyz'
    for char1 in alph:
        for char2 in alph:
            fname = 'x' + char1 + char2
            print fname
            if fname + '.pep' not in os.listdir(PATH):
                sys.exit(0)
            if fname + '.res' in os.listdir(PATH):
                print fname, '  done!'
                continue
            pep_info, pep_data = load_dataset(PATH + fname + '.pep', -4)
            pep_post = calculate_posterior(pep_data, pep_modl)

            pot_info, pot_data = load_dataset(PATH + fname + '.pot', -6)
            pot_post = calculate_posterior(pot_data, pot_modl)

            com_prob = []
            for i, j in zip(pep_post, pot_post):
                prob = bayes(i, j)[0]
                com_prob.append([i[0], j[0], prob])

            com_prob = array(com_prob)
            save_text(PATH + fname + '.res',
                      column_stack((pep_info, com_prob)))
