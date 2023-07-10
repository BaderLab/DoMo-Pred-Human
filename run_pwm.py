import os
import sys
sys.path.append('Classifier/')
# import matplotlib.pyplot as plt
from numpy import loadtxt, array, column_stack, savetxt

os.environ["WORK_DIR_PEP"] = os.getcwd() + '/Peptide/'
os.environ["WORK_DIR_PRO"] = os.getcwd() + '/Protein/'
os.environ["WORK_DIR_CLS"] = os.getcwd() + '/Classifier/'

from Peptide import run_peptide
from Protein import run_protein
from Classifier import semi_nb


def load_datasets():
    '''
    '''
    gen = run_peptide.process_genome_file()
    pep = run_peptide.setup_peptide()
    pot = run_protein.setup_protein()
    seq = dict(loadtxt('domain.txt', delimiter='\t', dtype='S'))

    return gen, pep, pot, seq


def peptide_features(pwm, domain, genome, pval, data):
    '''
    '''
    res = run_peptide.run_pwm(pwm, domain, genome, pval, data[0], data[1],
                              data[2], data[3])

    return res


def protein_features(prot_set, data):
    '''
    '''
    res = run_protein.run_features(prot_set, data[0]["C"],
                                   data[0]["P"], data[0]["F"],
                                   data[1], data[2], data[3], data[4])

    return res


def load_classifiers(pep_path, pot_path):
    '''
    '''
    pep_modl = semi_nb.load_model(pep_path)
    pot_modl = semi_nb.load_model(pot_path)

    return pep_modl, pot_modl


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


def calculate_posterior(data, model):
    '''
    '''
    result = []
    for case in data:
        label, prob = model.posterior(case)
        result.append({label: prob, 1 - label: 1 - prob})

    return result


def save_text(file_path, data):
    '''
    '''
    head = 'Domain\tPeptide\tStart\tEnd\tSequence\tPeptide_prob\tProtein_prob\tProbability'
    resl = []
    for i in data:
        # print i
        if float(i[-1]) >= 0.5:
            resl.append(i)
    resl = array(resl)
    savetxt(file_path, resl, delimiter="\t", fmt="%s", header=head)



if __name__ == '__main__':
     '''
     '''
     MODL = 'Classifier/saved_models/'
     gen, pep, pot, seq = load_datasets()
     pep_modl, pot_modl = load_classifiers(MODL + 'peptide_unlbl_cls_2019-03-18.pck',
                                           MODL + 'protein_unlbl_cls_2019-03-18.pck')
 
     for pwm in os.listdir('SH3/'):
         if pwm[:-4] + '.res' in os.listdir('OUT/'):
             continue
         print 'Working with PWM: ', pwm, '...................'
         dom = pwm[:-4].split('--')[0].replace('_', '/')
         
         for cutoff in [1e-05, 1e-04, 1e-03]:
             print 'Running peptide features with PWM cutoff: ', cutoff
             peptide, pep_info = peptide_features('SH3/' + pwm, seq[dom], gen,
                                                  cutoff, pep)
             print 'No. of peptide hits: ', len(peptide)
             if len(peptide) > 0:
                 break
 
         if len(peptide) == 0:
             pep_post = []
             pot_post = []
         else:
             print 'Calculating posterior prob for peptide classifier'
             pep_post = calculate_posterior(peptide[:, 1:], pep_modl)
 
             print 'Running protein features'
             protein, pot_info = protein_features(pep_info, pot)
             print 'Calculate posterior prob for protein classifier'
             pot_post = calculate_posterior(protein, pot_modl)
 
         print 'Calculating final probabilities'
         com_prob = []
         for i, j in zip(pep_post, pot_post):
             prob = bayes(i, j)[0]
             com_prob.append([i[0], j[0], prob])
 
         com_prob = array(com_prob)
 
         print 'Saving output to file'
         save_text('OUT/' + pwm[:-4] + '.res',
                   column_stack((pep_info, com_prob)))
         print('.......................................................')
