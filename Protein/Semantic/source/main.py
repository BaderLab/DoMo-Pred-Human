'''
Created on 2010-07-19

@author: Shobhit Jain

@contact: shobhit@cs.toronto.edu
'''
from ontology import GOGraph
from cPickle import dump, load
# from dill import dump, load
import os


def save_semantic_similarity(ontology_file, gene_file, ontology, code, gpath):
    '''
    Calls functions for loading and processing data.
    '''
    objs = {}
    ontology = ontology.split(",")
    g = GOGraph()
    write = False
    onto_v = g._obo_parser(ontology_file)
    anno_v = g._go_annotations(gene_file, code)
    run = {'C': g._cellular_component, 'P': g._biological_process,
           'F': g._molecular_function}
    ont = {'C': "cc", 'P': "bp", 'F': "mf"}
    for i in ontology:
        fd = ont[i[0]] + ".pck"
        if fd in os.listdir(gpath):
            prompt = raw_input('Do you want to over-write existing ' + fd +
                               ' file? (y/n): ')
            if prompt == 'y':
                write = True
        else:
            write = True

        if write:
            objs_file = open(gpath + "/%s" % (fd), "w")
            i = i.split(":")
            print 'Run started...'
            objs[i[0]] = run[i[0]]()
            print 'Species...'
            objs[i[0]]._species()
            print 'Clustering...'
            objs[i[0]]._clustering(float(i[1]))
            print 'Writing file...'
            dump(objs[i[0]], objs_file)

    return objs


def load_semantic_similarity(CC, BP, MF):
    '''
    Calls functions for loading and processing data.
    '''
    objs = {}
    objs['C'] = load(open(CC))
    objs['P'] = load(open(BP))
    objs['F'] = load(open(MF))

    return objs


def return_details(result, geneA, geneB, detail):
    '''
    Formats the output for printing on screen or on file.
    '''
    domain_def = {'C':'Cellular Component', 'P':'Biological Process', 'F':'Molecular Function'}
    r = "\nSemantic similarity between " + geneA + " and " + geneB + " is:\n\n"
    for domain in result:
        r += " " + domain_def[domain] + ": " + str(result[domain][0]) + "\n"
        if detail:
            for data in result[domain][1]:
                r += "  GO id assigned to " + geneA + " is: " + data[0] + \
                     "\n  GO id assigned to " + geneB + " is: " + data[1] + \
                     "\n  LCA of assigned GO ids is: " + "|".join(result[domain][1][data]['lca']) + "\n\n"
    return r + "\n\n\n"


def return_detail(result, geneA, geneB, detail):
    '''
    '''
    for domain in result:
        r = geneA + '\t' + geneB + '\t' + str(result[domain][0]) + '\n'
    return r


def calculate_semantic_similarity(objs, geneA, geneB):
    '''
    Calls the function for calculating semantic similarity between
    genesA and genesB.
    '''

    return objs._semantic_similarity(geneA, geneB)


if __name__ == '__main__':
    path = '/Users/shobhit/Project/Yeast-SH3-Project/DoMo-Pred/Protein'
    objs = save_semantic_similarity(path + "/Db/gene_ontology.1_2.obo.txt", path +\
            "/Db/gene_association.sgd", "C:2.4,P:3.5,F:3.3", "", path + '/Db/')
