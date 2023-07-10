''' Load protein pipeline '''

from Semantic import *
from Expression import *
from Sequence import *
from ProteinExp import *
from Topology import *
from numpy import array

PATH = os.getenv("WORK_DIR_PRO") + 'Db'

def setup_protein():
    '''
    Setup feature datasets.
    '''
    semt, expn, pexp, sign, topo = (None, None, None, None, None)

    print "\nSetting up protein features .... \n"
    print "  Setting up cellular location .... "
    print "  Setting up biological process .... "
    print "  Setting up molecular function .... "
    semt = load_semantic(PATH + '/GO')

    print "  Setting up expression ...."
    expn = load_expression(PATH + '/EX')

    print "  Setting up protein expression ...."
    pexp = load_proteinx(PATH + '/PX')

    print "  Setting up sequence signature ...."
    sign = load_sequence(PATH + '/SS')

    print "  Setting up topology ....\n\n"
    topo = load_topology(PATH + '/NT')


    return semt, expn, pexp, sign, topo


def run_features(prot_set, *args):
    '''
    fetaure data
    '''

    features = [semantic_score, semantic_score, semantic_score, expression, pexpression, sequence, topology]
    defaults = [None, None, None, None, None, None, None]

    results = []
    int_set = []
    count = 0
    for prot1, prot2, s, e, seq in prot_set:
        # print count
        count += 1
        int_set.append([prot1, prot2, s, e, seq])
        if prot1 == 'None' or prot2 == 'None':
            results.append(defaults)
        else:
            tmp = []
            for idx, data in enumerate(args):
                if data:
                    # print(features[idx])
                    val = features[idx](prot1, prot2, data)
                    tmp.append(val)
            results.append(tmp)

    return array(results), array(int_set)
