''' Load genomic pipeline '''

from Disordered import *
from Surface import *
from Peptide import *
from PWMsearch import *
from Structure import *
from numpy import loadtxt, log2, array, savetxt, column_stack, string_

PATH = os.getenv("WORK_DIR_PEP") + 'Db'


def process_genome_file(file_name='human_protein.fasta'):
    '''process the protein file and generate a dictionary of sequences'''

    file = open(PATH + '/SQ/' + file_name)
    uid, seq, data = ("", "", {})
    for line in file:
        if line.startswith(">"):
            if uid in data:
                data[uid] = seq
                seq = ""
            uid = line.strip().strip(">")
            if uid not in data:
                data[uid] = ""
            else:
                print "duplicate %s" % uid
        else:
            seq += line.strip()
    data[uid] = seq
    return data


def setup_peptide():
    '''
    Setup feature datasets.
    '''
    disd, surf, pepd, strt = (None, None, None, None)

    print "\nSetting up peptide features .... \n"
    print "   Setting up disorder .... "
    disd = load_disorder(PATH + '/DR')

    print "   Setting up surface .... "
    surf = load_surface(PATH + '/SA')

    print "   Setting up peptide .... "
    pepd = load_peptide(PATH + '/PC')

    print "   Setting up structure .... "
    strt = load_structure(PATH + '/SC')

    return disd, surf, pepd, strt


# def information_content(column):
#    '''Return information cotent of column'''
#
#    return abs(sum(column * log(column)))

def before_entropy(background):
    '''
    '''
    return -sum(background * log2(background))


def after_entropy(column):
    '''
    '''
    return -sum(column * log2(column))


def read_pwms(pwm_path, bck_path, edge_threshold=0.5, inner_threshold=0.4, extra=False):
    '''
    Return PWM as array and its significant positions
    Note: made changes in entropy
    '''

    barray = loadtxt(bck_path)
    farray = loadtxt(pwm_path, delimiter='\t', skiprows=1, dtype=str)
    column = farray[:, 0]
    farray = farray[:, 1:-1] # extra space in pwm file at end of line 
    farray = farray.astype(float)
    farray = farray / 20
    #farray = farray[:, 3: -1]
    #return farray, range(10)

    background = before_entropy(barray)

    ic = []
    for i in range(farray.shape[1]):
        ic.append(background - after_entropy(farray[:, i]))


    max_ic = max(ic)
    ic = array(ic) / max_ic

    i = 0
    while ic[i] < edge_threshold:
        i += 1

    j = 0
    len_ic = len(ic)
    while ic[len_ic - 1 - j] < edge_threshold:
        j += 1

    #end = farray.shape[1] + j + 1

    trimmed = farray[:, i: (len_ic - j)]
    if trimmed.shape[1] < 4:
        trimmed = farray
    else:
        ic = ic[i: (len_ic - j)]

    pos = []
    for idx, val in enumerate(ic):
        if val >= inner_threshold:
            pos.append(idx)

    if extra:
        sig = '-'.join(map(str, pos))
        savetxt('../trimmed_pwm/%s-(%s)'%(pwm_path.split('/')[-1], sig), column_stack((column, trimmed)), fmt='%s')

    return trimmed, tuple(pos)


def run_features(prot_set, *args):
    '''
    '''

    features = [disorder_score, surface_score, peptide_score, structure_score]

    results = []
    int_set = []

    for prot1, prot2, dseq, pseq, start, end, pos in prot_set:
        tmp = []
        if type(pos) is str or type(pos) is string_: 
            pos = map(int, pos.split(','))
        for idx, data in enumerate(args):
            if data:
                if idx == 3:
                    val = features[idx](dseq, pseq, data[0], data[1], data[2])
                else:
                    val = features[idx](prot2, int(start), int(end), data, pos)
                tmp.append(val)

        results.append(tmp)
        int_set.append([prot1, prot2, int(start) + 1, int(end), pseq])

    return array(results), array(int_set)


def run_pwm(pwm_file, dseq, genome, tol, *args):
    '''
    '''

    features = [disorder_score, surface_score, peptide_score, structure_score]

    pwm, pos = read_pwms(pwm_file, PATH + '/SQ/human.background')
    # print pwm, pos
    initialize_database(PATH + '/SQ/human_protein.pwm',
                        PATH + '/SQ/human.background')
    length = initialize_pwm(pwm)
    tol = get_threshold(tol)
    hits = search_database(tol)
    # print hits
    pwm_name = open(pwm_file).readline().strip('#').strip()

    results = []
    int_set = []
    for hit in hits:
        query = hit
        for start, end in hits[hit]:
            score = hits[hit][(start, end)] / 100.0
            pseq = genome[query][start:end]

            tmp = [score]
            for idx, data in enumerate(args):
                if data:
                    if idx == 3:
                        val = features[idx](dseq, pseq, data[0], data[1], data[2])
                    else:
                        val = features[idx](query, start, end, data, pos)
                    tmp.append(val)

            results.append(tmp)
            int_set.append([pwm_name, query, start + 1, end, pseq])

    results = array(results)
    if results != []:
        max_scr = max(results[:, 0])
        results[:, 0] /= max_scr ##### CHECK #####

    # print results
    return results, array(int_set)


def pwm_file_search(pwm_file, pwm_name, genome, tol):
    '''
    '''
    initialize_database(PATH + '/SQ/dk.pwm',
                        PATH + '/SQ/human.background')
    length = initialize_pwm_file(pwm_file)
    print length
    tol = get_threshold(tol)
    hits = search_database(tol)
    int_set = []
    results = []
    for hit in hits:
        query = hit
        for start, end in hits[hit]:
            score = hits[hit][(start, end)] / 100.0
            pseq = genome[query][start:end]
            results.append([score])
            int_set.append([pwm_name, query, start + 1, end, pseq, score])
    # results = array(results)
    # max_scr = max(results[:, 0])
    # print results
    return array(int_set)


def pwm_search(pwm_file, pwm_name, genome, tol):
    '''
    '''
    pwm, pos = read_pwms(pwm_file, PATH + '/SQ/human.background')
    initialize_database(PATH + '/SQ/dk.pwm',
                        PATH + '/SQ/human.background')
    length = initialize_pwm(pwm)
    print length
    tol = get_threshold(tol)
    hits = search_database(tol)
    int_set = []
    results = []
    for hit in hits:
        query = hit
        for start, end in hits[hit]:
            score = hits[hit][(start, end)] / 100.0
            pseq = genome[query][start:end]
            results.append([score])
            int_set.append([pwm_name, query, start + 1, end, pseq, score])
    # results = array(results)
    # max_scr = max(results[:, 0])
    # print results
    return array(int_set)




if __name__ == '__main__':

    for pwm in os.listdir('../pwm_dir/'):
        print pwm
        read_pwms('../pwm_dir/' + pwm, 'Db/yeast.background')
