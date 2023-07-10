import os, sys, shutil
import numpy as np


def max_sample_size(path1, path2):
    '''
    '''
    length = []
    for fname in os.listdir(path1):
        samples = len(open(path1 + fname).readline().strip("\n").split("\t"))
        if samples < 301:
            length.append(samples)
            shutil.move(path1 + fname, path2 + fname)
    return max(length)


def mapping(path):
    '''
    '''
    fdata = open(path).readlines()
    map_data = {}
    for line in fdata[1:]:
        line = line.strip().split()
        if line[0] not in map_data:
            map_data[line[0]] = set()
        map_data[line[0]].add(line[1])
    return map_data


def temporary_files(in_path, entrez_unip, out_path):
    '''
    '''
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    else:
        shutil.rmtree(out_path)
        os.makedirs(out_path)

    i = 0
    for fname in os.listdir(in_path):
        print i, fname
        fdata = open(in_path + fname)
        for line in fdata:
            line = line.strip().split("\t")
            if line[0] in entrez_unip:
                for unip in entrez_unip[line[0]]:
                    resl = open(out_path + unip, "a+")
                    resl.write(str(i) + "\t" + "\t".join(line[1:]) + "\n")
                    resl.close()
        i += 1
        fdata.close()


def process(data, file, cnt, sample_size):
    line = []
    for i in range(len(data[0])):
        col = []
        count = 0
        for j in range(len(data)):
            try:
                col.append(float(data[j][i]))
            except:
                pass
        if len(col) != 0:
            line.append(sum(col) / len(col))
        else:
            print cnt, file
            return [np.nan] * sample_size
    return line



def create_numpy_files(ipath, opath, sample_size):
    '''
    '''

    if not os.path.exists(opath):
        os.makedirs(opath)
    else:
        shutil.rmtree(opath)
        os.makedirs(opath)


    fcnt = 1
    file_count = len(os.listdir('../db/profile'))
    for fname in os.listdir(ipath):
        print fcnt, fname
        f = open(ipath + fname)
        data = []
        fdata = {}
        for oline in f:
            oline = oline.strip().split("\t")
            if len(oline) > 1:
                if int(oline[0]) not in fdata:
                    fdata[int(oline[0])] = []
                fdata[int(oline[0])].append(oline[1:])

        for i in range(file_count):
            if i not in fdata:
                line = [np.nan] * sample_size
                aray = np.array(line, dtype = 'f')
            else:
                line = process(fdata[i], fname, i, sample_size)
                try:
                    aray = np.array(line, dtype = 'f')
                except:
                    continue
                add = [float(line[-1])] * (sample_size - len(line))
                aray = np.append(aray, add)
            data.append(aray)
            del aray
        data = np.array(data)
        np.save(opath + fname + ".npy", data)
        del data
        f.close()
        fcnt += 1

if __name__ == '__main__':
    path = "/Volumes/External/Shobhit_work/shobhit/"
    #sample_size = max_sample_size(path + 'human/', path + 'profile/')
    entrez_unip = mapping("../db/gene_unip.txt")
    temporary_files(path + 'profile/', entrez_unip, path + 'profile_tmp/')
    #create_numpy_files("../db/profile_tmp/", "../db/profile_binary/",
    #        sample_size)
