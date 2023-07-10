from numpy import array, savetxt, random

def load_feature_data(fname):
    data = []
    file_data = open(fname)
    for line in file_data:
        line = line.strip().split()
        if 'None' not in line and 'nan' not in line:
            data.append(line[:])
    # data = array(data, dtype=float)
    data = array(data)

    return data


data = load_feature_data('irefweb_hc_no_hubs.txt') 
print data.shape
data = data[random.choice(data.shape[0],
                          544, replace=False), :]
savetxt('positive_semi_latest_hc_noNaN_544.txt', data, fmt='%s',
        delimiter='\t')