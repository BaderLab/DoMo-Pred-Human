import networkx as nx
from numpy import loadtxt, log, column_stack, row_stack, array, savetxt, NaN, isnan
import cPickle as pickle
from scipy.misc import comb
from sklearn import linear_model
from sklearn import preprocessing

'''
graph_tool.clustering.global_clustering
graph_tool.stats.vertex_average

'''


def load_network_data(path, out_path):
    '''
    '''
    edges = loadtxt(path, dtype=str)
    graph = nx.Graph()
    graph.add_edges_from(edges)

    return graph


def graph_cliques(graph):
    '''
    '''
    cliques = list(nx.find_cliques(graph))

    return cliques


def jaccard(n1, n2):
    '''
    '''
    return len(n1 & n2) / float(len(n1 | n2))


def meetmin(n1, n2):
    '''
    '''
    return len(n1 & n2) / float(min(len(n1), len(n2)))


def geometric(n1, n2):
    '''
    '''
    return len(n1 & n2) ** 2 / float(len(n1) * len(n2))


def hyper_geometric(n1, n2, total):
    '''
    '''
    sat = len(n1 & n2)
    end = min(len(n1), len(n2))
    res = 0
    for i in range(sat, end + 1):
        num = comb(len(n1), i) * comb((total - len(n1)), (len(n2) - 1))
        den = comb(total, len(n2))
        res = num / float(den)
    print res
    return -log(res)


def graph_properties(graph):
    '''
    Return a list of subgraph/cluster properties.
    '''

    density = nx.density(graph)
    degrees = nx.degree(graph)
    avg_deg = sum(degrees.values()) / float(nx.number_of_nodes(graph))
    edg_con = nx.edge_connectivity(graph)
    cls_val = nx.transitivity(graph)

    return [density, avg_deg, edg_con, cls_val]


def save_graph(graph):
    '''
    '''
    nds = nx.nodes(graph)
    data = {}
    subs = {}
    prop = {}
    count = 1
    for node in nds:
        data[node] = set(graph.neighbors(node))
        data[node].add(node)
        subs[node] = graph.subgraph(data[node])
        prop[node] = graph_properties(subs[node])
        print count
        count += 1
    with open('nt_neighbours.pck', 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    # with open('nt_subgraphs.pck', 'wb') as f:
    #    pickle.dump(subs, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open('nt_properties.pck', 'wb') as f:
        pickle.dump(prop, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_graph(graph_file, neighbours='', subgraphs='', properties=''):
    '''
    '''
    print 'Loading graph information...'
    # subgs = None
    graph = load_network_data(graph_file, '')
    with open('nt_neighbours_copy.pck', 'rb') as f:
        neigh = pickle.load(f)
    # with open('nt_subgraphs.pck', 'rb') as f:
    #     subgs = pickle.load(f)
    with open('nt_properties_copy.pck', 'rb') as f:
        props = pickle.load(f)
    print 'Finished loading...'
    return graph, neigh, props


def explore_neighbourhood(a, b, graph, neigh, props):
    '''
    '''
    if a not in neigh or b not in neigh:
        return NaN
    n1 = neigh[a]
    n2 = neigh[b]
    sub3 = graph.subgraph(n1 | n2)
    prop = props[a][:]
    prop.extend(props[b][:])
    prop.extend(graph_properties(sub3)[:])
    prop.append(meetmin(n1, n2))

    return array(prop)


def classifier(x, y):
    '''
    '''
    clf = linear_model.LogisticRegression(solver='lbfgs', C=1.0)
    clf.fit(x, y)

    return clf


def scaled_matrix(data):
    scaled_mat = array([]).reshape(data.shape[0], 0)
    for i in range(data.shape[1]):
        minmax = preprocessing.MinMaxScaler()
        scaled_i = minmax.fit_transform(data[:, i].reshape(-1, 1))
        scaled_mat = column_stack((scaled_mat, scaled_i))

    return scaled_mat


def train_classifier(positive, negative):
    '''
    '''
    train_pos = loadtxt(positive, dtype='str')[:, 2:].astype(float)
    train_neg = loadtxt(negative, dtype='str')[:, 2:].astype(float)
    # train_pos = scaled_matrix(train_pos)
    # train_neg = scaled_matrix(train_neg)
    data1 = column_stack((train_pos, [[0]] * len(train_pos)))
    data2 = column_stack((train_neg, [[1]] * len(train_neg)))
    training = row_stack((data1, data2))
    clf = classifier(training[:, :-1], training[:, -1])

    return clf


def load_topology():
    '''
    '''
    graph, neigh, props = load_graph('irefweb_all.txt')
    clf = train_classifier('pos.txt', 'neg.txt')

    return graph, neigh, props, clf


def topology(a, b, graph, neigh, props, clf):
    '''
    '''
    case = explore_neighbourhood(a, b, graph, neigh, props)
    if isnan(case).all():
        return NaN
    return clf.predict_proba(case.reshape(1, -1))[0][0]


if __name__ == '__main__':
    '''
    '''
    graph, neigh, props, clf = load_topology()
    data = loadtxt('unlbl_pos_semi.txt', dtype='str')
    res = []
    count = 1
    for i in data:
        print count
        if count % 500 == 0:
            savetxt('unlbl_%s_pos_semi.txt' % count, res, fmt='%s')
        count += 1
        val = topology(i[0], i[1], graph, neigh, props, clf)
        res.append(val)
    data = column_stack((data, res))
    savetxt('unlbl_new_pos_semi.txt', data, fmt='%s')

