import numpy as np
import sys
import cPickle
import dill


class gnb(object):
    '''
    '''

    def __init__(self, num_features, num_classes, prior=None):
        '''
        prior: number of classes
        means: num features x num classes
        varns: num features x num classes
        '''
        self.__feats = num_features
        self.__class = num_classes
        self.__prior = np.ones(num_classes)
        self.__means = np.zeros((num_features, num_classes))
        self.__varns = np.zeros((num_features, num_classes))

    def __mean(self, data, delta):
        '''
        '''
        # return np.nanmean(col)

        return sum(data * delta) / float(sum(delta))

    def __varn(self, data, delta, mean):
        '''
        '''
        # return np.nanvar(col, ddof=1)
        # res = sum(((data - mean) ** 2) * delta) / max(float(sum(delta) - 1), 1)
        v1 = float(sum(delta))
        # n = len(delta)
        v2 = float(sum(delta ** 2))
        res = sum(((data - mean) ** 2) * delta) / (v1 - (v2 / v1))
        return res

    def train(self, data, delta):
        '''
        '''
        assert delta.shape[1] == self.__class, 'Mismatch in num of classes'
        assert data.shape[1] == self.__feats, 'Mismatch in num of features'
        assert data.shape[0] == delta.shape[0], 'Data - delta mismatch'

        self.__prior = np.sum(delta, axis=0).astype(float)
        self.__prior /= sum(self.__prior)
        print self.__prior

        for j in range(self.__class):
            for i in range(self.__feats):
                feat = np.column_stack((data[:, i], delta[:, j]))
                feat = feat[~np.isnan(feat).any(axis=1)]
                self.__means[i][j] = self.__mean(feat[:, 0], feat[:, 1])
                self.__varns[i][j] = self.__varn(feat[:, 0], feat[:, 1],
                                                 self.__means[i][j])

        # print self.__means
        # print self.__varns

    def __gauss(self, x, mean, var):
        '''
        '''
        return (1 / (np.sqrt(2 * np.pi * var)))\
            * np.exp(-(x - mean) ** 2 / (2 * var))

    def posterior(self, case, semi=False, norm=True):
        '''
        '''
        prob = np.zeros(self.__class)
        for j in range(self.__class):
            temp = np.ones(self.__feats)
            for i in range(self.__feats):
                if not np.isnan(case[i]):
                    temp[i] = self.__gauss(case[i], self.__means[i][j],
                                           self.__varns[i][j])
            prob[j] = np.prod(temp) * self.__prior[j]
        if norm:
            prob = prob / sum(prob)
        if semi:
            return prob
        return np.argmax(prob), np.max(prob)

    def reset(self):
        '''
        '''
        self.__prior = np.ones(self.__class)
        self.__means = np.zeros((self.__feats, self.__class))
        self.__varns = np.zeros((self.__feats, self.__class))

    def print_params(self):
        '''
        '''
        print 'Prior:'
        print self.__prior
        print 'Means:'
        print self.__means
        print 'Varns:'
        print self.__varns


class seminb(gnb):
    '''
    '''

    def __init__(self, num_features, num_classes, prior=None):
        '''
        '''
        gnb.__init__(self, num_features, num_classes, prior)

    def posterior_all(self, data, norm=True):
        '''
        '''
        delta = np.zeros((data.shape[0], 2))
        for i in range(data.shape[0]):
            delta[i, :] = self.posterior(data[i, :], True, norm)
        return delta

    def posterior_labelled(self, ld, delta):
        '''
        '''
        prob = self.posterior_all(ld, norm=False)
        assert prob.shape == delta.shape, 'Some Issue 1'
        return np.ma.masked_equal(prob * delta, 0)

    def posterior_unlabelled(self, uld):
        '''
        '''
        prob = self.posterior_all(uld, norm=False)
        return np.ma.masked_equal(np.sum(prob, axis=1), 0)

    def logsum(self, loga, k=-np.inf):
        """
        Compute a sum of logs without underflow.
        \log \sum_c e^{b_c} = log [(\sum_c e^{b_c}) e^{-B}e^B]
                            = log [(\sum_c e^{b_c-B}) e^B]
                            = [log(\sum_c e^{b_c-B}) + B
        where B = max_c b_c
        """
        B = np.max(loga)
        logaB = aB = loga - B
        sup = logaB > k
        inf = np.logical_not(sup)
        aB[sup] = np.exp(logaB[sup])
        aB[inf] = 0.0
        return (np.log(np.sum(aB)) + B)

    def likelihood(self, ld, delta, uld):
        '''
        '''
        lik = 0.0
        # labelled
        prob = self.posterior_labelled(ld, delta)
        lik_l = np.ma.sum(np.ma.log(prob))

        # unlabelled
        prob = self.posterior_unlabelled(uld)
        lik_u = np.ma.sum(np.ma.log(prob))

        lik = lik_l + (1 * lik_u)
        return lik

    def train_semi(self, data, delta, unlabeled=None, maxiter=100, eps=0.01,
                   debug=True):
        '''
        '''
        self.train(data, delta)

        if unlabeled is None:
            print "No unlabeled data"
            return

        like = self.likelihood(data, delta, unlabeled)

        for iteration in range(1, maxiter + 1):
            # E-step
            delta_u = self.posterior_all(unlabeled)

            # M-step
            data_combn = np.row_stack((data, unlabeled))
            delta_u = delta_u * 1
            delta_comb = np.row_stack((delta, delta_u))
            self.reset()
            self.train(data_combn, delta_comb)

            new_like = self.likelihood(data, delta, unlabeled)
            diff = new_like - like
            like = new_like

            if diff < eps:
                print "Stopping EM at iteration ", iteration
                break

            if debug:
                print "Iteration: ", iteration, " Difference: ", diff


def save_model(obj, filename):
    '''
    '''
    with open(filename, 'wb') as output:
        # cPickle.dump(obj, output, cPickle.HIGHEST_PROTOCOL)
        dill.dump(obj, output)


def load_model(filename):
    '''
    '''
    with open(filename, 'rb') as input:
        # model = cPickle.load(input)
        model = dill.load(input)
    return model


if __name__ == '__main__':

    cls = seminb(3, 2)
    data = np.array([[1, 1, 1], [0.9, 0.9, 0.9], [0.1, 0, 1], [0, 1, 0]])
    delta = np.array([[1, 0], [1, 0], [0, 1], [0, 1]])
    unlbl = np.array([[1, 1, 0.9], [0.5, 0.5, 0.5], [0.5, 0.1, 0.1]])
    cls.train_semi(data, delta, unlbl)
    print cls.posterior([1, 1, 1])
    save_model(cls, 'test.pck')
    del cls
    cls1 = load_model('test.pck')
    print cls1.posterior([1, 1, 1])
