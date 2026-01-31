import unittest

import numpy
import scipy.sparse
import scipy.spatial.distance

from igua.hca import manhattan, manhattan_pair


class TestManhattan(unittest.TestCase):

    def test_unweighted(self):
        rng = numpy.random.Generator(numpy.random.MT19937(42))
        for i in range(20):
            # create a somewhat sparse matrix 
            X = rng.integers(1, 100, size=(100, 30), dtype="int32")
            X *= (rng.random(size=X.shape) < 0.2)
            S = scipy.sparse.csr_matrix(X)
            # compute expected pdist with scipy
            d_exp = scipy.spatial.distance.pdist(X, metric="cityblock")
            # compute result with custom sparse implementation
            d_out = numpy.zeros(d_exp.shape, dtype=numpy.float64)
            manhattan(S.data, S.indices, S.indptr, None, d_out)
            # ensure equality
            self.assertEqual(d_exp.tolist(), d_out.tolist())

    def test_weighted(self):
        rng = numpy.random.Generator(numpy.random.MT19937(42))
        for i in range(20):
            # create a somewhat sparse matrix 
            X = rng.integers(1, 100, size=(100, 30), dtype="int32")
            X *= (rng.random(size=X.shape) < 0.2)
            S = scipy.sparse.csr_matrix(X)
            # generate random weights
            weights = rng.integers(1, 10, X.shape[1], dtype="int64")
            # compute expected pdist with scipy
            d_exp = scipy.spatial.distance.pdist(X, metric="cityblock", w=weights)
            # compute result with custom sparse implementation
            d_out = numpy.zeros(d_exp.shape, dtype=numpy.float64)
            manhattan(S.data, S.indices, S.indptr, weights, d_out)
            # ensure equality
            self.assertEqual(d_exp.tolist(), d_out.tolist())
            # compute result as if unweighted w/ manually applied weights
            X2 = X * weights
            S2 = scipy.sparse.csr_matrix(X2, dtype="int32")
            # compute using non-weighted implementation and weighted matrix
            d_out2 = numpy.zeros(d_exp.shape, dtype=numpy.float64)
            manhattan(S2.data, S2.indices, S2.indptr, None, d_out2)
            # ensure equality
            self.assertEqual(d_out.tolist(), d_out2.tolist())