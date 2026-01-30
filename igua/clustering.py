import abc
import typing
from typing import Literal

import numpy
import scipy.sparse
from scipy.cluster.hierarchy import fcluster, DisjointSet

from .hca import linkage, manhattan


class BaseClustering(abc.ABC):

    @abc.abstractmethod
    def cluster(self, X: scipy.sparse.spmatrix) -> numpy.ndarray:
        """Cluster the given observations.
        """


class HierarchicalClustering(BaseClustering):
    """A clustering method implemeting hierarchical clustering.
    """

    def __init__(
        self,
        *,
        method: Literal["average", "single", "complete", "weighted", "centroid", "median", "ward"],
        distance: float = 0.8,
        precision: Literal["half", "single", "double"] = "double",
        jobs: int = 1,
    ):
        self.method = method
        self.distance = distance
        self.precision = precision
        self.jobs = jobs

    def _compute_distances(
        self,
        X: scipy.sparse.csr_matrix,
    ) -> numpy.ndarray:
        # check matrix format
        if not isinstance(X, scipy.sparse.csr_matrix):
            raise TypeError(f"expected csr_matrix, got {type(X).__name__}")

        r = X.shape[0]
        # compute the number of amino acids per cluster
        clusters_aa = numpy.zeros(r, dtype=numpy.int32)
        clusters_aa[:] = X.sum(axis=1).A1
        # make sure the sparse matrix has sorted indices (necessary for
        # the distance algorithm to work efficiently)
        if not X.has_sorted_indices:
            X.sort_indices()
        # compute manhattan distance on sparse matrix
        distance_vector = numpy.zeros(r * (r - 1) // 2, dtype=self.precision)
        manhattan(
            X.data,
            X.indices,
            X.indptr,
            distance_vector,
            threads=self.jobs,
        )
        # ponderate by sum of amino-acid distance
        n = 0
        for i in range(r - 1):
            l = r - (i + 1)
            distance_vector[n : n + l] /= (clusters_aa[i + 1 :] + clusters_aa[i]).clip(
                min=1
            )
            n += l
        # enforce distances to be in [0, 1] (slight possibility of >1 due
        # to handling of alignment weights but we can round off)
        return numpy.clip(distance_vector, 0.0, 1.0, out=distance_vector)

    def cluster(
        self,
        X: scipy.sparse.spmatrix,
    ) -> numpy.ndarray:
        pdist = self._compute_distances(X)
        Z = linkage(pdist, method=self.method)
        return fcluster(Z, criterion="distance", t=self.distance)


class LinearClustering(BaseClustering):
    """A clustering method similar to MMseqs2 linear clustering.
    """

    def __init__(
        self,
        *,
        distance: float = 0.8,
    ):
        self.distance = distance

    def cluster(
        self,
        X: scipy.sparse.csr_matrix
    ):
        # get a csr_matrix form the input to support fast selection
        # of columns
        X_csc = X.tocsc()

        # compute the number of amino acids per gene cluster
        # FIXME: may not be reliable if we can disable weighting, replace
        #        by an additional argument most likely
        r = X.shape[0]
        clusters_aa = numpy.zeros(r, dtype=numpy.int32)
        clusters_aa[:] = X.sum(axis=1).A1

        # use scipy DisjointSet / UnionFind to record clustering
        ds = DisjointSet(range(r))

        # extract columns with more than one row so we only iterate on those
        n_obs = (X_csc != 0).sum(axis=0)
        cols = numpy.where(n_obs > 1)[1]

        # iterate on subsets formed by protein membership and cluster together
        # all observations with relative weighted Manhattan under the distance
        # threshold, taking the largest observation (in aa.) as the centroid
        for y in cols:
            # extract which gene clusters contain protein 'y'
            # (indexing is faster with the compressed-sparse-column matrix)
            mask = X_csc[:, y]
            if mask.nnz <= 1:
                continue
            # extract the row subset indices directly from the CSC (r x 1) array
            indices = mask.indices
            # find the largest gene cluster of the row subset
            centroid = indices[numpy.argmax(clusters_aa[indices])]
            # compare every other gene clusters to the reference
            for query in indices:
                d = numpy.abs((X[query] - X[centroid]).toarray()).sum()
                d /= (clusters_aa[query] + clusters_aa[centroid]).clip(min=1.0)
                if d <= self.distance:
                    ds.merge(query, centroid)

        # extract the cluster indices from the disjoint set
        # NOTE: implementation not super efficient at the moment as
        #       `ds.subsets` returns a list of set that we need to convert
        #       again, maybe prohibitive for very large datasets
        flat = numpy.zeros(r, dtype=numpy.int32)
        for i, subset in enumerate(ds.subsets(), 1):
            flat[list(subset)] = i

        return flat
