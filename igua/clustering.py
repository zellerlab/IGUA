import abc
import typing
from typing import Literal

import numpy
import scipy.sparse
from scipy.cluster.hierarchy import fcluster, DisjointSet

from .hca import linkage, manhattan, manhattan_pair


class ClusteringStrategy(abc.ABC):
    """An abstract clustering strategy to cluster compositional data.
    """

    @abc.abstractmethod
    def cluster(
        self,
        X: scipy.sparse.csr_matrix,
        weights: typing.Optional[numpy.ndarray] = None,
    ) -> numpy.ndarray:
        r"""Cluster the given observations.

        Arguments:
            X (`scipy.sparse.csr_matrix`): A matrix of shape :math:`m \times n`
                with compositional data.
            weights (`numpy.ndarray` or `None`): The weights (of shape
                :math:`n`) to use for computing distances. If `None`,
                use uniform weights.

        Returns:
            `numpy.ndarray`: A flat array of shape :math:`m` which assigns
            an arbitrary cluster number to each observation of ``X``
            (see `scipy.cluster.hierarchy.fcluster` documentation).

        """


class HierarchicalClustering(ClusteringStrategy):
    """A clustering strategy implementing hierarchical clustering.
    """

    def __init__(
        self,
        *,
        method: Literal["average", "single", "complete", "weighted", "centroid", "median", "ward"],
        distance: float = 0.8,
        precision: Literal["half", "single", "double"] = "double",
        jobs: int = 1,
    ):
        """Create a new hierarchical clustering strategy.

        Arguments:
            method (`str`): The name of the linkage method to use: *average*,
                *single*, *complete*, *weighted*, *centroid*, *median* or
                *ward*.
            distance (`float`): The distance cutoff for the flat clusters
                created with `scipy.cluster.hierarchy.fcluster`.
            precision (`str`): The floating-point precision to use: *half*,
                *single* or *double*. Note that changing precision (in
                particular with half-precision) may lead to differences in
                the generated clustering.
            jobs (`int`): The number of parallel threads to use to perform
                the pairwise distance computation.

        """
        self.method = method
        self.distance = distance
        self.precision = precision
        self.jobs = jobs

    def _compute_distances(
        self,
        X: scipy.sparse.csr_matrix,
        weights: typing.Optional[numpy.ndarray] = None,
    ) -> numpy.ndarray:
        # check matrix format
        if not isinstance(X, scipy.sparse.csr_matrix):
            raise TypeError(f"expected csr_matrix, got {type(X).__name__}")

        # make sure the sparse matrix has sorted indices (necessary for
        # the distance algorithm to work efficiently)
        if not X.has_sorted_indices:
            X.sort_indices()

        # compute the total length of each observation
        # (used to rescale the absolute distances to a relative range)
        r = X.shape[0]
        if weights is None:
            total = X.sum(axis=1, dtype=numpy.int32).A1
        else:
            total = X @ weights

        # compute manhattan distance on sparse matrix
        distance_vector = numpy.zeros(r * (r - 1) // 2, dtype=self.precision)
        manhattan(
            X.data,
            X.indices,
            X.indptr,
            weights,
            distance_vector,
            threads=self.jobs,
        )
        # ponderate by sum of amino-acid distance
        n = 0
        for i in range(r - 1):
            l = r - (i + 1)
            maxdist = numpy.clip(total[i + 1 :] + total[i], min=1)
            distance_vector[n : n + l] /= maxdist
            n += l
        # enforce distances to be in [0, 1] (slight possibility of >1 due
        # to handling of alignment weights but we can round off)
        return numpy.clip(distance_vector, 0.0, 1.0, out=distance_vector)

    def cluster(
        self,
        X: scipy.sparse.spmatrix,
        weights: typing.Optional[numpy.ndarray] = None,
    ) -> numpy.ndarray:
        if weights is not None and X.shape[1] != weights.shape[0]:
            raise ValueError("inconsistent shapes between X and weights")
        pdist = self._compute_distances(X, weights)
        Z = linkage(pdist, method=self.method)
        return fcluster(Z, criterion="distance", t=self.distance)


class LinearClustering(ClusteringStrategy):
    """A clustering strategy similar to MMseqs2 linear clustering.
    """

    def __init__(
        self,
        *,
        distance: float = 0.8,
    ):
        """Create a new linear clustering strategy.

        Arguments:
            distance (`float`): The distance cutoff to use for clustering
                observations together.

        """
        self.distance = distance

    def cluster(
        self,
        X: scipy.sparse.csr_matrix,
        weights: typing.Optional[numpy.ndarray] = None,
    ) -> numpy.ndarray:
        # check matrix format
        if not isinstance(X, scipy.sparse.csr_matrix):
            raise TypeError(f"expected csr_matrix, got {type(X).__name__}")
        if weights is not None and X.shape[1] != weights.shape[0]:
            raise ValueError("inconsistent shapes between X and weights")

        # make sure the sparse matrix has sorted indices (necessary for
        # the distance algorithm to work without returning bogus values)
        if not X.has_sorted_indices:
            X.sort_indices()

        # also compute a csc_matrix form the input to support fast
        # selection of columns in the loop
        X_csc = X.tocsc()

        # compute the total length of each observation
        # (used to rescale the absolute distances to a relative [0, 1] range)
        r = X.shape[0]
        if weights is None:
            total = X.sum(axis=1, dtype=numpy.int32).A1
        else:
            total = X @ weights

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
            # find t he largest gene cluster of the row subset
            centroid = indices[numpy.argmax(total[indices])]
            # compare every other gene clusters to the reference
            # TODO: implement a Rust function to for manhattan distance
            #       a la cdist so we can do distance of centroid to every
            #       other query efficiently in a single call?
            for query in indices:
                d = manhattan_pair(X.data, X.indices, X.indptr, weights, query, centroid)
                d /= numpy.clip(total[query] + total[centroid], min=1.0)
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
