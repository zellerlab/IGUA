import abc
import typing
from typing import Literal

import numpy
import scipy.sparse
from scipy.cluster.hierarchy import fcluster

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
        X: scipy.sparse.spmatrix,
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
