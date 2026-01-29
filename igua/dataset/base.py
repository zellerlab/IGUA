import abc
import dataclasses
import pathlib
import typing

import Bio.Seq
import pandas as pd
import rich.progress

from ..mmseqs import Database, MMSeqs


@dataclasses.dataclass
class Cluster:
    """A gene cluster.

    Attributes:
        id (`str`): The identifier of the gene cluster. Must be unique 
            across all datasets, as the identifier will be used to refer
            to each cluster in the resulting gene cluster families.
        sequence (`str`): The genomic sequence of the dataset.
        source (`str` or `None`): The source of the gene cluster, 
            usually the path to the file where this cluster record
            was extracted from.

    """
    id: str
    sequence: str
    source: typing.Optional[str] = dataclasses.field(default=None)

    @property
    def length(self):
        """`int`: The length of the cluster sequence.
        """
        return len(self.sequence)


@dataclasses.dataclass
class Protein:
    """A protein inside a gene cluster.

    Attributes:
        id (`str`): The identifier of the protein. Must be unique 
            across all datasets, as the identifier will be used to refer
            to each protein while clustering.
        sequence (`str`): The amino-acid sequence of the protein.
        cluster_id (`str`): The identifier of the cluster this protein
            belongs to.

    """
    id: str
    sequence: str
    cluster_id: str

    @property
    def length(self):
        """`int`: The length of the protein sequence.
        """
        return len(self.sequence)


class BaseDataset(abc.ABC):
    """An abstract dataset to provide clusters to a `ClusteringPipeline`.

    IGUA needs the genomic sequences of the gene clusters to process, and
    the individual proteins of each gene cluster. To provide them to 
    a clustering pipeline, the `BaseDataset` interface provides two 
    functions: `BaseDataset.extract_clusters` and 
    `BaseDataset.extract_proteins` which yield `Cluster` and `Protein`
    objects, respectively.

    """

    @abc.abstractmethod
    def extract_clusters(
        self, 
        progress: rich.progress.Progress,
    ) -> typing.Iterable[Cluster]:
        """Extract the clusters from the dataset.

        Arguments:
            progress (`rich.progress.Progress`): A `Progress` instance
                that can be used for tracking progress.

        Yields:
            `Cluster`: A cluster object for each gene cluster to be 
            processed in the dataset.

        """

    @abc.abstractmethod
    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        cluster_ids: typing.Container[str],
    ) -> typing.Iterable[Protein]:
        """Extracts protein sequences from GenBank files.

        Arguments:
            progress (`rich.progress.Progress`): A `Progress` instance
                that can be used for tracking progress.
            clusters (`collections.abc.Container` of `str`): A container
                storing which clusters to extract proteins from.

        Yields:
            `Protein`: A protein object for each protein of the gene 
            clusters to be processed in the dataset.

        """

