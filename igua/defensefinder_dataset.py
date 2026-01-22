import pathlib
import typing

from .fasta_gff_dataset import ClusterDataAdapter
from .fasta_gff_dataset import FastaGFFDataset


class DefenseFinderDataset(FastaGFFDataset):
    """DefenseFinder-specific dataset class.

    Inherits all functionality from FastaGFFDataset but provides
    explicit type identity for format-specific protein ID parsing.
    Uses double-underscore delimiter (cluster_id__protein_id) for protein IDs.
    """

    def __init__(
        self,
        inputs: typing.List[pathlib.Path],
        adapter: ClusterDataAdapter,
    ) -> None:
        """Initialize DefenseFinder dataset.

        Args:
            inputs: List of input file paths (typically metadata TSV).
            adapter: Adapter for handling DefenseFinder-specific logic.
        """
        super().__init__(inputs, adapter)
