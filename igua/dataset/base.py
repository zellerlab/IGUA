import abc
import pathlib
import typing

import Bio.Seq
import pandas as pd
import rich.progress

from ..mmseqs import Database, MMSeqs


class BaseDataset(abc.ABC):
    """Base class for dataset extraction."""

    def __init__(self, inputs: typing.List[pathlib.Path]):
        """Initialize dataset with input file paths.

        Args:
            inputs: List of input file paths to process.
        """
        self.inputs = inputs

    @abc.abstractmethod
    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
    ) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> pd.DataFrame:
        pass

    def write_fasta(self, file: typing.TextIO, name: str, sequence: str) -> None:
        file.write(">{}\n".format(name))
        file.write(sequence)
        file.write("\n")
        return None

    def translate_orf(
        self, sequence: typing.Union[str, bytes], translation_table: int = 11
    ) -> str:
        return str(Bio.Seq.Seq(sequence).translate(translation_table))

    def create_sequence_database(
        self,
        mmseqs: MMSeqs,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        output_db_path: pathlib.Path,
    ) -> Database:
        """Default implementation creates a temporary file then a database"""
        tmp_fasta = output_db_path.with_suffix(".fna")
        self.extract_sequences(progress, inputs, tmp_fasta)
        return Database.create(mmseqs, tmp_fasta)

    def create_protein_database(
        self,
        mmseqs: MMSeqs,
        progress: rich.progress.Progress,
        inputs: typing.List[pathlib.Path],
        representatives: typing.Container[str],
        output_db_path: pathlib.Path,
    ) -> typing.Tuple[Database, typing.Dict[str, int]]:
        """Default implementation creates a temporary file then a database"""
        tmp_fasta = output_db_path.with_suffix(".faa")
        protein_sizes = self.extract_proteins(
            progress, inputs, tmp_fasta, representatives
        )
        return Database.create(mmseqs, tmp_fasta), protein_sizes
