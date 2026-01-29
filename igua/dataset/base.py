import abc
import pathlib
import typing

import Bio.Seq
import pandas as pd
import rich.progress

from ..mmseqs import Database, MMSeqs
from ..sink import BaseRecordSink


class BaseDataset(abc.ABC):
    """Base class for dataset extraction."""

    @abc.abstractmethod
    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: BaseRecordSink,
    ) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        output: BaseRecordSink,
        representatives: typing.Container[str],
    ) -> pd.DataFrame:
        pass

    def translate_orf(
        self, sequence: typing.Union[str, bytes], translation_table: int = 11
    ) -> str:
        return str(Bio.Seq.Seq(sequence).translate(translation_table))
