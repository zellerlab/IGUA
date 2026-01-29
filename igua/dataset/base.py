import abc
import dataclasses
import pathlib
import typing

import Bio.Seq
import pandas as pd
import rich.progress

from ..mmseqs import Database, MMSeqs
from ..sink import BaseRecordSink


@dataclasses.dataclass
class Cluster:
    id: str
    sequence: str
    source: typing.Optional[str] = dataclasses.field(default=None)

    @property
    def length(self):
        return len(self.sequence)


@dataclasses.dataclass
class Protein:
    id: str
    sequence: str
    cluster_id: str

    @property
    def length(self):
        return len(self.sequence)


class BaseDataset(abc.ABC):
    """Base class for dataset extraction."""

    @abc.abstractmethod
    def extract_clusters(
        self, 
        progress: rich.progress.Progress,
    ) -> typing.Iterable[Cluster]:
        pass

    @abc.abstractmethod
    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        cluster_ids: typing.Container[str],
    ) -> typing.Iterable[Protein]:
        pass

    def translate_orf(
        self, sequence: typing.Union[str, bytes], translation_table: int = 11
    ) -> str:
        return str(Bio.Seq.Seq(sequence).translate(translation_table))
