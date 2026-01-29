import gzip
import io
import pathlib
import typing

import gb_io
import Bio.Seq
import pandas as pd
import rich.progress

from .base import BaseDataset, Cluster, Protein
from ..sink import BaseRecordSink

_GZIP_MAGIC = b"\x1f\x8b"


class GenBankDataset(BaseDataset):
    """A dataset composed of gene clusters in a GenBank file.

    GenBank files are commonly used to distribute a genomic sequence
    with associated metadata inside a single file.

    Note:
        This method treats each GenBank record as an independent gene
        cluster and extracts the full record sequence and all the 
        annotated genes. For GenBank files obtained with antiSMASH,
        please use the `AntiSMASHGenBankDataset` class to enable
        additional processing of the regions.

    """

    def __init__(self, path: pathlib.Path):
        """Create a new GenBank dataset.

        Args:
            input (`pathlib.Path`): The path to a GenBank file.

        """
        super().__init__()
        self.path = path

    def extract_clusters(
        self,
        progress: rich.progress.Progress,
    ) -> typing.Iterable[Cluster]:
        task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
        with io.BufferedReader(progress.open(self.path, "rb", task_id=task)) as reader:  # type: ignore
            if reader.peek().startswith(_GZIP_MAGIC):
                reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
            for record in gb_io.iter(reader):
                yield Cluster(
                    id=record.name,
                    sequence=record.sequence.decode('ascii'),
                    source=str(self.path),
                )
        progress.remove_task(task)

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        clusters: typing.Container[str],
    ) -> typing.Iterable[Protein]:
        task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
        with io.BufferedReader(progress.open(self.path, "rb", task_id=task)) as reader:  # type: ignore
            if reader.peek()[:2] == b"\x1f\x8b":
                reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
            for record in gb_io.iter(reader):
                if record.name in clusters:
                    for i, feat in enumerate(
                        filter(lambda f: f.kind == "CDS", record.features)
                    ):
                        qualifier = next(
                            (
                                qualifier
                                for qualifier in feat.qualifiers
                                if qualifier.key == "translation"
                            ),
                            None,
                        )
                        if qualifier is None:
                            rich.print(
                                f"[bold yellow]{'Warning':>12}[/] no "
                                "'translation' qualifier found in CDS "
                                f"feature of {record.name!r}"
                            )
                            translation = self._translate_orf(
                                record.sequence[
                                    feat.location.start : feat.location.end
                                ]
                            )
                        else:
                            translation = qualifier.value.rstrip("*")
                        yield Protein(
                            id="{}_{}".format(record.name, i),
                            sequence=translation,
                            cluster_id=record.name,
                        )
        progress.remove_task(task)

    def _translate_orf(
        self, sequence: typing.Union[str, bytes], translation_table: int = 11
    ) -> str:
        return str(Bio.Seq.Seq(sequence).translate(translation_table))