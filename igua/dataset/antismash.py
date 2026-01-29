import gzip
import io
import re
import pathlib
import typing
import zipfile

import gb_io
import pandas as pd
import rich.progress

from .base import BaseDataset, Cluster, Protein
from ..sink import BaseRecordSink

_GZIP_MAGIC = b"\x1f\x8b"


def _extract_clusters_from_record(
    record: gb_io.Record,
    source: str,
) -> typing.Iterable[Cluster]:
    for region in filter(lambda feat: feat.kind == "region", record.features):
        start = region.location.start
        end = region.location.end
        number = next(q.value for q in region.qualifiers if q.key == "region_number")
        region_id = f"{record.name}.region{number:>03}"
        yield Cluster(
            id=region_id,
            sequence=record.sequence[start-1:end].decode('ascii'),
            source=source,
        )


def _extract_proteins_from_record(
    record: gb_io.Record,
    clusters: typing.Container[str],
) -> typing.Iterable[Protein]:
    for region in filter(lambda feat: feat.kind == "region", record.features):
        number = next(q.value for q in region.qualifiers if q.key == "region_number")
        region_id = f"{record.name}.region{number:>03}"
        if region_id not in clusters:
            continue
        for cds in filter(
            lambda feat: feat.kind == "CDS" and feat.location.start >= region.location.start and feat.location.end <= region.location.end,
            record.features,
        ):
            locus_tag = next(q.value for q in cds.qualifiers if q.key == "locus_tag")
            translation = next(q.value for q in cds.qualifiers if q.key == "translation")
            yield Protein(
                f"{region_id}_{locus_tag}",
                translation,
                cluster_id=region_id
            )


class AntiSMASHGenBankDataset(BaseDataset):
    """A dataset composed of antiSMASH regions in a GenBank file.

    antiSMASH reports regions in a GenBank file but not necessarily
    with one region per record. This class supports extracting regions
    independently and recovering the right cluster ID per region.

    """

    def __init__(
        self,
        path: pathlib.Path,
        *,
        mode: str = "region",
    ):
        super().__init__()

        if mode not in {"region", "protocluster", "core"}:
            raise ValueError(f"invalid mode: {mode!r}")
        if mode in {"protocluster", "core"}:
            raise NotImplementedError(f"not implemented mode: {mode!r}")

        self.mode = mode
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
                yield from _extract_clusters_from_record(record, source=str(self.path))
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
                yield from _extract_proteins_from_record(record, clusters)
        progress.remove_task(task)


class AntiSMASHZipDataset(BaseDataset):
    """A dataset composed of antiSMASH results in a Zip file.

    antiSMASH can be configured to report all results inside a Zip file.
    This class supports reading the antiSMASH-predicted regions from a
    Zip archive, including handling of region IDs, without requiring the
    archive to be decompressed.

    """

    _REGION_RX = re.compile(r"region(\d){3,}.gbk$")

    def __init__(
        self,
        path: pathlib.Path,
        *,
        mode: str = "region",
    ):
        super().__init__()

        if mode not in {"region", "protocluster", "core"}:
            raise ValueError(f"invalid mode: {mode!r}")
        if mode in {"protocluster", "core"}:
            raise NotImplementedError(f"not implemented mode: {mode!r}")

        self.mode = mode
        self.path = path

    def extract_clusters(
        self,
        progress: rich.progress.Progress,
    ) -> typing.Iterable[Cluster]:
        with zipfile.ZipFile(self.path, mode="r") as archive:
            for file in archive.filelist:
                if self._REGION_RX.search(file.filename) is None:
                    continue
                with archive.open(file, mode="r") as src:
                    records = gb_io.iter(src)
                    record = next(records)
                    if next(records, None) is not None:
                        raise ValueError(f"more than one record found in {file.filename}")
                    yield from _extract_clusters_from_record(record, source=str(self.path))

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        clusters: typing.Container[str],
    ) -> typing.Iterable[Protein]:
        with zipfile.ZipFile(self.path, mode="r") as archive:
            files = { file.filename for file in archive.filelist }
            for file in archive.filelist:
                if self._REGION_RX.search(file.filename) is None:
                    continue
                with archive.open(file, mode="r") as src:
                    records = gb_io.iter(src)
                    record = next(records)
                    if next(records, None) is not None:
                        raise ValueError(f"more than one record found in {file.filename}")
                    yield from _extract_proteins_from_record(record, clusters)

