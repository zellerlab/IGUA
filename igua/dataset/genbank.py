import gzip
import io
import pathlib
import typing

import gb_io
import Bio.Seq
import pandas as pd
import rich.progress

from .base import BaseDataset, Cluster, Protein

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

    def __init__(
        self,
        path: pathlib.Path,
        *,
        gene_feature: str = "CDS",
    ):
        """Create a new GenBank dataset.

        Args:
            input (`pathlib.Path`): The path to a GenBank file.
            gene_feature (`str`): The GenBank feature to extract gene
                sequences from. Defaults to *CDS* which is used by
                most annotation tools.

        """
        super().__init__()
        self.path = path
        self.gene_feature = gene_feature

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
                        filter(lambda f: f.kind == self.gene_feature, record.features)
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
                            rich.print(  # FIXME: use warning
                                f"[bold yellow]{'Warning':>12}[/] no "
                                f"'translation' qualifier found in "
                                f"{self.gene_feature!r} feature of "
                                f"{record.name!r}"
                            )
                            loc = feat.location
                            gene = record.sequence[loc.start : loc.end]
                            translation = self._translate_orf(gene)
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
        seq = Bio.Seq.MutableSeq(sequence)
        rest = len(seq) % 3
        if rest != 0:
            rich.print(  # FIXME: use warning
                f"[bold yellow]{'Warning':>12}[/] CDS feature location "
                "contains an incomplete codon"
            )
            seq.extend('N' * (3 - rest))
        return str(seq.translate(translation_table))