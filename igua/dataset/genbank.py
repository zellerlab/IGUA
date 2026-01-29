import gzip
import io
import pathlib
import typing

import gb_io
import pandas as pd
import rich.progress

from .base import BaseDataset
from ..sink import BaseRecordSink

_GZIP_MAGIC = b"\x1f\x8b"


class GenBankDataset(BaseDataset):
    """GenBank dataset class."""

    def __init__(self, path: pathlib.Path):
        """Initialize GenBank dataset.

        Args:
            inputs: List of GenBank file paths.
        """
        super().__init__()
        self.path = path

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: BaseRecordSink,
    ) -> pd.DataFrame:
        """Extracts nucleotide sequences from GenBank files.
        
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted sequences.
        
        Returns:
            `pandas.DataFrame`: DataFrame containing the extracted sequences 
            and metadata with columns ``cluster_id``, ``cluster_length`` and
            ``filename``.

        """
        data = []
        done = set()
        n_duplicate = 0

        task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
        with io.BufferedReader(progress.open(self.path, "rb", task_id=task)) as reader:  # type: ignore
            if reader.peek().startswith(_GZIP_MAGIC):
                reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
            for record in gb_io.iter(reader):
                if record.name in done:
                    n_duplicate += 1
                else:
                    output.add_record(record.name, record.sequence.decode('ascii'))
                    data.append((record.name, len(record.sequence), self.path))
                    done.add(record.name)
        progress.remove_task(task)
        
        if n_duplicate > 0:
            progress.console.print(
                f"[bold yellow]{'Skipped':>12}[/] {n_duplicate} "
                f"clusters with duplicate identifiers"
            )
        return pd.DataFrame(
            data=data, 
            columns=["cluster_id", "cluster_length", "filename"]
        ).set_index("cluster_id")

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        output: BaseRecordSink,
        representatives: typing.Container[str],
    ) -> pd.DataFrame:
        """Extracts protein sequences from GenBank files.
        
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted protein sequences.
            representatives (typing.Container[str]): Set of representative cluster IDs.
        
        Returns:
            `pandas.DataFrame`: DataFrame containing the extracted proteins 
            and metadata with columns ``cluster_id``, ``protein_id``,
            ``protein_length`` and ``filename``.

        """
        data = []
        done = set()
        
        task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
        with io.BufferedReader(progress.open(self.path, "rb", task_id=task)) as reader:  # type: ignore
            if reader.peek()[:2] == b"\x1f\x8b":
                reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
            for record in gb_io.iter(reader):
                if record.name in representatives:
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
                                f"[bold yellow]{'Warning':>12}[/] no 'translation' qualifier found in CDS feature of {record.name!r}"
                            )
                            translation = self.translate_orf(
                                record.sequence[
                                    feat.location.start : feat.location.end
                                ]
                            )
                        else:
                            translation = qualifier.value.rstrip("*")
                        protein_id = "{}_{}".format(record.name, i)
                        if protein_id not in done:
                            output.add_record(protein_id, translation)
                            data.append((record.name, protein_id, len(translation), self.path))
                            done.add(protein_id)
        progress.remove_task(task)

        df = pd.DataFrame(
            data=data, 
            columns=["cluster_id", "protein_id", "protein_length", "filename"],
        )
        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(df):,} proteins from "
            f"{len(representatives):,} nucleotide representative"
        )
        return df
