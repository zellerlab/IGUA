import gzip
import io
import pathlib
import typing

import gb_io
import pandas as pd
import rich.progress

from .base_dataset import BaseDataset


_GZIP_MAGIC = b"\x1f\x8b"


class GenBankDataset(BaseDataset):
    """GenBank dataset class."""

    def __init__(self, inputs: typing.List[pathlib.Path]):
        """Initialize GenBank dataset.

        Args:
            inputs: List of GenBank file paths.
        """
        super().__init__(inputs)

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
    ) -> pd.DataFrame:
        """Extracts nucleotide sequences from GenBank files.
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted sequences.
        Returns:
            pd.DataFrame: DataFrame containing the extracted sequences.
        """
        data = []
        done = set()
        n_duplicate = 0
        with open(output, "w") as dst:
            task1 = progress.add_task(f"[bold blue]{'Working':>9}[/]")
            for input_path in progress.track(self.inputs, task_id=task1):
                task2 = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
                with io.BufferedReader(progress.open(input_path, "rb", task_id=task2)) as reader:  # type: ignore
                    if reader.peek().startswith(_GZIP_MAGIC):
                        reader = gzip.GzipFile(mode="rb", fileobj=reader)  # type: ignore
                    for record in gb_io.iter(reader):
                        if record.name in done:
                            n_duplicate += 1
                        else:
                            self.write_fasta(
                                dst, record.name, record.sequence.decode("ascii")
                            )
                            data.append((record.name, len(record.sequence), input_path))
                            done.add(record.name)
                progress.remove_task(task2)
            progress.remove_task(task1)
        if n_duplicate > 0:
            progress.console.print(
                f"[bold yellow]{'Skipped':>12}[/] {n_duplicate} clusters with duplicate identifiers"
            )
        return pd.DataFrame(
            data=data, columns=["cluster_id", "cluster_length", "filename"]
        ).set_index("cluster_id")

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
        representatives: typing.Container[str],
    ) -> typing.Dict[str, int]:
        """Extracts protein sequences from GenBank files.
        Args:
            progress (rich.progress.Progress): Progress bar for tracking progress.
            inputs (typing.List[pathlib.Path]): List of input GenBank files.
            output (pathlib.Path): Output file path for the extracted protein sequences.
            representatives (typing.Container[str]): Set of representative cluster IDs.
        Returns:
            typing.Dict[str, int]: Dictionary containing protein IDs and their sizes.
        """
        protein_sizes = {}
        with output.open("w") as dst:
            for input_path in self.inputs:
                task = progress.add_task(f"[bold blue]{'Reading':>9}[/]")
                with io.BufferedReader(progress.open(input_path, "rb", task_id=task)) as reader:  # type: ignore
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
                                if protein_id not in protein_sizes:
                                    self.write_fasta(dst, protein_id, translation)
                                    protein_sizes[protein_id] = len(translation)
                progress.remove_task(task)
        progress.console.print(
            f"[bold green]{'Extracted':>12}[/] {len(protein_sizes):,} proteins from {len(representatives):,} nucleotide representative"
        )
        return protein_sizes
