import pathlib
import typing
from typing import Sequence, Iterable

import pandas
import rich.progress

from .base import BaseDataset


class DatasetList(BaseDataset, Sequence[BaseDataset]):
    """A dataset consisting of a list of other datasets.
    """

    def __init__(self, datasets: Iterable[BaseDataset] = ()):
        self.datasets = list(datasets)

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, index):
        return self.datasets[index]

    def __iter__(self):
        return iter(self.datasets)

    def extract_sequences(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
    ):
        task = progress.add_task(f"[bold blue]{'Working':>9}[/]")
        for dataset in progress.track(self.datasets, task_id=task):
            dataset.extract_sequences(progress, output)
        progress.remove_task(task)

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        output: pathlib.Path,
        representatives: typing.Container[str],
    ):
        task = progress.add_task(f"[bold blue]{'Working':>9}[/]")
        for dataset in progress.track(self.datasets, task_id=task):
            dataset.extract_proteins(progress, output, representatives)
        progress.remove_task(task)
                    