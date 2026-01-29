import pathlib
import typing
from typing import Sequence, Iterable

import pandas
import rich.progress

from .base import BaseDataset, Cluster, Protein


class DatasetList(BaseDataset, Sequence[BaseDataset]):
    """A dataset consisting in a list of other datasets.
    """

    def __init__(self, datasets: Iterable[BaseDataset] = ()):
        super().__init__()
        self.datasets = list(datasets)

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, index):
        return self.datasets[index]

    def __iter__(self):
        return iter(self.datasets)

    def extract_clusters(
        self,
        progress: rich.progress.Progress,
    ) -> typing.Iterator[Cluster]:
        task = progress.add_task(f"[bold blue]{'Working':>9}[/]")
        for dataset in progress.track(self.datasets, task_id=task):
            yield from dataset.extract_clusters(progress)
        progress.remove_task(task)

    def extract_proteins(
        self,
        progress: rich.progress.Progress,
        representatives: typing.Container[str],
    ) -> typing.Iterator[Protein]:
        task = progress.add_task(f"[bold blue]{'Working':>9}[/]")
        for dataset in progress.track(self.datasets, task_id=task):
            yield from dataset.extract_proteins(progress, representatives)
        progress.remove_task(task)
                    