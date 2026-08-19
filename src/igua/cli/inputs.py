import abc
import argparse
import pathlib
import typing

from ..dataset.base import BaseDataset
from ..dataset.genbank import GenBankDataset
from ..dataset.antismash import AntiSMASHGenBankDataset, AntiSMASHZipDataset
from ..dataset.defensefinder import DefenseFinderDataset
from ..dataset.list import DatasetList


class BaseInput(abc.ABC):

    def __init__(self, filename: str):
        self.filename = filename

    @abc.abstractmethod
    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        pass


class ListInput(BaseInput):

    @abc.abstractmethod
    def single_dataset(self, filename: pathlib.Path):
        pass

    def to_dataset(self, args: argparse.Namespace) -> DatasetList[BaseDataset]:
        datasets = []
        with open(self.filename) as files:
            for file in map(str.strip, files):
                datasets.append(self.single_dataset(pathlib.Path(file)))
        return DatasetList(datasets)


class GenBankInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return GenBankDataset(self.filename)


class GenBankListInput(ListInput):

    def single_dataset(self, filename: pathlib.Path) -> BaseDataset:
        return GenBankDataset(self.filename)


class AntiSMASHGenBankInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return AntiSMASHGenBankDataset(self.filename)


class AntiSMASHGenBankListInput(ListInput):

    def single_dataset(self, filename: pathlib.Path) -> BaseDataset:
        return AntiSMASHGenBankDataset(filename)


class AntiSMASHZipInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return AntiSMASHZipDataset(self.filename)


class AntiSMASHZipListInput(ListInput):

    def single_dataset(self, filename: pathlib.Path) -> BaseDataset:
        return AntiSMASHZipDataset(filename)


class DefenseFinderTSV(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return DefenseFinderDataset([self.filename])
