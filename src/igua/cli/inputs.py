import argparse
import pathlib
import typing

from ..dataset.base import BaseDataset
from ..dataset.genbank import GenBankDataset
from ..dataset.antismash import AntiSMASHGenBankDataset, AntiSMASHZipDataset
from ..dataset.defensefinder import DefenseFinderDataset


class BaseInput:
    
    def __init__(self, filename: str):
        self.filename = filename


class GenBankInput(BaseInput):
    
    def __init__(self, filename: str):
        self.filename = pathlib.Path(filename)

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return GenBankDataset(self.filename)


class GenBankListInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        datasets = []
        with open(self.filename) as files:
            for file in files:
                datasets.append(GenBankDataset(pathlib.Path(file.strip())))
        return DatasetList(datasets)


class AntiSMASHGenBankInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return AntiSMASHGenBankDataset(self.filename)


class AntiSMASHZipInput(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return AntiSMASHZipDataset(self.filename)


class DefenseFinderTSV(BaseInput):

    def to_dataset(self, args: argparse.Namespace) -> BaseDataset:
        return DefenseFinderDataset([self.filename])