import abc
import typing


class BaseRecordSink(abc.ABC):

    @abc.abstractmethod
    def add_record(self, name: str, sequence: str) -> None:
        pass


class FASTASink(BaseRecordSink):

    def __init__(self, file: typing.TextIO) -> None:
        self.file = file

    def add_record(self, name: str, sequence: str) -> None:
        self.file.write(">{}\n".format(name))
        self.file.write(sequence)
        self.file.write("\n")