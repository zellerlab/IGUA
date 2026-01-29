import abc
import typing

import pandas


class BaseRecordSink(abc.ABC):

    def __init__(self):
        self.stats = []
        self.done = set()

    def add_record(self, name: str, sequence: str, **kwargs) -> bool:
        if name not in self.done:
            self.stats.append({
                "id": name,
                "length": len(sequence),
                **kwargs,
            })
            self.done.add(name)
            return True
        else:
            return False

    def report_statistic(self) -> pandas.DataFrame:
        return pandas.DataFrame(self.stats)


class FASTASink(BaseRecordSink):

    def __init__(self, file: typing.TextIO) -> None:
        super().__init__()
        self.file = file

    def add_record(self, name: str, sequence: str, **kwargs) -> None:
        if not super().add_record(name, sequence, **kwargs):
            return False

        self.file.write(">{}\n".format(name))
        self.file.write(sequence)
        self.file.write("\n")
        return True