import contextlib
import datetime
import io
import pathlib
import time
from typing import BinaryIO, Iterator

try:
    from isal import igzip as gzip
except ImportError:
    import gzip

try:
    import lz4.frame
except ImportError as err:
    lz4 = err

try:
    import bz2
except ImportError as err:
    bz2 = err

try:
    import lzma
except ImportError as err:
    lzma = err


class Stopwatch:
    """A stopwatch class to time execution of a context block.
    """

    def __init__(self):
        self.start_time = None
        self.start_datetime = None
        self.stop_time = None
        self.stop_datetime = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.stop()
        return False

    def start(self):
        """Start the stopwatch.
        """
        self.start_time = time.monotonic()
        self.start_datetime = datetime.datetime.now()

    def stop(self):
        """Stop the stopwatch.
        """
        self.stop_time = time.monotonic()
        self.stop_datetime = datetime.datetime.now()

    def total(self):
        """Compute the total execution time.

        Returns: 
            `float`: The total time the stopwatch has been active
            for. May have subsecond resolution if `time.monotonic`
            supports it on the current operating system.

        """
        if self.start_time is None:
            raise RuntimeError("stopwatch has not been started")
        if self.stop_time is None:
            raise RuntimeError("stopwatch has not been stopped")
        return self.stop_time - self.start_time

    def total_human(self):
        """Format the total execution time in a human-readable format.

        Returns:
            `str`: The total execution time summarized in hours, minutes,
            seconds depending on the total runtime.
            
        """
        total_seconds = int(self.total())

        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)

        if hours > 0:
            return f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            return f"{minutes}m {seconds}s"
        else:
            return f"{seconds}s"


_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"

@contextlib.contextmanager
def zopen(path: Union[str, pathlib.Path, BinaryIO]) -> Iterator[BinaryIO]:
    """Open a file with optional compression in binary mode.
    """
    with contextlib.ExitStack() as ctx:
        if isinstance(path, (str, pathlib.Path)):
            file = ctx.enter_context(open(path, "rb"))
        else:
            file = ctx.enter_context(io.BufferedReader(path))
        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            file = ctx.enter_context(gzip.open(file, mode="rb"))
        elif peek.startswith(_BZ2_MAGIC):
            if isinstance(bz2, ImportError):
                raise RuntimeError("File compression is LZMA but lzma is not available") from lz4
            file = ctx.enter_context(bz2.open(file, mode="rb"))
        elif peek.startswith(_XZ_MAGIC):
            if isinstance(lzma, ImportError):
                raise RuntimeError("File compression is LZMA but lzma is not available") from lz4
            file = ctx.enter_context(lzma.open(file, mode="rb"))
        elif peek.startswith(_LZ4_MAGIC):
            if isinstance(lz4, ImportError):
                raise RuntimeError("File compression is LZ4 but python-lz4 is not installed") from lz4
            file = ctx.enter_context(lz4.frame.open(file))
        yield file