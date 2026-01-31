import time
import datetime


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