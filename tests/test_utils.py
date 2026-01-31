import datetime
import unittest

from igua._utils import Stopwatch


class TestStopwatch(unittest.TestCase):

    def test_error_not_started(self):
        stopwatch = Stopwatch()
        self.assertRaises(RuntimeError, stopwatch.total)

    def test_error_not_stopped(self):
        stopwatch = Stopwatch()
        stopwatch.start()
        self.assertRaises(RuntimeError, stopwatch.total)

    def test_total_human(self):
        stopwatch = Stopwatch()
        stopwatch.start_time = 1
        # stopwatch.start_datetime = datetime.datetime.now()
        stopwatch.stop_time = stopwatch.start_time + 2*3600 + 1
        # stopwatch.stop_datetime = datetime.datetime.now()
        self.assertEqual(stopwatch.total_human(), "2h 0m 1s")

        stopwatch.start_time = 1
        # stopwatch.start_datetime = datetime.datetime.now()
        stopwatch.stop_time = stopwatch.start_time + 2*60 + 1
        # stopwatch.stop_datetime = datetime.datetime.now()
        self.assertEqual(stopwatch.total_human(), "2m 1s")