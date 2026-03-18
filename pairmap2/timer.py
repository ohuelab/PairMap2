"""Stage timer context manager for PairMap2 pipeline profiling."""
import time
import os
import resource
from contextlib import contextmanager

from .types import StageTimings


def _get_rss_mb() -> float:
    try:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        # On macOS ru_maxrss is in bytes; on Linux it is in kilobytes
        if os.uname().sysname == "Darwin":
            return usage.ru_maxrss / 1024 / 1024
        else:
            return usage.ru_maxrss / 1024
    except Exception:
        return 0.0


@contextmanager
def stage_timer(stage_name: str, results_list: list):
    """Context manager that measures wall time, CPU time, and peak RSS."""
    t0_wall = time.perf_counter()
    t0_cpu = time.process_time()
    try:
        yield
    finally:
        wall = time.perf_counter() - t0_wall
        cpu = time.process_time() - t0_cpu
        rss = _get_rss_mb()
        results_list.append(
            StageTimings(
                stage=stage_name,
                wall_time=wall,
                cpu_time=cpu,
                peak_rss_mb=rss,
            )
        )
