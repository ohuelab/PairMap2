"""Global ProcessPoolExecutor for CPU-bound scoring tasks."""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

_executor: ProcessPoolExecutor | None = None


def init(max_workers: int | None = None) -> None:
    global _executor
    workers = max_workers or int(os.environ.get("PAIRMAP_WORKERS", "0")) or min(4, os.cpu_count() or 2)
    _executor = ProcessPoolExecutor(max_workers=workers)


def shutdown() -> None:
    global _executor
    if _executor is not None:
        _executor.shutdown(wait=False)
        _executor = None


def get() -> ProcessPoolExecutor:
    if _executor is None:
        raise RuntimeError("Executor not initialized — call executor.init() first")
    return _executor
