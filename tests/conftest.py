from pathlib import Path
import pytest

BENCHMARK_DATA_DIR = Path(__file__).parents[1] / "benchmarks" / "data"

ALL_CASES = [
    "Bace1_00",
    "P38_00",
    "P38_01",
    "P38_02",
    "PTP1B_00",
    "PTP1B_01",
    "PTP1B_02",
    "PTP1B_03",
    "PTP1B_04",
]


@pytest.fixture
def benchmark_data_dir():
    return BENCHMARK_DATA_DIR


@pytest.fixture
def all_cases():
    return ALL_CASES
