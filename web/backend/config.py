import os
from pathlib import Path

JOBS_DIR = Path(os.environ.get("PAIRMAP_JOBS_DIR", Path(__file__).parent.parent / "jobs"))
CACHE_DB = os.environ.get("PAIRMAP_CACHE_DB", None)
AVAILABLE_ENGINES = ["v2"]
DEFAULT_ENGINE = "v2"
