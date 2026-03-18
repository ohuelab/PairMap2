from .base import PairMapEngine, EngineResult
from .v2 import PairMapV2Engine


def get_engine(name: str) -> PairMapEngine:
    engines = {
        "v2": PairMapV2Engine,
    }
    if name not in engines:
        raise ValueError(f"Unknown engine: {name!r}. Available: {list(engines)}")
    return engines[name]()
