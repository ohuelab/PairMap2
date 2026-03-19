from .base import PairMapEngine, EngineResult
from .default import PairMapDefaultEngine


def get_engine(name: str) -> PairMapEngine:
    engines = {
        "default": PairMapDefaultEngine,
    }
    if name not in engines:
        raise ValueError(f"Unknown engine: {name!r}. Available: {list(engines)}")
    return engines[name]()
