from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any


@dataclass
class EngineResult:
    graphs: list  # list of nx.Graph
    node_mols: dict  # node_id -> RDKit Mol
    timings: dict = field(default_factory=dict)


class PairMapEngine(ABC):
    @abstractmethod
    def run(self, input_dir: str, config: dict) -> EngineResult:
        """Run from SDF files in input_dir."""
        ...

    @abstractmethod
    def run_from_moldf(self, mols: list, df: Any, config: dict) -> EngineResult:
        """Run from pre-built mol list and DataFrame."""
        ...
