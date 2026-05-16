# PairMap2

PairMap2 is a web application and Python package for intermediate insertion in
relative binding free energy (RBFE/FEP) workflows. Given source and target
ligands, it generates intermediate molecules, scores candidate transformations
with RDKit and LOMAP, and builds perturbation networks that can be inspected and
exported from a browser.

A public instance is available at <https://pairmap.yumizsui.com>.

PairMap2 prepares intermediate molecules and perturbation networks. It does not
run downstream FEP calculations.

## Features

- Pair mode for a single source-target transformation from SMILES or SDF input.
- Map mode for a multi-molecule SDF compound series.
- MCS-based edge inspection with common, deleted, and inserted atoms highlighted.
- Browser-based network visualization with Cytoscape.js and molecular rendering
  with RDKit.js and 3Dmol.js.
- Python API for scripted use.
- Export of generated intermediates and network links for downstream workflows.

## Requirements

- Python 3.12 or later
- A recent web browser for the web UI
- Linux, macOS, or Windows with a working Python environment

Main Python dependencies are listed in `pyproject.toml` and include RDKit,
LOMAP, NetworkX, FastAPI, pandas, NumPy, and Open Babel.

## Installation

Clone the repository:

```bash
git clone https://github.com/ohuelab/PairMap2.git
cd PairMap2
```

With `uv`:

```bash
uv sync
uv run uvicorn web.backend.main:app --host 127.0.0.1 --port 8000
```

With `pip`:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install ./local/gufe_stub -e .
uvicorn web.backend.main:app --host 127.0.0.1 --port 8000
```

Then open <http://127.0.0.1:8000>.

## Python API

Place input molecules in an input directory as SDF files, then run:

```python
from pairmap2 import Pipeline, PipelineConfig

config = PipelineConfig(
    input_dir="./input",
    output_dir="./output",
    save_output=True,
    similarity_threshold=0.6,
    max_path_length=4,
)

result = Pipeline(config).run()
print(result.timings)
```

When `save_output=True`, PairMap2 writes:

- `intermediate_graph.pkl`: graph history and RDKit molecule objects
- `intermediate_mols.sdf`: source, target, and generated intermediate molecules
- `intermediate_links.txt`: perturbation-network edges and similarity scores

## Configuration

Common `PipelineConfig` options include:

| Option | Default | Description |
| --- | ---: | --- |
| `similarity_threshold` | `0.6` | Edges below this score are treated as candidates for intermediate insertion. |
| `max_path_length` | `4` | Maximum allowed path length in the generated perturbation network. |
| `max_optimal_path_length` | `4` | Path-length limit used during map construction. |
| `max_intermediate` | `-1` | Maximum number of intermediate candidates; `-1` means no explicit limit. |
| `jobs` | `-1` | Number of worker processes for scoring; `-1` uses all available cores. |
| `ionize` | `True` | Apply ionization handling during intermediate generation. |
| `cache_db_path` | `None` | Optional SQLite score-cache path; `None` uses in-memory caching. |

## Citation

If you use the intermediate-insertion method, cite:

```bibtex
@article{Furui2025PairMap,
  author  = {Furui, Kairi and Shimizu, Takafumi and Akiyama, Yutaka and
             Kimura, S Roy and Terada, Yoh and Ohue, Masahito},
  title   = {{PairMap: An intermediate insertion approach for improving the
             accuracy of relative free energy perturbation calculations for
             distant compound transformations}},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2025},
  volume  = {65},
  number  = {2},
  pages   = {705--721},
  doi     = {10.1021/acs.jcim.4c01634}
}
```

## License

PairMap2 is distributed under the MIT License. See `LICENSE` for details.
