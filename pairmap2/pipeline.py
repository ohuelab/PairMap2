"""PairMap2 Pipeline – orchestrates the full intermediate-insertion workflow.

Usage::

    from pairmap2 import Pipeline, PipelineConfig

    cfg = PipelineConfig(input_dir="./input", output_dir="./output")
    result = Pipeline(cfg).run()
    print(result.timings)
"""
import logging
from typing import Optional

from .intermediate_graph import IntermediateGraphManager

from .score_engine import ScoreEngine
from .score_cache import ScoreCache
from .timer import stage_timer
from .types import PipelineConfig, PipelineResult

logger = logging.getLogger(__name__)


class Pipeline:
    """Full PairMap2 pipeline with improved scoring and timing.

    Parameters
    ----------
    config:
        ``PipelineConfig`` instance.  If not given, a default config is built
        from any keyword arguments whose names match ``PipelineConfig`` fields.
    **kwargs:
        Keyword arguments forwarded to ``PipelineConfig`` when *config* is
        ``None``.
    """

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        **kwargs,
    ):
        if config is None:
            config_fields = {
                k: v
                for k, v in kwargs.items()
                if k in PipelineConfig.__dataclass_fields__
            }
            config = PipelineConfig(**config_fields)
        self.config = config

        cache = ScoreCache(config.cache_db_path)
        self.score_engine = ScoreEngine(
            cache=cache,
            atom_count_diff_threshold=config.atom_count_diff_threshold,
            tanimoto_prefilter=config.tanimoto_prefilter,
            jobs=config.jobs,
        )

        mgr_kwargs = dict(
            similarity_threshold=config.similarity_threshold,
            max_intermediate=config.max_intermediate,
            jobs=config.jobs,
            max=config.max_path_length,
            max_dist_from_actives=config.max_dist_from_actives,
            allow_tree=config.allow_tree,
            maxOptimalPathLength=config.max_optimal_path_length,
            roughScoreThreshold=config.rough_score_threshold,
            optimal_path_mode=config.optimal_path_mode,
            minScoreThreshold=config.min_score_threshold,
            verbose=config.verbose,
            lomap_options=config.lomap_options,
            output_dir=config.output_dir,
            save_output=config.save_output,
            input_dir=config.input_dir,
            is_atom_modfication_enabled=config.is_atom_modification_enabled,
            cap_ring_with_carbon=config.cap_ring_with_carbon,
            cap_ring_with_hydrogen=config.cap_ring_with_hydrogen,
            no_backward_search=config.no_backward_search,
            use_seed=config.use_seed,
            ionize=config.ionize,
        )

        self.mgr = IntermediateGraphManager(
            custom_get_score_matrix=self._score_matrix_wrapper,
            custom_get_similarity=self._similarity_wrapper,
            **mgr_kwargs,
        )

    def _score_matrix_wrapper(self, mols, options, jobs=None):
        """Adaptor between IntermediateGraphManager and ScoreEngine."""
        return self.score_engine.get_score_matrix(mols, self.config.lomap_options)

    def _similarity_wrapper(self, mol_a, mol_b, options=None):
        """Route get_similarity through ScoreEngine for shared cache."""
        return self.score_engine.get_score(mol_a, mol_b, self.config.lomap_options)

    def run(self, input_dir: Optional[str] = None) -> PipelineResult:
        """Run the full pipeline from an input directory."""
        timings: list = []

        if input_dir:
            self.mgr.config["input_dir"] = input_dir

        with stage_timer("generate_moldf", timings):
            mols, df = self.mgr.generate_moldf(self.mgr.config["input_dir"])

        logger.info(f"Number of bad edges: {len(df[df['BadEdge']])}")

        with stage_timer("run_from_moldf", timings):
            new_graphs, node_mols = self.mgr.run_from_moldf(mols, df)

        if self.mgr.config.get("save_output"):
            with stage_timer("save_output", timings):
                self.mgr.save_output(new_graphs, node_mols)

        return PipelineResult(graphs=new_graphs, node_mols=node_mols, timings=timings)

    def run_from_moldf(self, mols, df) -> PipelineResult:
        """Run from an already-prepared mols list and DataFrame (e.g., from benchmark scripts)."""
        timings: list = []
        with stage_timer("run_from_moldf", timings):
            new_graphs, node_mols = self.mgr.run_from_moldf(mols, df)
        return PipelineResult(graphs=new_graphs, node_mols=node_mols, timings=timings)
