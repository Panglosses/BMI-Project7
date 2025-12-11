"""
Accessibility evaluator
"""

from abc import ABC, abstractmethod
import numpy as np

from .data_models import ResidueInfo, AccessibilityResult, AnalysisConfig, MethodType


class AccessibilityEvaluator(ABC):
    """Accessibility evaluator abstract base class"""

    @abstractmethod
    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """
        Evaluate residue accessibility

        Args:
            residues: Residue list
            min_distances: Minimum distance array
            water_counts: Water molecule count array
            config: Analysis configuration

        Returns:
            list[AccessibilityResult]: Accessibility result list
        """
        pass


class CentroidEvaluator(AccessibilityEvaluator):
    """Centroid method evaluator"""

    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """Centroid method evaluation"""
        results = []
        for i, residue in enumerate(residues):
            accessible = min_distances[i] <= config.threshold
            result = AccessibilityResult(
                residue=residue,
                min_distance=float(min_distances[i]),
                water_count=int(water_counts[i]),
                accessible=accessible,
                method=MethodType.CENTROID,
            )
            results.append(result)
        return results


class PerAtomEvaluator(AccessibilityEvaluator):
    """Per-atom method evaluator"""

    def __init__(self):
        self._atom_distances_cache: dict[tuple, np.ndarray] = {}

    def set_atom_distances(self, atom_distances: dict[tuple, np.ndarray]):
        """Set atom distance cache"""
        self._atom_distances_cache = atom_distances

    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """Per-atom method evaluation"""
        results = []
        for i, residue in enumerate(residues):
            key = (residue.chain, str(residue.resnum))
            atom_dists = self._atom_distances_cache.get(key, np.array([np.inf]))

            # Calculate per-atom accessibility
            accessible = self._evaluate_per_atom(
                residue=residue,
                atom_distances=atom_dists,
                centroid_distance=min_distances[i],
                config=config,
            )

            result = AccessibilityResult(
                residue=residue,
                min_distance=float(min_distances[i]),
                water_count=int(water_counts[i]),
                accessible=accessible,
                method=MethodType.PERATOM,
            )
            results.append(result)
        return results

    def _evaluate_per_atom(
        self,
        residue: ResidueInfo,
        atom_distances: np.ndarray,
        centroid_distance: float,
        config: AnalysisConfig,
    ) -> bool:
        """
        Per-atom accessibility determination logic

        Args:
            residue: Residue information
            atom_distances: Atom distance array
            centroid_distance: Centroid distance
            config: Analysis configuration

        Returns:
            bool: Whether accessible
        """
        n_atoms = len(atom_distances)
        if n_atoms == 0:
            return False

        # If centroid distance is too large, directly judge as inaccessible
        if centroid_distance > (config.threshold + config.margin):
            return False

        # Count hit atoms
        n_hits = int((atom_distances <= config.threshold).sum())

        # Determine if it's a small residue
        is_small = (
            residue.resname.upper() in config.small_residues
            or n_atoms <= config.small_residue_size
        )

        if is_small:
            # Small residue: only need to meet minimum hit count
            return n_hits >= config.min_hits
        else:
            # Normal residue: need to satisfy both fraction and minimum hit count
            fraction = float(n_hits) / float(n_atoms)
            return fraction >= config.fraction_threshold and n_hits >= config.min_hits


class EvaluatorFactory:
    """Evaluator factory"""

    @staticmethod
    def create_evaluator(
        method: MethodType,
        atom_distances: dict[tuple, np.ndarray] | None = None,
    ) -> AccessibilityEvaluator:
        """
        Create evaluator

        Args:
            method: Method type
            atom_distances: Atom distance dictionary (only needed for peratom method)

        Returns:
            AccessibilityEvaluator: Evaluator instance
        """
        if method == MethodType.CENTROID:
            return CentroidEvaluator()
        elif method == MethodType.PERATOM:
            evaluator = PerAtomEvaluator()
            if atom_distances is not None:
                evaluator.set_atom_distances(atom_distances)
            return evaluator
        else:
            raise ValueError(f"Unknown method type: {method}")
