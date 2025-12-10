"""
Per-atom method implementation
"""

import numpy as np

from core.data_models import (
    ResidueInfo,
    WaterInfo,
    AccessibilityResult,
    AnalysisConfig,
    MethodType,
)
from core.distance_calculator import PerAtomDistanceCalculator
from core.accessibility_evaluator import PerAtomEvaluator


class PerAtomMethod:
    """Per-atom method"""

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: Analysis configuration
        """
        self.config = config or AnalysisConfig()
        self.distance_calculator = PerAtomDistanceCalculator(
            chunk_size=self.config.chunk_size, num_processes=self.config.num_processes
        )
        self.evaluator = PerAtomEvaluator()

    def analyze(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        structure,
    ) -> list[AccessibilityResult]:
        """
        Perform per-atom method analysis

        Args:
            residues: Residue list
            waters: Water molecule information
            structure: BioPython structure object

        Returns:
            list[AccessibilityResult]: Accessibility results
        """
        if not residues:
            return []

        # Calculate centroid distances (for fast filtering)
        min_distances = self.distance_calculator.compute_min_distances(residues, waters)

        # Count water molecules within radius
        water_counts = self.distance_calculator.count_waters_within_radius(
            residues, waters, self.config.radius
        )

        # Collect atom distances
        atom_distances = self.distance_calculator.collect_atom_distances(
            residues, waters, structure
        )

        # Set atom distance cache
        self.evaluator.set_atom_distances(atom_distances)

        # Evaluate accessibility
        results = self.evaluator.evaluate(
            residues, min_distances, water_counts, self.config
        )

        return results

    def get_method_type(self) -> MethodType:
        """Get method type"""
        return MethodType.PERATOM

    def get_atom_distances(self, residue: ResidueInfo) -> np.ndarray:
        """
        Get atom distances for specified residue

        Args:
            residue: Residue information

        Returns:
            np.ndarray: Atom distance array
        """
        key = (residue.chain, str(residue.resnum))
        return self.evaluator._atom_distances_cache.get(key, np.array([np.inf]))
