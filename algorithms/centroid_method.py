"""
Centroid method implementation
"""

from core.data_models import (
    ResidueInfo,
    WaterInfo,
    AccessibilityResult,
    AnalysisConfig,
    MethodType,
)
from core.distance_calculator import ChunkedDistanceCalculator
from core.accessibility_evaluator import CentroidEvaluator


class CentroidMethod:
    """Centroid method"""

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: Analysis configuration
        """
        self.config = config or AnalysisConfig()
        self.distance_calculator = ChunkedDistanceCalculator(
            chunk_size=self.config.chunk_size, num_processes=self.config.num_processes
        )
        self.evaluator = CentroidEvaluator()

    def analyze(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        structure=None,  # Keep interface consistent, but centroid method does not need structure
    ) -> list[AccessibilityResult]:
        """
        Perform centroid method analysis

        Args:
            residues: Residue list
            waters: Water molecule information
            structure: BioPython structure object (optional)

        Returns:
            list[AccessibilityResult]: Accessibility results
        """
        # Validate input
        if not residues:
            return []

        # Calculate minimum distances
        min_distances = self.distance_calculator.compute_min_distances(residues, waters)

        # Count water molecules within radius
        water_counts = self.distance_calculator.count_waters_within_radius(
            residues, waters, self.config.radius
        )

        # Evaluate accessibility
        results = self.evaluator.evaluate(
            residues, min_distances, water_counts, self.config
        )

        return results

    def get_method_type(self) -> MethodType:
        """Get method type"""
        return MethodType.CENTROID
