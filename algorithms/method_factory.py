"""
Method factory
"""

from core.data_models import MethodType, AnalysisConfig
from algorithms.centroid_method import CentroidMethod
from algorithms.peratom_method import PerAtomMethod


class MethodFactory:
    """Method factory"""

    @staticmethod
    def create_method(
        method_type: MethodType | str,
        config: AnalysisConfig | None = None,
    ) -> CentroidMethod | PerAtomMethod:
        """
        Create analysis method

        Args:
            method_type: Method type (enum or string)
            config: Analysis configuration

        Returns:
            CentroidMethod | PerAtomMethod: Analysis method instance
        """
        # Handle string input
        if isinstance(method_type, str):
            method_type = MethodType(method_type.lower())

        # Ensure config is provided, because underlying constructor requires AnalysisConfig
        if config is None:
            raise ValueError("Analysis configuration config cannot be empty.")

        if method_type == MethodType.CENTROID:
            return CentroidMethod(config)
        elif method_type == MethodType.PERATOM:
            return PerAtomMethod(config)
        else:
            raise ValueError(f"Unknown method type: {method_type}")

    @staticmethod
    def get_available_methods() -> list:
        """Get available method list"""
        return [method.value for method in MethodType]
