"""
溶剂可及性分析工具包
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Project"

from .core.data_models import (
    ResidueInfo,
    WaterInfo,
    AccessibilityResult,
    AnalysisConfig,
    MethodType
)

__all__ = [
    "ResidueInfo",
    "WaterInfo",
    "AccessibilityResult",
    "AnalysisConfig",
    "MethodType",
]