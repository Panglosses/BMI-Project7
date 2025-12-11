# import solvent_analysis

import sys
import os

# Add the current directory to the path so that submodules can be imported.
current_dir = os.path.dirname(__file__)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Import and re-export the core data model.
try:
    from core.data_models import (
        ResidueInfo,
        WaterInfo,
        AccessibilityResult,
        AnalysisConfig,
        MethodType,
    )

    # Make these classes available at the module level.
    __all__ = [
        "ResidueInfo",
        "WaterInfo",
        "AccessibilityResult",
        "AnalysisConfig",
        "MethodType",
    ]

except ImportError as e:
    print(f"Warning: Failed to import the core data model: {e}")
    raise
