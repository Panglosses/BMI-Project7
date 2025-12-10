# import solvent_analysis

import sys
import os

# 将当前目录添加到路径，以便导入子模块
current_dir = os.path.dirname(__file__)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# 导入并重新导出核心数据模型
try:
    from core.data_models import (
        ResidueInfo,
        WaterInfo,
        AccessibilityResult,
        AnalysisConfig,
        MethodType,
    )

    # 使这些类在模块级别可用
    __all__ = [
        "ResidueInfo",
        "WaterInfo",
        "AccessibilityResult",
        "AnalysisConfig",
        "MethodType",
    ]

except ImportError as e:
    print(f"警告：导入核心数据模型失败: {e}")
    raise
