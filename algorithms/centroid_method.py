"""
质心法实现
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
    """质心法"""

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: 分析配置
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
        structure=None,  # 保持接口一致，但质心法不需要structure
    ) -> list[AccessibilityResult]:
        """
        执行质心法分析

        Args:
            residues: 残基列表
            waters: 水分子信息
            structure: BioPython结构对象（可选）

        Returns:
            list[AccessibilityResult]: 可及性结果
        """
        # 验证输入
        if not residues:
            return []

        # 计算最小距离
        min_distances = self.distance_calculator.compute_min_distances(residues, waters)

        # 统计半径内的水分子数量
        water_counts = self.distance_calculator.count_waters_within_radius(
            residues, waters, self.config.radius
        )

        # 评估可及性
        results = self.evaluator.evaluate(
            residues, min_distances, water_counts, self.config
        )

        return results

    def get_method_type(self) -> MethodType:
        """获取方法类型"""
        return MethodType.CENTROID
