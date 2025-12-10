"""
原子级方法实现
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
    """原子级方法"""

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: 分析配置
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
        执行原子级方法分析

        Args:
            residues: 残基列表
            waters: 水分子信息
            structure: BioPython结构对象

        Returns:
            list[AccessibilityResult]: 可及性结果
        """
        if not residues:
            return []

        # 计算质心距离（用于快速筛选）
        min_distances = self.distance_calculator.compute_min_distances(residues, waters)

        # 统计半径内的水分子数量
        water_counts = self.distance_calculator.count_waters_within_radius(
            residues, waters, self.config.radius
        )

        # 收集原子距离
        atom_distances = self.distance_calculator.collect_atom_distances(
            residues, waters, structure
        )

        # 设置原子距离缓存
        self.evaluator.set_atom_distances(atom_distances)

        # 评估可及性
        results = self.evaluator.evaluate(
            residues, min_distances, water_counts, self.config
        )

        return results

    def get_method_type(self) -> MethodType:
        """获取方法类型"""
        return MethodType.PERATOM

    def get_atom_distances(self, residue: ResidueInfo) -> np.ndarray:
        """
        获取指定残基的原子距离

        Args:
            residue: 残基信息

        Returns:
            np.ndarray: 原子距离数组
        """
        key = (residue.chain, str(residue.resnum))
        return self.evaluator._atom_distances_cache.get(key, np.array([np.inf]))
