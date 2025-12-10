"""
可及性评估器
"""

from abc import ABC, abstractmethod
import numpy as np

from .data_models import ResidueInfo, AccessibilityResult, AnalysisConfig, MethodType


class AccessibilityEvaluator(ABC):
    """可及性评估器抽象基类"""

    @abstractmethod
    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """
        评估残基的可及性

        Args:
            residues: 残基列表
            min_distances: 最小距离数组
            water_counts: 水分子数量数组
            config: 分析配置

        Returns:
            list[AccessibilityResult]: 可及性结果列表
        """
        pass


class CentroidEvaluator(AccessibilityEvaluator):
    """质心法评估器"""

    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """质心法评估"""
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
    """原子级方法评估器"""

    def __init__(self):
        self._atom_distances_cache: dict[tuple, np.ndarray] = {}

    def set_atom_distances(self, atom_distances: dict[tuple, np.ndarray]):
        """设置原子距离缓存"""
        self._atom_distances_cache = atom_distances

    def evaluate(
        self,
        residues: list[ResidueInfo],
        min_distances: np.ndarray,
        water_counts: np.ndarray,
        config: AnalysisConfig,
    ) -> list[AccessibilityResult]:
        """原子级方法评估"""
        results = []
        for i, residue in enumerate(residues):
            key = (residue.chain, str(residue.resnum))
            atom_dists = self._atom_distances_cache.get(key, np.array([np.inf]))

            # 计算原子级可及性
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
        原子级可及性判断逻辑

        Args:
            residue: 残基信息
            atom_distances: 原子距离数组
            centroid_distance: 质心距离
            config: 分析配置

        Returns:
            bool: 是否可及
        """
        n_atoms = len(atom_distances)
        if n_atoms == 0:
            return False

        # 质心距离过大时直接判定为不可及
        if centroid_distance > (config.threshold + config.margin):
            return False

        # 统计命中原子数
        n_hits = int((atom_distances <= config.threshold).sum())

        # 判断是否为小残基
        is_small = (
            residue.resname.upper() in config.small_residues
            or n_atoms <= config.small_residue_size
        )

        if is_small:
            # 小残基：只需满足最小命中数
            return n_hits >= config.min_hits
        else:
            # 普通残基：需同时满足比例和最小命中数
            fraction = float(n_hits) / float(n_atoms)
            return fraction >= config.fraction_threshold and n_hits >= config.min_hits


class EvaluatorFactory:
    """评估器工厂"""

    @staticmethod
    def create_evaluator(
        method: MethodType,
        atom_distances: dict[tuple, np.ndarray] | None = None,
    ) -> AccessibilityEvaluator:
        """
        创建评估器

        Args:
            method: 方法类型
            atom_distances: 原子距离字典（仅peratom方法需要）

        Returns:
            AccessibilityEvaluator: 评估器实例
        """
        if method == MethodType.CENTROID:
            return CentroidEvaluator()
        elif method == MethodType.PERATOM:
            evaluator = PerAtomEvaluator()
            if atom_distances is not None:
                evaluator.set_atom_distances(atom_distances)
            return evaluator
        else:
            raise ValueError(f"未知的方法类型: {method}")
