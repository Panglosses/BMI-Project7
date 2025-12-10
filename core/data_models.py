"""
核心数据模型定义
"""

from dataclasses import dataclass, field
from enum import Enum
import numpy as np


class MethodType(str, Enum):
    """分析方法类型"""

    CENTROID = "centroid"
    PERATOM = "peratom"


@dataclass
class ResidueInfo:
    """残基信息"""

    chain: str
    resnum: int
    resname: str
    coord: np.ndarray  # 质心坐标

    def __post_init__(self):
        """数据验证"""
        if not isinstance(self.chain, str):
            raise ValueError(f"chain必须是字符串，得到: {type(self.chain)}")
        if not isinstance(self.resnum, int):
            raise ValueError(f"resnum必须是整数，得到: {type(self.resnum)}")
        if not isinstance(self.resname, str):
            raise ValueError(f"resname必须是字符串，得到: {type(self.resname)}")
        if not isinstance(self.coord, np.ndarray):
            raise ValueError(f"coord必须是numpy数组，得到: {type(self.coord)}")
        if self.coord.shape != (3,):
            raise ValueError(f"coord必须是3维向量，形状: {self.coord.shape}")


@dataclass
class WaterInfo:
    """水分子信息"""

    coords: np.ndarray  # 水分子氧原子坐标，形状 (n, 3)
    names: list[str] = field(default_factory=list)  # 水分子名称

    def __post_init__(self):
        """数据验证"""
        if not isinstance(self.coords, np.ndarray):
            raise ValueError(f"coords必须是numpy数组，得到: {type(self.coords)}")
        if len(self.coords.shape) != 2 or self.coords.shape[1] != 3:
            raise ValueError(f"coords形状必须是(n, 3)，得到: {self.coords.shape}")
        if self.names and len(self.names) != len(self.coords):
            raise ValueError(
                f"names长度必须与coords匹配: {len(self.names)} != {len(self.coords)}"
            )

    @property
    def count(self) -> int:
        """水分子数量"""
        return len(self.coords)

    def is_empty(self) -> bool:
        """是否为空"""
        return self.count == 0


@dataclass
class AccessibilityResult:
    """可及性分析结果"""

    residue: ResidueInfo
    min_distance: float  # 最小距离（Å）
    water_count: int  # 半径R内的水分子数量
    accessible: bool  # 是否可及
    method: MethodType  # 使用方法

    def to_dict(self) -> dict[str, object]:
        """转换为字典"""
        return {
            "chain": self.residue.chain,
            "resnum": self.residue.resnum,
            "resname": self.residue.resname,
            "min_distance": self.min_distance,
            "water_count": self.water_count,
            "accessible": self.accessible,
            "method": self.method.value,
        }


@dataclass
class AnalysisConfig:
    """分析配置"""

    # 距离阈值
    threshold: float = 3.5  # 可及性判断阈值（Å）
    margin: float = 2.0  # 质心法的额外裕度（Å）
    radius: float = 5.0  # 统计水分子的半径（Å）

    # 原子级方法参数
    fraction_threshold: float = 0.20  # 原子可及比例阈值
    min_hits: int = 1  # 最小命中原子数
    small_residue_size: int = 5  # 小残基的原子数阈值

    # 计算参数
    chunk_size: int = 5000  # 分块计算大小
    num_processes: int = 1  # 并行进程数

    # 小残基集合
    small_residues: tuple[str, ...] = ("GLY", "ALA", "SER", "THR", "CYS", "PRO")

    # FreeSASA参数
    sasa_threshold: float = 10.0  # FreeSASA可及性阈值

    def validate(self):
        """验证配置参数"""
        if self.threshold <= 0:
            raise ValueError(f"threshold必须大于0: {self.threshold}")
        if self.radius <= 0:
            raise ValueError(f"radius必须大于0: {self.radius}")
        if not 0 <= self.fraction_threshold <= 1:
            raise ValueError(
                f"fraction_threshold必须在0-1之间: {self.fraction_threshold}"
            )
        if self.min_hits < 0:
            raise ValueError(f"min_hits必须非负: {self.min_hits}")
        if self.chunk_size <= 0:
            raise ValueError(f"chunk_size必须大于0: {self.chunk_size}")
        if self.num_processes <= 0:
            raise ValueError(f"num_processes必须大于0: {self.num_processes}")

        # 验证小残基名称
        for res in self.small_residues:
            if not isinstance(res, str) or len(res) != 3:
                raise ValueError(f"小残基名称必须是3字母代码: {res}")
