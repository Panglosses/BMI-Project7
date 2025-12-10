"""
方法工厂
"""

from core.data_models import MethodType, AnalysisConfig
from algorithms.centroid_method import CentroidMethod
from algorithms.peratom_method import PerAtomMethod


class MethodFactory:
    """方法工厂"""

    @staticmethod
    def create_method(
        method_type: MethodType | str,
        config: AnalysisConfig | None = None,
    ) -> CentroidMethod | PerAtomMethod:
        """
        创建分析方法

        Args:
            method_type: 方法类型（枚举或字符串）
            config: 分析配置

        Returns:
            CentroidMethod | PerAtomMethod: 分析方法实例
        """
        # 处理字符串输入
        if isinstance(method_type, str):
            method_type = MethodType(method_type.lower())

        # 确保提供了 config，因为底层构造函数需要 AnalysisConfig
        if config is None:
            raise ValueError("分析配置 config 不能为空。")

        if method_type == MethodType.CENTROID:
            return CentroidMethod(config)
        elif method_type == MethodType.PERATOM:
            return PerAtomMethod(config)
        else:
            raise ValueError(f"未知的方法类型: {method_type}")

    @staticmethod
    def get_available_methods() -> list:
        """获取可用方法列表"""
        return [method.value for method in MethodType]
