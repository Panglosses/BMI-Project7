"""
结果格式化器
"""

from typing import List, Dict, Any
from core.data_models import AccessibilityResult


class ResultFormatter:
    """结果格式化器"""

    @staticmethod
    def to_dict_list(results: List[AccessibilityResult]) -> List[Dict[str, Any]]:
        """转换为字典列表"""
        return [result.to_dict() for result in results]

    @staticmethod
    def to_simple_table(results: List[AccessibilityResult]) -> List[List[Any]]:
        """转换为简单表格格式（用于CSV）"""
        table = []
        for result in results:
            row = [
                result.residue.chain,
                result.residue.resnum,
                result.residue.resname,
                f"{result.min_distance:.3f}",
                result.water_count,
                "Yes" if result.accessible else "No",
            ]
            table.append(row)
        return table

    @staticmethod
    def create_comparison_table(
        custom_results: List[AccessibilityResult],
        sasa_results: List[Dict[str, Any]],
        match_ratio: float,
    ) -> List[List[Any]]:
        """
        创建对比表格

        Args:
            custom_results: 自定义方法结果
            sasa_results: FreeSASA结果（字典列表）
            match_ratio: 匹配比例

        Returns:
            List[List[Any]]: 对比表格
        """
        # 构建SASA结果映射
        sasa_map = {}
        for item in sasa_results:
            chain = item.get("chain", "").strip() or "A"
            resnum = str(item.get("resnum", ""))
            accessible = item.get("Accessible", "No")
            sasa_map[(chain, resnum)] = accessible

        # 构建对比表格
        comparison = []
        for result in custom_results:
            key = (result.residue.chain, str(result.residue.resnum))
            sasa_accessible = sasa_map.get(key, "No")
            match = (
                "Match"
                if result.accessible == (sasa_accessible == "Yes")
                else "Mismatch"
            )

            comparison.append(
                [
                    result.residue.chain,
                    result.residue.resnum,
                    result.residue.resname,
                    "Yes" if result.accessible else "No",
                    sasa_accessible,
                    match,
                ]
            )

        # 添加空行和统计信息
        comparison.append(["", "", "", "", "", ""])
        comparison.append(["Match_Ratio", f"{match_ratio:.4f}"])

        return comparison

    @staticmethod
    def format_summary(results: List[AccessibilityResult]) -> str:
        """格式化摘要信息"""
        total = len(results)
        accessible = sum(1 for r in results if r.accessible)
        ratio = accessible / total if total > 0 else 0.0

        summary = [
            "=== 分析结果摘要 ===",
            f"总残基数: {total}",
            f"可及残基数: {accessible}",
            f"可及比例: {ratio:.2%}",
            f"使用方法: {results[0].method if results else 'N/A'}",
        ]

        return "\n".join(summary)
