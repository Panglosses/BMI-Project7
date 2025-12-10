"""
结果格式化器
"""

from core.data_models import AccessibilityResult


class ResultFormatter:
    """结果格式化器"""

    @staticmethod
    def _normalize_chain(chain):
        """标准化链标识（与旧代码完全一致）"""
        if isinstance(chain, str):
            chain = chain.strip()
        else:
            chain = str(chain).strip()
        return chain if chain else "A"  # 空链转为"A"

    @staticmethod
    def to_dict_list(results: list[AccessibilityResult]) -> list[dict[str, object]]:
        """转换为字典列表"""
        return [result.to_dict() for result in results]

    @staticmethod
    def to_simple_table(results: list[AccessibilityResult]) -> list[list[object]]:
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
        custom_results: list[AccessibilityResult],
        sasa_results: list[dict[str, object]],
        match_ratio: float,
    ) -> list[list[object]]:
        """
        创建对比表格 - 使用旧代码逻辑
        """
        # 构建SASA结果映射（与旧代码一致）
        sasa_map = {}
        for item in sasa_results:
            chain = ResultFormatter._normalize_chain(item.get("chain", ""))
            resnum = str(item.get("resnum", ""))
            accessible = str(item.get("Accessible", "No"))
            sasa_map[(chain, resnum)] = accessible

        # 构建对比表格
        comparison = []
        for result in custom_results:
            # 使用相同的标准化函数
            chain = ResultFormatter._normalize_chain(result.residue.chain)
            resnum = str(result.residue.resnum)
            key = (chain, resnum)
            
            sasa_accessible = sasa_map.get(key, "No")
            custom_accessible = result.accessible
            
            # 对比逻辑（与旧代码一致）
            match = "Match" if custom_accessible == (sasa_accessible == "Yes") else "Mismatch"

            comparison.append(
                [
                    chain,  # 使用标准化后的chain
                    resnum,
                    result.residue.resname,
                    "Yes" if custom_accessible else "No",
                    sasa_accessible,
                    match,
                ]
            )

        # 添加空行和统计信息（与旧代码一致）
        comparison.append(["", "", "", "", "", ""])
        comparison.append(["Match_Ratio", f"{match_ratio:.4f}"])

        return comparison

    @staticmethod
    def format_summary(results: list[AccessibilityResult]) -> str:
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