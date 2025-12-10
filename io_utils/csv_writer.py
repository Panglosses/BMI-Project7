import csv
from pathlib import Path

from core.data_models import AccessibilityResult


class CSVWriter:
    """CSV文件写入器"""

    @staticmethod
    def write_results(
        filepath: str,
        results: list[AccessibilityResult],
        include_header: bool = True,
    ) -> None:
        """
        写入可及性分析结果

        Args:
            filepath: 输出文件路径
            results: 可及性结果列表
            include_header: 是否包含表头
        """
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)

            if include_header:
                header = [
                    "chain",
                    "resnum",
                    "resname",
                    "minDist_A",
                    "nWaterWithinR",
                    f"accessible_{results[0].method if results else 'unknown'}",
                ]
                writer.writerow(header)

            for result in results:
                row = [
                    result.residue.chain,
                    result.residue.resnum,
                    result.residue.resname,
                    f"{result.min_distance:.3f}",
                    result.water_count,
                    "Yes" if result.accessible else "No",
                ]
                writer.writerow(row)

    @staticmethod
    def write_comparison(
        filepath: str,
        comparison_data: list[list[object]],
        header: list[str],
    ) -> None:
        """
        写入对比结果

        Args:
            filepath: 输出文件路径
            comparison_data: 对比数据（每行一个列表）
            header: 表头
        """
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(comparison_data)

    @staticmethod
    def write_generic(
        filepath: str,
        data: list[list[object]],
        header: list[str] | None = None,
    ) -> None:
        """
        通用CSV写入

        Args:
            filepath: 输出文件路径
            data: 数据行列表
            header: 可选的表头
        """
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if header:
                writer.writerow(header)
            writer.writerows(data)
