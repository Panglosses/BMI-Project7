import csv
from pathlib import Path

from core.data_models import AccessibilityResult


class CSVWriter:
    @staticmethod
    def write_results(
        filepath: str,
        results: list[AccessibilityResult],
        include_header: bool = True,
    ) -> None:
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
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if header:
                writer.writerow(header)
            writer.writerows(data)
