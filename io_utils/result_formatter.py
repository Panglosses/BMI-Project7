from core.data_models import AccessibilityResult


class ResultFormatter:
    @staticmethod
    def to_dict_list(results: list[AccessibilityResult]) -> list[dict[str, object]]:
        return [result.to_dict() for result in results]

    @staticmethod
    def to_simple_table(results: list[AccessibilityResult]) -> list[list[object]]:
        """Convert to a simple table format (for CSV)"""
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
        Create a comparison table

        Args:
            custom_results: Results from the custom method
            sasa_results: FreeSASA results (a list of dictionaries)
            match_ratio: Matching ratio

        Returns:
            list[list[object]]: Comparison table
        """
        # Constructing the SASA result mapping
        sasa_map = {}
        for item in sasa_results:
            chain = str(item.get("chain", "")).strip() or "A"
            resnum = str(item.get("resnum", ""))
            accessible = str(item.get("Accessible", "No"))
            sasa_map[(chain, resnum)] = accessible

        # Create a comparison table.
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

        # Add blank lines and statistics
        comparison.append(["", "", "", "", "", ""])
        comparison.append(["Match_Ratio", f"{match_ratio:.4f}"])

        return comparison

    @staticmethod
    def format_summary(results: list[AccessibilityResult]) -> str:
        total = len(results)
        accessible = sum(1 for r in results if r.accessible)
        ratio = accessible / total if total > 0 else 0.0

        summary = [
            "=== Result Summary ===",
            f"Total residues: {total}",
            f"Accessible residues: {accessible}",
            f"Accessible ratio: {ratio:.2%}",
            f"Use case: {results[0].method if results else 'N/A'}",
        ]

        return "\n".join(summary)
