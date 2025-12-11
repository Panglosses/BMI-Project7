from core.data_models import AccessibilityResult


class ResultFormatter:
    """Formats analysis results for output and comparison."""
    
    @staticmethod
    def _normalize_chain(chain):
        """
        Standardizes chain identifiers for consistent matching.
        
        Args:
            chain: Input chain identifier (can be string, None, or other types).
        
        Returns:
            str: A cleaned chain identifier. Empty or whitespace-only inputs are
                 defaulted to 'A'.
        """
        if isinstance(chain, str):
            chain = chain.strip()
        else:
            chain = str(chain).strip()
        # Ensure a non-empty default chain for matching purposes
        return chain if chain else "A"

    @staticmethod
    def to_dict_list(results: list[AccessibilityResult]) -> list[dict[str, object]]:
        """Converts a list of AccessibilityResult objects to a list of dictionaries."""
        return [result.to_dict() for result in results]

    @staticmethod
    def to_simple_table(results: list[AccessibilityResult]) -> list[list[object]]:
        """Converts results to a simple table format suitable for CSV export."""
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
        Creates a unified table comparing custom method results with FreeSASA results.
        
        This method ensures residues are matched correctly by standardizing chain
        identifiers from both data sources before performing the comparison.
        
        Args:
            custom_results: Results from the custom accessibility method.
            sasa_results: Results from FreeSASA (a list of dictionaries).
            match_ratio: The calculated match ratio between the two methods.
        
        Returns:
            list[list[object]]: A table where each row contains the residue info,
                                results from both methods, and a match status.
        """
        # Map FreeSASA results using standardized chain identifiers
        sasa_map = {}
        for item in sasa_results:
            chain = ResultFormatter._normalize_chain(item.get("chain", ""))
            resnum = str(item.get("resnum", ""))
            accessible = str(item.get("Accessible", "No"))
            sasa_map[(chain, resnum)] = accessible

        # Build the comparison table
        comparison = []
        for result in custom_results:
            # Standardize chain identifier for consistent lookup
            chain = ResultFormatter._normalize_chain(result.residue.chain)
            resnum = str(result.residue.resnum)
            lookup_key = (chain, resnum)
            
            sasa_accessible = sasa_map.get(lookup_key, "No")
            # Determine if the accessibility calls from both methods agree
            match_status = "Match" if result.accessible == (sasa_accessible == "Yes") else "Mismatch"

            comparison.append(
                [
                    chain,  # Use the standardized chain identifier
                    resnum,
                    result.residue.resname,
                    "Yes" if result.accessible else "No",
                    sasa_accessible,
                    match_status,
                ]
            )

        # Append a separator and the final match ratio
        comparison.append(["", "", "", "", "", ""])
        comparison.append(["Match_Ratio", f"{match_ratio:.4f}"])

        return comparison

    @staticmethod
    def format_summary(results: list[AccessibilityResult]) -> str:
        """Generates a textual summary of the analysis results."""
        total = len(results)
        accessible = sum(1 for r in results if r.accessible)
        ratio = accessible / total if total > 0 else 0.0

        summary = [
            "=== Result Summary ===",
            f"Total residues: {total}",
            f"Accessible residues: {accessible}",
            f"Accessible ratio: {ratio:.2%}",
            f"Method used: {results[0].method if results else 'N/A'}",
        ]
        return "\n".join(summary)