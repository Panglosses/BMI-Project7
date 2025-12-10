"""
FreeSASA wrapper
"""

import freesasa

from core.data_models import AnalysisConfig


class FreeSASAWrapper:
    """FreeSASA calculation wrapper"""

    # Water molecule name set (consistent with PDB loader)
    WATER_NAMES = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: Analysis configuration
        """
        self.config = config or AnalysisConfig()

    def compute_residue_sasa(self, pdb_file: str) -> list[dict[str, object]]:
        """
        Compute residue solvent accessible surface area

        Args:
            pdb_file: PDB file path

        Returns:
            list[dict[str, object]]: Residue SASA result list
                - chain: chain identifier
                - resnum: residue number
                - resname: residue name
                - SASA: solvent accessible surface area
                - Accessible: whether accessible (based on threshold)
        """
        try:
            structure = freesasa.Structure(pdb_file)
            result = freesasa.calc(structure)
            residue_areas = result.residueAreas()

            output = []
            for chain, chain_dict in residue_areas.items():
                for resnum, area_obj in chain_dict.items():
                    resname = area_obj.residueType

                    # Skip water molecules
                    if resname.upper() in self.WATER_NAMES:
                        continue

                    sasa = area_obj.total
                    accessible = "Yes" if sasa >= self.config.sasa_threshold else "No"

                    output.append(
                        {
                            "chain": chain,
                            "resnum": str(resnum),
                            "resname": resname,
                            "SASA": sasa,
                            "Accessible": accessible,
                        }
                    )

            return output

        except Exception as e:
            raise RuntimeError(f"FreeSASA calculation failed: {e}")

    @staticmethod
    def compute_simple(
        pdb_file: str, water_names: set, access_threshold: float
    ) -> list[dict[str, object]]:
        """
        Simple interface, maintaining compatibility with original compute_residue_sasa.py

        Args:
            pdb_file: PDB file path
            water_names: Water molecule name set
            access_threshold: Accessibility threshold

        Returns:
            list[dict[str, object]]: Residue SASA results
        """
        wrapper = FreeSASAWrapper()
        wrapper.WATER_NAMES = water_names
        wrapper.config.sasa_threshold = access_threshold
        return wrapper.compute_residue_sasa(pdb_file)
