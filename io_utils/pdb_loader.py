import numpy as np
from Bio.PDB import PDBParser  # type: ignore

from core.data_models import ResidueInfo, WaterInfo


class PDBLoader:
    WATER_NAMES = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}
    WATER_OXYGEN_NAMES = {"O", "OW", "OH", "OW1", "O1"}

    def __init__(self, quiet: bool = False):
        """
        Args:
            quiet: suspend BioPython warning
        """
        self.quiet = quiet

    def load(self, pdb_path: str) -> tuple[list[ResidueInfo], WaterInfo, object | None]:
        """
        load pdb file

        Returns:
            tuple[list[ResidueInfo], WaterInfo, object | None]:
                - Residue list
                - Water molecule information
                - BioPython structure object (may be None)
        """
        parser = PDBParser(QUIET=self.quiet)
        structure = parser.get_structure("prot", pdb_path)

        if structure is None:
            return [], WaterInfo(coords=np.empty((0, 3), dtype=float)), None

        residues = []
        water_coords = []
        water_names = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().upper().strip()
                    het_flag = residue.id[0]

                    if resname in self.WATER_NAMES:
                        self._extract_water_oxygen(residue, water_coords, water_names)
                        continue

                    if het_flag.strip():
                        continue

                    residue_info = self._extract_residue_info(chain, residue)
                    if residue_info is not None:
                        residues.append(residue_info)

        water_info = WaterInfo(
            coords=np.array(water_coords, dtype=float),
            names=water_names,
        )

        return residues, water_info, structure

    def _extract_water_oxygen(
        self,
        residue,
        water_coords: list,
        water_names: list,
    ) -> None:
        for atom in residue:
            element = getattr(atom, "element", "").upper()
            atom_name = atom.get_name().strip().upper()

            if atom_name in self.WATER_OXYGEN_NAMES or element == "O":
                water_coords.append(atom.coord)
                water_names.append(residue.get_resname())

    def _extract_residue_info(self, chain, residue) -> ResidueInfo | None:
        # Collecting non-H atom location
        atom_coords = []
        for atom in residue:
            element = getattr(atom, "element", "").upper()
            atom_name = atom.get_name().strip().upper()

            if element == "H" or atom_name.startswith("H"):
                continue

            atom_coords.append(atom.coord)

        if not atom_coords:
            return None

        centroid = np.mean(np.array(atom_coords, dtype=float), axis=0)

        return ResidueInfo(
            chain=chain.id,
            resnum=residue.id[1],
            resname=residue.get_resname().upper().strip(),
            coord=centroid,
        )


# Compatibility
def load_pdb(
    pdb_path: str, quiet: bool = False
) -> tuple[list[ResidueInfo], WaterInfo, object | None]:
    """
    Adding compatibility, maintaining the same interface as the original `load_pdb.py`

    Args:
        pdb_path: Path to the PDB file
        quiet: Whether to operate in silent mode

    Returns:
        tuple[list[ResidueInfo], WaterInfo, object | None]:
            - List of residues
            - Information about water molecules
            - BioPython structure object
    """
    loader = PDBLoader(quiet=quiet)
    return loader.load(pdb_path)
