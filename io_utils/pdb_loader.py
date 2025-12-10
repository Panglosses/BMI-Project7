"""
PDB文件加载器
"""

import numpy as np
from Bio.PDB import PDBParser  # type: ignore

from core.data_models import ResidueInfo, WaterInfo


class PDBLoader:
    """PDB文件加载器"""

    # 水分子名称集合
    WATER_NAMES = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

    # 水分子氧原子名称
    WATER_OXYGEN_NAMES = {"O", "OW", "OH", "OW1", "O1"}

    def __init__(self, quiet: bool = False):
        """
        Args:
            quiet: 是否静默模式（抑制BioPython警告）
        """
        self.quiet = quiet

    def load(self, pdb_path: str) -> tuple[list[ResidueInfo], WaterInfo, object | None]:
        """
        加载PDB文件

        Args:
            pdb_path: PDB文件路径

        Returns:
            tuple[list[ResidueInfo], WaterInfo, object | None]:
                - 残基列表
                - 水分子信息
                - BioPython结构对象（可能为None）
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

                    # 水分子检测
                    if resname in self.WATER_NAMES:
                        self._extract_water_oxygen(residue, water_coords, water_names)
                        continue

                    # 跳过异质残基（配体等）
                    if het_flag.strip():
                        continue

                    # 提取残基信息
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
        """提取水分子氧原子坐标"""
        for atom in residue:
            element = getattr(atom, "element", "").upper()
            atom_name = atom.get_name().strip().upper()

            if atom_name in self.WATER_OXYGEN_NAMES or element == "O":
                water_coords.append(atom.coord)
                water_names.append(residue.get_resname())

    def _extract_residue_info(self, chain, residue) -> ResidueInfo | None:
        """提取残基信息"""
        # 收集非氢原子坐标
        atom_coords = []
        for atom in residue:
            element = getattr(atom, "element", "").upper()
            atom_name = atom.get_name().strip().upper()

            if element == "H" or atom_name.startswith("H"):
                continue

            atom_coords.append(atom.coord)

        if not atom_coords:
            return None

        # 计算质心
        centroid = np.mean(np.array(atom_coords, dtype=float), axis=0)

        return ResidueInfo(
            chain=chain.id,
            resnum=residue.id[1],
            resname=residue.get_resname().upper().strip(),
            coord=centroid,
        )


# 兼容性函数
def load_pdb(
    pdb_path: str, quiet: bool = False
) -> tuple[list[ResidueInfo], WaterInfo, object | None]:
    """
    兼容性函数，保持与原load_pdb.py相同的接口

    Args:
        pdb_path: PDB文件路径
        quiet: 是否静默模式

    Returns:
        tuple[list[ResidueInfo], WaterInfo, object | None]:
            - 残基列表
            - 水分子信息
            - BioPython结构对象
    """
    loader = PDBLoader(quiet=quiet)
    return loader.load(pdb_path)
