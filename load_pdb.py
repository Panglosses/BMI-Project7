import numpy as np
from Bio.PDB import PDBParser  # pyright: ignore[reportPrivateImportUsage]
from collections import namedtuple

ResidueInfo = namedtuple("ResidueInfo", ["chain", "resnum", "resname", "coord"])


def load_pdb(pdb_path):
    # parser = PDBParser(QUIET=True)
    parser = PDBParser(QUIET=False)
    struct = parser.get_structure("prot", pdb_path)
    if struct is None:
        return [], np.empty((0, 3), dtype=float), None

    residues = []
    water_coords = []
    water_names = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

    for model in struct:
        for chain in model:
            for residue in chain:
                rn = residue.get_resname().upper().strip()
                het = residue.id[0]

                # water detection
                if rn in water_names:
                    for atom in residue:
                        elem = getattr(atom, "element", "").upper()
                        aname = atom.get_name().strip().upper()
                        if aname in ("O", "OW", "OH", "OW1", "O1") or elem == "O":
                            water_coords.append(atom.coord)
                    continue

                # skip hetero (ligands etc.)
                if het.strip():
                    continue

                coords = [
                    atom.coord
                    for atom in residue
                    if getattr(atom, "element", "").upper() != "H"
                    and not atom.get_name().strip().startswith("H")
                ]
                if not coords:
                    continue
                centroid = np.mean(np.array(coords, dtype=float), axis=0)
                residues.append(ResidueInfo(chain.id, residue.id[1], rn, centroid))

    return residues, np.array(water_coords, dtype=float), struct


if __name__ == "__main__":

    pdb_path = "./pdb/KRAS_water.pdb"
    residues, water_coords, structure = load_pdb(pdb_path)

    print(f"残基数量: {len(residues)}")
    print(f"水分子坐标数量: {len(water_coords)}")
    print(f"结构对象: {'存在' if structure else '不存在'}")

    # 显示前几个残基的信息
    if residues:
        print(f"\n前5个残基信息:")
        for i, residue in enumerate(residues[:5]):
            print(
                f"  链: {residue.chain}, 残基号: {residue.resnum}, "
                f"残基名: {residue.resname}, 质心坐标: [{residue.coord[0]:.2f}, "
                f"{residue.coord[1]:.2f}, {residue.coord[2]:.2f}]"
            )

    # 显示水分子坐标
    if len(water_coords) > 0:
        print(f"\n前5个水分子坐标:")
        for i, coord in enumerate(water_coords[:5]):
            print(f"  {i+1}. [{coord[0]:.2f}, {coord[1]:.2f}, {coord[2]:.2f}]")

    # 显示结构信息
    if structure:
        print(f"\n结构信息:")
        print(f"模型数量: {len(list(structure.get_models()))}")
        print(f"链数量: {len(list(structure.get_chains()))}")
        print(f"原子总数: {len(list(structure.get_atoms()))}")
