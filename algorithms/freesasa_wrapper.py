"""
FreeSASA包装器
"""

import freesasa

from core.data_models import AnalysisConfig


class FreeSASAWrapper:
    """FreeSASA计算包装器"""

    # 水分子名称集合（与PDB加载器保持一致）
    WATER_NAMES = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

    def __init__(self, config: AnalysisConfig | None = None):
        """
        Args:
            config: 分析配置
        """
        self.config = config or AnalysisConfig()

    def compute_residue_sasa(self, pdb_file: str) -> list[dict[str, object]]:
        """
        计算残基的溶剂可及表面积

        Args:
            pdb_file: PDB文件路径

        Returns:
            list[dict[str, object]]: 残基SASA结果列表
                - chain: 链标识
                - resnum: 残基编号
                - resname: 残基名称
                - SASA: 溶剂可及表面积
                - Accessible: 是否可及（基于阈值）
        """
        try:
            structure = freesasa.Structure(pdb_file)
            result = freesasa.calc(structure)
            residue_areas = result.residueAreas()

            output = []
            for chain, chain_dict in residue_areas.items():
                for resnum, area_obj in chain_dict.items():
                    resname = area_obj.residueType

                    # 跳过水分子
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
            raise RuntimeError(f"FreeSASA计算失败: {e}")

    @staticmethod
    def compute_simple(
        pdb_file: str, water_names: set, access_threshold: float
    ) -> list[dict[str, object]]:
        """
        简单接口，保持与原compute_residue_sasa.py兼容

        Args:
            pdb_file: PDB文件路径
            water_names: 水分子名称集合
            access_threshold: 可及性阈值

        Returns:
            list[dict[str, object]]: 残基SASA结果
        """
        wrapper = FreeSASAWrapper()
        wrapper.WATER_NAMES = water_names
        wrapper.config.sasa_threshold = access_threshold
        return wrapper.compute_residue_sasa(pdb_file)
