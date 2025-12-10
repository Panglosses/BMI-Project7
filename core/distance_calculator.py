"""
距离计算器接口定义
"""

from abc import ABC, abstractmethod
import numpy as np
from scipy.spatial import KDTree

from .data_models import ResidueInfo, WaterInfo
from .parallel import ParallelKDTreeQuery  # 并行查询组件


class DistanceCalculator(ABC):
    """距离计算器抽象基类"""

    @abstractmethod
    def compute_min_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
    ) -> np.ndarray:
        """
        计算每个残基到最近水分子的最小距离

        Args:
            residues: 残基列表
            waters: 水分子信息

        Returns:
            np.ndarray: 最小距离数组，形状 (n_residues,)
        """
        pass

    @abstractmethod
    def count_waters_within_radius(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        radius: float,
    ) -> np.ndarray:
        """
        统计每个残基半径R内的水分子数量

        Args:
            residues: 残基列表
            waters: 水分子信息
            radius: 统计半径（Å）

        Returns:
            np.ndarray: 水分子数量数组，形状 (n_residues,)
        """
        pass


class ChunkedDistanceCalculator(DistanceCalculator):
    """
    分块距离计算器
    使用分块计算优化内存使用
    """

    def __init__(self, chunk_size: int = 5000, num_processes: int = 1):
        """
        Args:
            chunk_size: 分块大小
            num_processes: 并行进程数（使用线程），1表示串行
        """
        self.chunk_size = chunk_size
        self.num_processes = num_processes
        self._water_tree: KDTree | None = None
        self._water_coords: np.ndarray | None = None
        self._parallel_query: ParallelKDTreeQuery | None = None

    def compute_min_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
    ) -> np.ndarray:
        """计算最小距离（分块版本）"""
        if waters.is_empty():
            return np.full(len(residues), np.inf)

        res_coords = np.vstack([r.coord for r in residues])
        water_coords = waters.coords

        n_res = len(residues)
        min_d2 = np.full(n_res, np.inf)
        n_w = len(water_coords)

        # 分块计算
        for start in range(0, n_w, self.chunk_size):
            end = min(start + self.chunk_size, n_w)
            diff = res_coords[:, None, :] - water_coords[None, start:end, :]
            d2 = np.sum(diff * diff, axis=2)
            min_d2 = np.minimum(min_d2, np.min(d2, axis=1))

        return np.sqrt(min_d2)

    def count_waters_within_radius(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        radius: float,
    ) -> np.ndarray:
        """统计半径内的水分子数量"""
        if waters.is_empty():
            return np.zeros(len(residues), dtype=int)

        res_coords = np.vstack([r.coord for r in residues])
        water_coords = waters.coords

        # 构建或重用KDTree
        if (
            self._water_tree is None
            or self._water_coords is None
            or not np.array_equal(self._water_coords, water_coords)
        ):
            self._water_tree = KDTree(water_coords)
            self._water_coords = water_coords.copy()
            # 重置并行查询器（树已改变）
            self._parallel_query = None

        # 根据并行进程数选择查询方式
        if self.num_processes <= 1:
            # 串行查询
            counts = np.zeros(len(residues), dtype=int)
            for i, coord in enumerate(res_coords):
                indices = self._water_tree.query_ball_point(coord, radius)
                counts[i] = len(indices)
        else:
            # 并行查询
            if (
                self._parallel_query is None
                or self._parallel_query.tree is not self._water_tree
            ):
                self._parallel_query = ParallelKDTreeQuery(
                    self._water_tree, self.num_processes
                )

            # 执行并行半径查询
            neighbor_lists = self._parallel_query.query_ball_point_parallel(
                res_coords, radius
            )
            counts = np.array(
                [len(neighbors) for neighbors in neighbor_lists], dtype=int
            )

        return counts


class PerAtomDistanceCalculator(DistanceCalculator):
    """
    原子级距离计算器
    计算每个非氢原子到最近水分子的距离
    """

    def __init__(self, chunk_size: int = 5000, num_processes: int = 1):
        """
        Args:
            chunk_size: 分块大小
            num_processes: 并行进程数（使用线程），1表示串行
        """
        self.chunk_size = chunk_size
        self.num_processes = num_processes

    def compute_min_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
    ) -> np.ndarray:
        """
        原子级最小距离计算
        注意：这里返回的是残基质心到水分子的最小距离
        原子级距离在collect_peratom_dists中单独处理
        """
        # 使用质心距离作为基础
        calculator = ChunkedDistanceCalculator(self.chunk_size, self.num_processes)
        return calculator.compute_min_distances(residues, waters)

    def count_waters_within_radius(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        radius: float,
    ) -> np.ndarray:
        """统计半径内的水分子数量"""
        calculator = ChunkedDistanceCalculator(self.chunk_size, self.num_processes)
        return calculator.count_waters_within_radius(residues, waters, radius)

    def collect_atom_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        structure,
    ) -> dict:
        """
        收集每个原子的距离

        Args:
            residues: 残基列表
            waters: 水分子信息
            structure: BioPython结构对象

        Returns:
            dict: 键为 (chain, resnum)，值为原子距离数组
        """
        from scipy.spatial import KDTree
        import concurrent.futures

        if waters.is_empty():
            water_tree = None
        else:
            water_tree = KDTree(waters.coords)

        # 准备并行查询器（如果需要）
        parallel_query = None
        if self.num_processes > 1 and water_tree is not None:
            from .parallel import ParallelKDTreeQuery

            parallel_query = ParallelKDTreeQuery(water_tree, self.num_processes)

        # 处理单个残基的函数
        def process_residue(r: ResidueInfo):
            """处理单个残基的原子距离计算"""
            try:
                residue = structure[0][r.chain][(" ", r.resnum, " ")]
            except Exception:
                key = (r.chain, str(r.resnum))
                return key, np.array([np.inf])

            # 收集非氢原子坐标
            atom_coords = []
            for atom in residue:
                elem = getattr(atom, "element", "").upper()
                aname = atom.get_name().strip().upper()
                if elem == "H" or aname.startswith("H"):
                    continue
                atom_coords.append(atom.coord)

            if not atom_coords:
                key = (r.chain, str(r.resnum))
                return key, np.array([np.inf])

            atom_coords = np.array(atom_coords, dtype=float)

            # 计算距离
            if water_tree is not None:
                if parallel_query is not None:
                    # 使用并行查询（针对整个原子坐标数组）
                    dists, _ = parallel_query.query_nearest_parallel(atom_coords, k=1)
                    dists = np.array(dists.flatten(), dtype=float)
                else:
                    # 串行查询
                    dists, _ = water_tree.query(atom_coords, k=1)
                    dists = np.array(dists, dtype=float)
            else:
                dists = np.full(len(atom_coords), np.inf)

            key = (r.chain, str(r.resnum))
            return key, dists

        # 根据并行进程数选择执行方式
        dists_map = {}
        if self.num_processes <= 1:
            # 串行处理
            for r in residues:
                key, dists = process_residue(r)
                dists_map[key] = dists
        else:
            # 并行处理
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=self.num_processes, thread_name_prefix="atom_dist_worker"
            ) as executor:
                # 提交所有任务
                future_to_residue = {
                    executor.submit(process_residue, r): r for r in residues
                }

                # 收集结果
                for future in concurrent.futures.as_completed(future_to_residue):
                    try:
                        key, dists = future.result()
                        dists_map[key] = dists
                    except Exception as e:
                        r = future_to_residue[future]
                        raise RuntimeError(f"处理残基 {r.chain}:{r.resnum} 失败: {e}")

        return dists_map
