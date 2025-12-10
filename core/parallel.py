"""
并行计算组件
利用Python 3.14自由线程特性（PEP 703）实现KDTree查询并行化
"""

import numpy as np
from scipy.spatial import KDTree
import concurrent.futures
import threading  # 好像不能删
import time  # 好像不能删，被别人调用的时候有用？需要验证一下


class ParallelKDTreeQuery:
    """
    并行KDTree查询器

    在构建好KDTree后，对多个查询点进行并行最近邻搜索。
    利用Python 3.14自由线程特性（禁用GIL）实现真正的线程级并行。
    """

    def __init__(self, tree: KDTree, num_workers: int = 1):
        """
        Args:
            tree: 构建好的KDTree对象（只读，线程安全）
            num_workers: 并行工作线程数，默认为1（串行）
        """
        self.tree = tree
        self.num_workers = num_workers

    def query_ball_point_parallel(
        self, points: np.ndarray, radius: float
    ) -> list[list[int]]:
        """
        并行半径查询

        Args:
            points: 查询点数组，形状 (n_points, 3)
            radius: 查询半径

        Returns:
            每个查询点的邻居索引列表
        """
        if self.num_workers <= 1:
            # 串行版本
            return [self.tree.query_ball_point(p, radius) for p in points]

        n_points = len(points)
        results: list[list[int]] = [[] for _ in range(n_points)]

        # 使用ThreadPoolExecutor，利用Python 3.14自由线程特性
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_workers, thread_name_prefix="kdtree_worker"
        ) as executor:
            # 提交任务
            future_to_idx = {
                executor.submit(self.tree.query_ball_point, points[i], radius): i
                for i in range(n_points)
            }

            # 收集结果
            for future in concurrent.futures.as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    results[idx] = future.result()
                except Exception as e:
                    raise RuntimeError(f"并行查询失败 (点 {idx}): {e}")

        return results

    def query_nearest_parallel(self, points: np.ndarray, k: int = 1) -> tuple:
        """
        并行最近邻查询

        Args:
            points: 查询点数组，形状 (n_points, 3)
            k: 返回的最近邻数量

        Returns:
            (distances, indices): 距离和索引数组
        """
        if self.num_workers <= 1:
            # 串行版本
            return self.tree.query(points, k=k)

        n_points = len(points)
        # 根据k值确定数组形状
        if k == 1:
            distances = np.zeros(n_points)
            indices = np.zeros(n_points, dtype=int)
        else:
            distances = np.zeros((n_points, k))
            indices = np.zeros((n_points, k), dtype=int)

        # 分块并行处理
        chunk_size = max(1, n_points // self.num_workers)

        def process_chunk(start: int, end: int):
            """处理一个数据块"""
            chunk_points = points[start:end]
            if len(chunk_points) == 0:
                return start, end, None, None

            d, idx = self.tree.query(chunk_points, k=k)
            return start, end, d, idx

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_workers, thread_name_prefix="kdtree_nearest"
        ) as executor:
            # 提交分块任务
            futures = []
            for start in range(0, n_points, chunk_size):
                end = min(start + chunk_size, n_points)
                futures.append(executor.submit(process_chunk, start, end))

            # 收集结果
            for future in concurrent.futures.as_completed(futures):
                start, end, d, idx = future.result()
                if d is not None:
                    # 处理不同形状的数组
                    if d.ndim == 1:
                        # 一维数组 (k=1 情况)
                        distances[start:end] = d
                        indices[start:end] = idx
                    else:
                        # 二维数组 (k>1 情况)
                        distances[start:end, :] = d
                        indices[start:end, :] = idx

        return distances, indices


class ParallelDistanceMixin:
    """
    距离计算并行化混入类

    可拔插的并行化组件，通过继承或组合方式添加到现有计算器。
    """

    def __init__(self, num_processes: int = 1, **kwargs):
        """
        Args:
            num_processes: 并行进程数（实际使用线程）
            **kwargs: 传递给父类的其他参数
        """
        super().__init__(**kwargs)  # type: ignore
        self.num_processes = num_processes
        self._parallel_executor: ParallelKDTreeQuery | None = None

    def _get_parallel_executor(self, tree: KDTree) -> ParallelKDTreeQuery:
        """获取或创建并行查询器"""
        if self._parallel_executor is None or self._parallel_executor.tree is not tree:
            self._parallel_executor = ParallelKDTreeQuery(tree, self.num_processes)
        return self._parallel_executor

    def count_waters_within_radius_parallel(
        self,
        residues_coords: np.ndarray,
        water_tree: KDTree,
        radius: float,
    ) -> np.ndarray:
        """
        并行统计半径内的水分子数量

        Args:
            residues_coords: 残基坐标数组，形状 (n_residues, 3)
            water_tree: 水分子KDTree
            radius: 统计半径

        Returns:
            水分子数量数组
        """
        if self.num_processes <= 1:
            # 回退到串行版本
            counts = np.zeros(len(residues_coords), dtype=int)
            for i, coord in enumerate(residues_coords):
                indices = water_tree.query_ball_point(coord, radius)
                counts[i] = len(indices)
            return counts

        # 使用并行查询器
        executor = self._get_parallel_executor(water_tree)
        neighbor_lists = executor.query_ball_point_parallel(residues_coords, radius)

        # 转换为计数数组
        counts = np.array([len(neighbors) for neighbors in neighbor_lists], dtype=int)
        return counts

    def query_atom_distances_parallel(
        self,
        atom_coords_list: list[np.ndarray],
        water_tree: KDTree,
    ) -> list[np.ndarray]:
        """
        并行查询原子距离

        Args:
            atom_coords_list: 原子坐标列表，每个元素形状 (n_atoms, 3)
            water_tree: 水分子KDTree

        Returns:
            每个残基的原子距离列表
        """
        if self.num_processes <= 1 or len(atom_coords_list) < 2:
            # 串行版本
            results = []
            for atom_coords in atom_coords_list:
                if water_tree is not None:
                    dists, _ = water_tree.query(atom_coords, k=1)
                    results.append(np.array(dists, dtype=float))
                else:
                    results.append(np.full(len(atom_coords), np.inf))
            return results

        # 并行处理每个残基的原子
        executor = self._get_parallel_executor(water_tree)

        def process_single(atoms: np.ndarray) -> np.ndarray:
            """处理单个残基的原子"""
            if water_tree is not None:
                dists, _ = water_tree.query(atoms, k=1)
                return np.array(dists, dtype=float)
            else:
                return np.full(len(atoms), np.inf)

        # 使用线程池并行处理
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_processes, thread_name_prefix="atom_query"
        ) as pool:
            futures = [pool.submit(process_single, atoms) for atoms in atom_coords_list]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]

        # 保持原始顺序
        return results
