"""
Parallel computation components
Implements KDTree query parallelization using Python 3.14 free-threading feature (PEP 703)
"""

import numpy as np
from scipy.spatial import KDTree
import concurrent.futures
import threading  # Seems cannot be removed
import time  # Seems cannot be removed, might be used by others? Need to verify


class ParallelKDTreeQuery:
    """
    Parallel KDTree query

    Performs parallel nearest neighbor search for multiple query points after KDTree construction.
    Uses Python 3.14 free-threading feature (GIL disabled) for true thread-level parallelism.
    """

    def __init__(self, tree: KDTree, num_workers: int = 1):
        """
        Args:
            tree: Constructed KDTree object (read-only, thread-safe)
            num_workers: Number of parallel worker threads, default 1 (serial)
        """
        self.tree = tree
        self.num_workers = num_workers

    def query_ball_point_parallel(
        self, points: np.ndarray, radius: float
    ) -> list[list[int]]:
        """
        Parallel radius query

        Args:
            points: Query point array, shape (n_points, 3)
            radius: Query radius

        Returns:
            Neighbor index list for each query point
        """
        if self.num_workers <= 1:
            # Serial version
            return [self.tree.query_ball_point(p, radius) for p in points]

        n_points = len(points)
        results: list[list[int]] = [[] for _ in range(n_points)]

        # Use ThreadPoolExecutor, leveraging Python 3.14 free-threading feature
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_workers, thread_name_prefix="kdtree_worker"
        ) as executor:
            # Submit tasks
            future_to_idx = {
                executor.submit(self.tree.query_ball_point, points[i], radius): i
                for i in range(n_points)
            }

            # Collect results
            for future in concurrent.futures.as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    results[idx] = future.result()
                except Exception as e:
                    raise RuntimeError(f"Parallel query failed (point {idx}): {e}")

        return results

    def query_nearest_parallel(self, points: np.ndarray, k: int = 1) -> tuple:
        """
        Parallel nearest neighbor query

        Args:
            points: Query point array, shape (n_points, 3)
            k: Number of nearest neighbors to return

        Returns:
            (distances, indices): Distance and index arrays
        """
        if self.num_workers <= 1:
            # Serial version
            return self.tree.query(points, k=k)

        n_points = len(points)
        # Determine array shape based on k value
        if k == 1:
            distances = np.zeros(n_points)
            indices = np.zeros(n_points, dtype=int)
        else:
            distances = np.zeros((n_points, k))
            indices = np.zeros((n_points, k), dtype=int)

        # Chunked parallel processing
        chunk_size = max(1, n_points // self.num_workers)

        def process_chunk(start: int, end: int):
            """Process a data chunk"""
            chunk_points = points[start:end]
            if len(chunk_points) == 0:
                return start, end, None, None

            d, idx = self.tree.query(chunk_points, k=k)
            return start, end, d, idx

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_workers, thread_name_prefix="kdtree_nearest"
        ) as executor:
            # Submit chunk tasks
            futures = []
            for start in range(0, n_points, chunk_size):
                end = min(start + chunk_size, n_points)
                futures.append(executor.submit(process_chunk, start, end))

            # Collect results
            for future in concurrent.futures.as_completed(futures):
                start, end, d, idx = future.result()
                if d is not None:
                    # Handle arrays of different shapes
                    if d.ndim == 1:
                        # 1D array (k=1 case)
                        distances[start:end] = d
                        indices[start:end] = idx
                    else:
                        # 2D array (k>1 case)
                        distances[start:end, :] = d
                        indices[start:end, :] = idx

        return distances, indices


class ParallelDistanceMixin:
    """
    Distance calculation parallelization mixin class

    Pluggable parallelization component that can be added to existing calculators via inheritance or composition.
    """

    def __init__(self, num_processes: int = 1, **kwargs):
        """
        Args:
            num_processes: Number of parallel processes (actually uses threads)
            **kwargs: Other arguments passed to parent class
        """
        super().__init__(**kwargs)  # type: ignore
        self.num_processes = num_processes
        self._parallel_executor: ParallelKDTreeQuery | None = None

    def _get_parallel_executor(self, tree: KDTree) -> ParallelKDTreeQuery:
        """Get or create parallel query"""
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
        Parallel count of water molecules within radius

        Args:
            residues_coords: Residue coordinate array, shape (n_residues, 3)
            water_tree: Water molecule KDTree
            radius: Counting radius

        Returns:
            Water molecule count array
        """
        if self.num_processes <= 1:
            # Fall back to serial version
            counts = np.zeros(len(residues_coords), dtype=int)
            for i, coord in enumerate(residues_coords):
                indices = water_tree.query_ball_point(coord, radius)
                counts[i] = len(indices)
            return counts

        # Use parallel query
        executor = self._get_parallel_executor(water_tree)
        neighbor_lists = executor.query_ball_point_parallel(residues_coords, radius)

        # Convert to count array
        counts = np.array([len(neighbors) for neighbors in neighbor_lists], dtype=int)
        return counts

    def query_atom_distances_parallel(
        self,
        atom_coords_list: list[np.ndarray],
        water_tree: KDTree,
    ) -> list[np.ndarray]:
        if self.num_processes <= 1 or len(atom_coords_list) < 2:
            # Serial version
            results = []
            for atom_coords in atom_coords_list:
                if water_tree is not None:
                    dists, _ = water_tree.query(atom_coords, k=1)
                    results.append(np.array(dists, dtype=float))
                else:
                    results.append(np.full(len(atom_coords), np.inf))
            return results

        # Process the atoms of each residue in parallel.
        executor = self._get_parallel_executor(water_tree)

        def process_single(atoms: np.ndarray) -> np.ndarray:
            """处理单个残基的原子"""
            if water_tree is not None:
                dists, _ = water_tree.query(atoms, k=1)
                return np.array(dists, dtype=float)
            else:
                return np.full(len(atoms), np.inf)

        # Use a thread pool for parallel processing.
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_processes, thread_name_prefix="atom_query"
        ) as pool:
            futures = [pool.submit(process_single, atoms) for atoms in atom_coords_list]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]

        return results
