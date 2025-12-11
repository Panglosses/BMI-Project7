"""
Distance calculator interface definitions
"""

from abc import ABC, abstractmethod
import numpy as np
from scipy.spatial import KDTree

from .data_models import ResidueInfo, WaterInfo
from .parallel import ParallelKDTreeQuery  # Parallel query component


class DistanceCalculator(ABC):
    """Distance calculator abstract base class"""

    @abstractmethod
    def compute_min_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
    ) -> np.ndarray:
        """
        Calculate the minimum distance from each residue to the nearest water molecule

        Args:
            residues: Residue list
            waters: Water molecule information

        Returns:
            np.ndarray: Minimum distance array, shape (n_residues,)
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
        Count the number of water molecules within radius R of each residue

        Args:
            residues: Residue list
            waters: Water molecule information
            radius: Counting radius (Å)

        Returns:
            np.ndarray: Water molecule count array, shape (n_residues,)
        """
        pass


class ChunkedDistanceCalculator(DistanceCalculator):
    """
    Chunked distance calculator
    Uses chunked computation to optimize memory usage
    """

    def __init__(self, chunk_size: int = 5000, num_processes: int = 1):
        """
        Args:
            chunk_size: Chunk size
            num_processes: Number of parallel processes (using threads), 1 indicates serial
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
        """Calculate minimum distances (chunked version)"""
        if waters.is_empty():
            return np.full(len(residues), np.inf)

        res_coords = np.vstack([r.coord for r in residues])
        water_coords = waters.coords

        n_res = len(residues)
        min_d2 = np.full(n_res, np.inf)
        n_w = len(water_coords)

        # Chunked computation
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
        """Count water molecules within radius"""
        if waters.is_empty():
            return np.zeros(len(residues), dtype=int)

        res_coords = np.vstack([r.coord for r in residues])
        water_coords = waters.coords

        # Build or reuse KDTree
        if (
            self._water_tree is None
            or self._water_coords is None
            or not np.array_equal(self._water_coords, water_coords)
        ):
            self._water_tree = KDTree(water_coords)
            self._water_coords = water_coords.copy()
            # Reset parallel query (tree changed)
            self._parallel_query = None

        # Select query method based on parallel process count
        if self.num_processes <= 1:
            # Serial query
            counts = np.zeros(len(residues), dtype=int)
            for i, coord in enumerate(res_coords):
                indices = self._water_tree.query_ball_point(coord, radius)
                counts[i] = len(indices)
        else:
            # Parallel query
            if (
                self._parallel_query is None
                or self._parallel_query.tree is not self._water_tree
            ):
                self._parallel_query = ParallelKDTreeQuery(
                    self._water_tree, self.num_processes
                )

            # Execute parallel radius query
            neighbor_lists = self._parallel_query.query_ball_point_parallel(
                res_coords, radius
            )
            counts = np.array(
                [len(neighbors) for neighbors in neighbor_lists], dtype=int
            )

        return counts


class PerAtomDistanceCalculator(DistanceCalculator):
    """
    Per-atom distance calculator
    Calculates distance from each non-hydrogen atom to the nearest water molecule
    """

    def __init__(self, chunk_size: int = 5000, num_processes: int = 1):
        """
        Args:
            chunk_size: Chunk size
            num_processes: Number of parallel processes (using threads), 1 indicates serial
        """
        self.chunk_size = chunk_size
        self.num_processes = num_processes

    def compute_min_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
    ) -> np.ndarray:
        """
        Per-atom minimum distance calculation
        Note: Returns the minimum distance from residue centroid to water molecule
        Atom-level distances are handled separately in collect_atom_distances
        """
        # Use centroid distance as base
        calculator = ChunkedDistanceCalculator(self.chunk_size, self.num_processes)
        return calculator.compute_min_distances(residues, waters)

    def count_waters_within_radius(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        radius: float,
    ) -> np.ndarray:
        """Count water molecules within radius"""
        calculator = ChunkedDistanceCalculator(self.chunk_size, self.num_processes)
        return calculator.count_waters_within_radius(residues, waters, radius)

    def collect_atom_distances(
        self,
        residues: list[ResidueInfo],
        waters: WaterInfo,
        structure,
    ) -> dict:
        """
        Collect distances for each atom

        Args:
            residues: Residue list
            waters: Water molecule information
            structure: BioPython structure object

        Returns:
            dict: Key is (chain, resnum), value is atom distance array
        """
        from scipy.spatial import KDTree
        import concurrent.futures

        if waters.is_empty():
            water_tree = None
        else:
            water_tree = KDTree(waters.coords)

        # Prepare parallel query (if needed)
        parallel_query = None
        if self.num_processes > 1 and water_tree is not None:
            from .parallel import ParallelKDTreeQuery

            parallel_query = ParallelKDTreeQuery(water_tree, self.num_processes)

        # Function to process a single residue
        def process_residue(r: ResidueInfo):
            """Atom distance calculation for a single residue"""
            try:
                residue = structure[0][r.chain][(" ", r.resnum, " ")]
            except Exception:
                key = (r.chain, str(r.resnum))
                return key, np.array([np.inf])

            # Collect non-hydrogen atom coordinates
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

            # Calculate distances
            if water_tree is not None:
                if parallel_query is not None:
                    # Use parallel query (for entire atom coordinate array)
                    dists, _ = parallel_query.query_nearest_parallel(atom_coords, k=1)
                    dists = np.array(dists.flatten(), dtype=float)
                else:
                    # Serial query
                    dists, _ = water_tree.query(atom_coords, k=1)
                    dists = np.array(dists, dtype=float)
            else:
                dists = np.full(len(atom_coords), np.inf)

            key = (r.chain, str(r.resnum))
            return key, dists

        # Select execution method based on parallel process count
        dists_map = {}
        if self.num_processes <= 1:
            # Serial processing
            for r in residues:
                key, dists = process_residue(r)
                dists_map[key] = dists
        else:
            # Parallel processing
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=self.num_processes, thread_name_prefix="atom_dist_worker"
            ) as executor:
                # Submit all tasks
                future_to_residue = {
                    executor.submit(process_residue, r): r for r in residues
                }

                # Collect results
                for future in concurrent.futures.as_completed(future_to_residue):
                    try:
                        key, dists = future.result()
                        dists_map[key] = dists
                    except Exception as e:
                        r = future_to_residue[future]
                        raise RuntimeError(f"Failed to process residue {r.chain}:{r.resnum}: {e}")

        return dists_map
