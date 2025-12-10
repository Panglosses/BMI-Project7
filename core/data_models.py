"""
Core data model definitions
"""

from dataclasses import dataclass, field
from enum import Enum
import numpy as np


class MethodType(str, Enum):
    """Analysis method type"""

    CENTROID = "centroid"
    PERATOM = "peratom"


@dataclass
class ResidueInfo:
    """Residue information"""

    chain: str
    resnum: int
    resname: str
    coord: np.ndarray  # Centroid coordinates

    def __post_init__(self):
        """Data validation"""
        if not isinstance(self.chain, str):
            raise ValueError(f"chain must be a string, got: {type(self.chain)}")
        if not isinstance(self.resnum, int):
            raise ValueError(f"resnum must be an integer, got: {type(self.resnum)}")
        if not isinstance(self.resname, str):
            raise ValueError(f"resname must be a string, got: {type(self.resname)}")
        if not isinstance(self.coord, np.ndarray):
            raise ValueError(f"coord must be a numpy array, got: {type(self.coord)}")
        if self.coord.shape != (3,):
            raise ValueError(f"coord must be a 3D vector, shape: {self.coord.shape}")


@dataclass
class WaterInfo:
    """Water molecule information"""

    coords: np.ndarray  # Water molecule oxygen atom coordinates, shape (n, 3)
    names: list[str] = field(default_factory=list)  # Water molecule names

    def __post_init__(self):
        """Data validation"""
        if not isinstance(self.coords, np.ndarray):
            raise ValueError(f"coords must be a numpy array, got: {type(self.coords)}")
        if len(self.coords.shape) != 2 or self.coords.shape[1] != 3:
            raise ValueError(f"coords shape must be (n, 3), got: {self.coords.shape}")
        if self.names and len(self.names) != len(self.coords):
            raise ValueError(
                f"names length must match coords: {len(self.names)} != {len(self.coords)}"
            )

    @property
    def count(self) -> int:
        """Water molecule count"""
        return len(self.coords)

    def is_empty(self) -> bool:
        """Whether empty"""
        return self.count == 0


@dataclass
class AccessibilityResult:
    """Accessibility analysis result"""

    residue: ResidueInfo
    min_distance: float  # Minimum distance (Å)
    water_count: int  # Water molecule count within radius R
    accessible: bool  # Whether accessible
    method: MethodType  # Method used

    def to_dict(self) -> dict[str, object]:
        """Convert to dictionary"""
        return {
            "chain": self.residue.chain,
            "resnum": self.residue.resnum,
            "resname": self.residue.resname,
            "min_distance": self.min_distance,
            "water_count": self.water_count,
            "accessible": self.accessible,
            "method": self.method.value,
        }


@dataclass
class AnalysisConfig:
    """Analysis configuration"""

    # Distance thresholds
    threshold: float = 3.5  # Accessibility threshold (Å)
    margin: float = 2.0  # Extra margin for centroid method (Å)
    radius: float = 5.0  # Radius for counting water molecules (Å)

    # Per-atom method parameters
    fraction_threshold: float = 0.20  # Atom accessibility fraction threshold
    min_hits: int = 1  # Minimum hit atoms
    small_residue_size: int = 5  # Small residue atom count threshold

    # Computation parameters
    chunk_size: int = 5000  # Chunk size for computation
    num_processes: int = 1  # Number of parallel processes

    # Small residue set
    small_residues: tuple[str, ...] = ("GLY", "ALA", "SER", "THR", "CYS", "PRO")

    # FreeSASA parameters
    sasa_threshold: float = 10.0  # FreeSASA accessibility threshold

    def validate(self):
        """Validate configuration parameters"""
        if self.threshold <= 0:
            raise ValueError(f"threshold must be greater than 0: {self.threshold}")
        if self.radius <= 0:
            raise ValueError(f"radius must be greater than 0: {self.radius}")
        if not 0 <= self.fraction_threshold <= 1:
            raise ValueError(
                f"fraction_threshold must be between 0 and 1: {self.fraction_threshold}"
            )
        if self.min_hits < 0:
            raise ValueError(f"min_hits must be non-negative: {self.min_hits}")
        if self.chunk_size <= 0:
            raise ValueError(f"chunk_size must be greater than 0: {self.chunk_size}")
        if self.num_processes <= 0:
            raise ValueError(f"num_processes must be greater than 0: {self.num_processes}")

        # Validate small residue names
        for res in self.small_residues:
            if not isinstance(res, str) or len(res) != 3:
                raise ValueError(f"Small residue name must be a 3-letter code: {res}")
