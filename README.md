# Solvent Accessibility Analysis Toolkit

High-efficiency calculation tool for protein residue solvent accessibility based on water molecule proximity.

## Project Introduction

Protein folding is driven by the hydrophobic effect, which buries hydrophobic residues in the protein core and exposes hydrophilic residues to the aqueous environment. Traditional methods based on geometric surface area calculation (e.g., FreeSASA), this method directly evaluates the proximity of residues to explicit water molecules, providing a more physically realistic solvent accessibility measure.

**Core Features**:
- Dual-algorithm support: centroid method (residue centroid distance) and per-atom method (atom contact ratio)
- High-performance computing: chunked computation, KDTree spatial indexing, vectorized operations
- Modular architecture: strategy pattern design, easy to extend with new algorithms
- Parallelization support: pluggable parallel components based on Python 3.14 free-threading feature

## Quick Start

### Installing Dependencies

1. Install [Python Install Manager](https://www.python.org/downloads/)

2. Install [Visual Studio 2026 Installer](https://visualstudio.microsoft.com/downloads/) and check the desktop C++ development environment

3. Run [setup.sh](./setup.sh) to automatically install the virtual environment (recommended to use Git Bash on Windows)

4. Activate the Python virtual environment `./bmi/`

### Adding Data

It is recommended to place the PDB files to be processed (e.g., `SUMO1_water.pdb`) inside `./pdb/`

```bash

Command line example:
python cli/main.py --wet-pdb "path/to/wet_structure.pdb" --dry-pdb "path/to/dry_structure.pdb" --method "peratom" --threshold 3.5 --margin 2.0 --R 5.0 --fraction-threshold 0.20 --min-hits 1 --small-residue-size 5 --chunk 5000 --nproc 4 --output-dir "./results" --verbose --no-comparison

explanation:
python cli/main.py \
  --wet-pdb  (Necessary) PDB file with water molecules
  --dry-pdb  (Necessary) Raw PDB file for FreeSASA
  --method  (Default peratom method) Analysis method: "centroid" or "peratom"
  --threshold 3.5 Accessibility threshold (Å): distance < threshold = accessible
  --margin 2.0 Centroid filter margin (Å): centroid distance > (threshold + margin) = filtered out
  --R 5.0  Radius for water counting (Å): count waters within the radius
  --fraction-threshold 0.20 Fraction threshold for peratom method (0-1): higher than this ratio of atoms within threshold is considered accessible
  --min-hits 1  Minimum atom hits for peratom method: at least this count of atoms within threshold is considered accessible
  --small-residue-size 5  Small residue threshold: residues with ≤ this count atoms apply 'min hit' rule while the bigger one applies 'fraction threshold' rule
  --chunk 5000  Chunk size for distance calculations: larger ---> more memory and faster
  --nproc 4  Number of parallel processes: use more cores for faster computation
  --output-dir  Output directory for results files
  --verbose Show detailed progress and statistics
  --no-comparison Skip FreeSASA comparison (dry-pdb still required but not used)
'''
```

### Python API

```python
from solvent_analysis import ResidueInfo, WaterInfo, AnalysisConfig, MethodType
from io_utils import PDBLoader
from algorithms import MethodFactory

loader = PDBLoader()
residues, waters, structure = loader.load("protein.pdb")

config = AnalysisConfig(threshold=3.5, radius=5.0, fraction_threshold=0.20)
method = MethodFactory.create_method(MethodType.PERATOM, config)
results = method.analyze(residues, waters, structure)

accessible = sum(1 for r in results if r.accessible)
print(f"Accessible residues: {accessible}/{len(results)}")
```

## Detailed Tutorial

For complete code explanation, algorithm principles, performance optimization, and parallelization design, please see:

**[tutorial.ipynb](tutorial.ipynb)** - Jupyter notebook detailed usage example, includes:
- Step-by-step code analysis and algorithm comparison
- Performance benchmark testing and parallelization analysis
- Interactive running and experimentation

## Project Structure

```
solvent_analysis/          # Main package
├── core/                  # Core interfaces and data models
├── io_utils/              # Input/output modules
├── algorithms/            # Algorithm implementations
├── utils/                 # Utility modules
└── cli/                   # Command line interface
```

**Core Modules**:
- `core/data_models.py` - Data class definitions (ResidueInfo, WaterInfo, etc.)
- `core/distance_calculator.py` - Distance calculation abstract interface
- `algorithms/centroid_method.py` - Centroid method implementation
- `algorithms/peratom_method.py` - Per-atom method implementation
- `io_utils/pdb_loader.py` - PDB file loading (BioPython integration)

## Testing and Validation

```bash
# Run unit tests
python test_refactored.py

# Performance benchmark testing
python benchmark_parallel.py --processes 1 2 4
```

## Performance Optimization

- **Chunked computation**: Controls memory peak (`chunk_size` parameter)
- **Spatial indexing**: KDTree accelerates nearest neighbor search
- **Parallel computation**: Configure `num_processes` parameter to utilize multiple cores
- **Vectorized operations**: NumPy broadcasting optimizes distance calculations

## Extension Development

The project uses strategy pattern and factory pattern design, supporting easy extension:
1. Implement `DistanceCalculator` abstract base class to add new distance algorithms
2. Implement `AccessibilityEvaluator` abstract base class to add new evaluation rules
3. Register new methods in `MethodFactory`

## Reference