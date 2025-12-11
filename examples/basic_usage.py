# Example

import sys
from pathlib import Path

# adding parent's path
sys.path.insert(0, str(Path(__file__).parent.parent))


def example_python_api():
    print("=== Python API Example ===")

    from solvent_analysis import ResidueInfo, WaterInfo, AnalysisConfig, MethodType
    from io_utils.pdb_loader import PDBLoader
    from algorithms import MethodFactory
    from io_utils.csv_writer import CSVWriter

    config = AnalysisConfig(
        threshold=3.5,  # Accessibility threshold (Å)
        radius=5.0,
        fraction_threshold=0.20,  # accessibility ratio
        min_hits=1,  # the minimun atom num that count as hit
        chunk_size=5000,
    )

    print(
        f"Config created successfully: threshold={config.threshold}, radius={config.radius}"
    )

    pdb_path = "../pdb/SUMO1_water.pdb"
    if not Path(pdb_path).exists():
        print(f"⚠ PDB file doesn't exist: {pdb_path}")
        return

    loader = PDBLoader(quiet=True)
    residues, waters, structure = loader.load(pdb_path)

    print(f"PDB load success:")
    print(f"  Num of residue: {len(residues)}")
    print(f"  Num of H2O: {waters.count}")
    print(f"  Example residue: {residues[0] if residues else '无'}")

    method = MethodFactory.create_method(MethodType.PERATOM, config)
    print(f"Analysis method created successfully: {method.get_method_type()}")

    print("run analysis...")
    results = method.analyze(residues, waters, structure)

    print(f"Analysis completed, number of results: {len(results)}")

    accessible = sum(1 for r in results if r.accessible)
    ratio = accessible / len(results) if results else 0

    print(f"Accessible residues: {accessible}/{len(results)} ({ratio:.1%})")

    print("\nThe first 5 results:")
    for i, result in enumerate(results[:5]):
        status = "Accessible" if result.accessible else "Not accessible"
        print(
            f"  {result.residue.chain}{result.residue.resnum} {result.residue.resname}: "
            f"Distance={result.min_distance:.2f}Å, H2O={result.water_count}, {status}"
        )

    output_file = "../output/example_results.csv"
    CSVWriter.write_results(output_file, results)
    print(f"\nResult saved to: {output_file}")

    return results


def example_configuration():
    """Usage example for configuration"""
    print("\n=== Configuration Usage Example ===")

    from solvent_analysis import AnalysisConfig

    # Default configuration
    default_config = AnalysisConfig()
    print(
        f"Default configuration: threshold={default_config.threshold}, radius={default_config.radius}"
    )

    # Custom configuration
    custom_config = AnalysisConfig(
        threshold=4.0,  # More lenient threshold
        radius=6.0,  # Larger statistical radius
        fraction_threshold=0.15,  # More lenient fraction threshold
        min_hits=2,  # Minimum of 2 atom hits
        small_residues=("GLY", "ALA", "SER"),  # Custom small residue set
    )

    print(
        f"Custom configuration: threshold={custom_config.threshold}, min_hits={custom_config.min_hits}"
    )

    # Validate configuration
    try:
        custom_config.validate()
        print("Configuration validation successful")
    except ValueError as e:
        print(f"Configuration validation failed: {e}")


def example_advanced_usage():
    """Advanced usage example"""
    print("\n=== Advanced Usage Example ===")

    from algorithms import FreeSASAWrapper
    from utils.progress import ProgressBar
    from utils.logger import setup_logger
    import logging

    # 1. Setup logging
    logger = setup_logger(level=logging.INFO, console=True)
    logger.info("Starting advanced example")

    # 2. Use FreeSASA
    pdb_path = "../pdb/SUMO1.pdb"
    if Path(pdb_path).exists():
        wrapper = FreeSASAWrapper()
        sasa_results = wrapper.compute_residue_sasa(pdb_path)
        logger.info(f"FreeSASA calculation completed: {len(sasa_results)} residues")
    else:
        logger.warning(f"PDB file does not exist: {pdb_path}")

    # 3. Use progress bar
    print("\nProgress bar example:")
    items = list(range(100))
    for item in ProgressBar.iterate(items, prefix="Processing", suffix="Complete"):
        # Simulate processing
        pass

    print("Advanced example completed")


def main():
    print("Solvent Accessibility Analysis Toolkit - Usage Examples")
    try:
        results = example_python_api()
        example_configuration()
        example_advanced_usage()

        print("\n" + "=" * 50)
        print("All examples executed successfully!")
        print("\nFor more usage, please refer to:")
        print("1. solvent_analysis/README.md - Module documentation")
        print("2. MIGRATION_GUIDE.md - Migration guide")
        print("3. test_refactored.py - Test script")

    except ImportError as e:
        print(f"Import error: {e}")
        print(
            "Please ensure you are running from the project root directory and the virtual environment is activated"
        )
    except Exception as e:
        print(f"Example execution error: {e}")

        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
