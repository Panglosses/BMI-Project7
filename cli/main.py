"""
Main command line interface
"""

import argparse
import sys
import os
from pathlib import Path

current_dir = os.path.dirname(os.path.abspath(__file__))  # Locate to /cli
project_root = os.path.dirname(current_dir)              # Locate to program path

sys.path.insert(0, project_root) # Temporarily add the project root directory to the front of the Python module search path (sys.path)
from core.data_models import AnalysisConfig, MethodType
from io_utils import PDBLoader, CSVWriter, ResultFormatter
from algorithms import MethodFactory, FreeSASAWrapper


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Solvent accessibility analysis tool - based on water molecule proximity"
    )

    # Input files
    parser.add_argument(
        "--wet-pdb", required=True, help="Hydrated PDB file (for custom method analysis)"
    )
    parser.add_argument(
        "--dry-pdb", required=True, help="Dehydrated PDB file (for FreeSASA analysis)"
    )

    # Method selection
    parser.add_argument(
        "--method",
        choices=["centroid", "peratom"],
        default="peratom",
        help="Analysis method: centroid (centroid method) or peratom (per-atom method)",
    )

    # Distance parameters
    parser.add_argument(
        "--threshold", type=float, default=3.5, help="Accessibility threshold (Å), default: 3.5"
    )
    parser.add_argument(
        "--margin", type=float, default=2.0, help="Extra margin for centroid method (Å), default: 2.0"
    )
    parser.add_argument(
        "--R", type=float, default=5.0, help="Radius for counting water molecules (Å), default: 5.0"
    )

    # Per-atom method parameters
    parser.add_argument(
        "--fraction-threshold",
        type=float,
        default=0.20,
        help="Atom accessibility fraction threshold (0-1), default: 0.20",
    )
    parser.add_argument(
        "--min-hits", type=int, default=1, help="Minimum hit atoms, default: 1"
    )
    parser.add_argument(
        "--small-residue-size", type=int, default=5, help="Small residue atom count threshold, default: 5"
    )

    # Computation parameters
    parser.add_argument(
        "--chunk", type=int, default=5000, help="Chunk size for computation, default: 5000"
    )
    parser.add_argument("--nproc", type=int, default=1, help="Number of parallel processes, default: 1")

    # Output control
    parser.add_argument(
        "--output-dir", default="./output", help="Output directory, default: ./output"
    )
    parser.add_argument("--verbose", action="store_true", help="Show verbose output")
    parser.add_argument(
        "--no-comparison", action="store_true", help="Skip FreeSASA comparison"
    )

    return parser.parse_args(args)


def create_config(args) -> AnalysisConfig:
    """Create configuration from command line arguments"""
    config = AnalysisConfig(
        threshold=args.threshold,
        margin=args.margin,
        radius=args.R,
        fraction_threshold=args.fraction_threshold,
        min_hits=args.min_hits,
        small_residue_size=args.small_residue_size,
        chunk_size=args.chunk,
        num_processes=args.nproc,
    )
    config.validate()
    return config


def run_custom_analysis(args, config: AnalysisConfig):
    """Run custom method analysis"""
    if args.verbose:
        print(f"Loading PDB file: {args.wet_pdb}")

    # Load PDB
    loader = PDBLoader(quiet=not args.verbose)
    residues, waters, structure = loader.load(args.wet_pdb)

    if args.verbose:
        print(f"  Number of residues: {len(residues)}")
        print(f"  Number of water molecules: {waters.count}")

    # Create analysis method
    method = MethodFactory.create_method(args.method, config)

    # Execute analysis
    if args.verbose:
        print(f"Executing {args.method} analysis...")

    results = method.analyze(residues, waters, structure)

    if args.verbose:
        summary = ResultFormatter.format_summary(results)
        print(summary)

    return residues, results


def run_freesasa_analysis(args, config: AnalysisConfig):
    """Run FreeSASA analysis"""
    if args.verbose:
        print(f"Running FreeSASA analysis: {args.dry_pdb}")

    wrapper = FreeSASAWrapper(config)
    sasa_results = wrapper.compute_residue_sasa(args.dry_pdb)

    if args.verbose:
        accessible = sum(1 for r in sasa_results if r["Accessible"] == "Yes")
        total = len(sasa_results)
        print(f"  FreeSASA results: {accessible}/{total} accessible")

    return sasa_results


def save_results(args, custom_results, sasa_results=None):
    """Save result files"""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save custom method results
    prefix = Path(args.wet_pdb).stem
    custom_file = output_dir / f"{prefix}_{args.method}.csv"
    CSVWriter.write_results(str(custom_file), custom_results)

    if args.verbose:
        print(f"Saving custom method results: {custom_file}")

    # Save FreeSASA results
    if sasa_results:
        sasa_file = output_dir / f"{Path(args.dry_pdb).stem}_freesasa.csv"
        CSVWriter.write_generic(
            str(sasa_file),
            [
                [r["chain"], r["resnum"], r["resname"], r["SASA"], r["Accessible"]]
                for r in sasa_results
            ],
            ["chain", "resnum", "resname", "SASA", "Accessible"],
        )

        if args.verbose:
            print(f"Saving FreeSASA results: {sasa_file}")

    return custom_file


def compare_results(custom_results, sasa_results):
    """Compares custom method results with FreeSASA results and calculates the match ratio."""
    
    def normalize_chain(c):
        """Standardizes chain identifiers (logic identical to ResultFormatter._normalize_chain)."""
        if isinstance(c, str):
            c = c.strip()
        else:
            c = str(c).strip()
        return c if c else "A"
    
    sasa_map = {}
    for item in sasa_results:
        chain = normalize_chain(item.get("chain", ""))
        resnum = str(item.get("resnum", ""))
        accessible = item.get("Accessible", "No")
        sasa_map[(chain, resnum)] = accessible == "Yes"

    match_count = 0
    total = 0

    for result in custom_results:
        chain = normalize_chain(result.residue.chain)
        resnum = str(result.residue.resnum)
        key = (chain, resnum)
        
        sasa_accessible = sasa_map.get(key, False)
        
        if result.accessible == sasa_accessible:
            match_count += 1
        total += 1

    match_ratio = match_count / total if total > 0 else 0.0
    return match_ratio


def main(args=None):
    """Main function"""
    if args is None:
        args = parse_args()

    try:
        # Create configuration
        config = create_config(args)

        # Run custom method analysis
        residues, custom_results = run_custom_analysis(args, config)

        # Run FreeSASA analysis
        if not args.no_comparison:
            sasa_results = run_freesasa_analysis(args, config)

            # Compare results
            match_ratio = compare_results(custom_results, sasa_results)

            # Save comparison results
            comparison_file = Path(args.output_dir) / "comparison.csv"
            comparison_table = ResultFormatter.create_comparison_table(
                custom_results, sasa_results, match_ratio
            )
            CSVWriter.write_comparison(
                str(comparison_file),
                comparison_table,
                ["chain", "resnum", "resname", "Custom", "FreeSASA", "Match"],
            )

            if args.verbose:
                print(f"\n=== Comparison completed ===")
                print(f"Match ratio: {match_ratio:.4f}")
                print(f"Comparison results: {comparison_file}")

        # Save results
        custom_file = save_results(
            args, custom_results, sasa_results if not args.no_comparison else None
        )

        if args.verbose:
            print(f"\nAnalysis completed!")
            print(f"Result file: {custom_file}")
            if not args.no_comparison:
                print(f"Comparison file: {Path(args.output_dir) / 'comparison.csv'}")

        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
