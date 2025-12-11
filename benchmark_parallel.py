"""
Parallel performance comparison script

Tests the time/space performance of the KDTree parallel query component:
1. Compare runtimes under different thread counts
2. Verify consistency between parallel and serial results
3. Analyze speedup and scalability
4. Measure memory usage changes

Proposition verification: The distance calculation between each amino acid and its nearest water molecule is independent,
thus can be parallelized after constructing the KDTree.
"""

import sys
import time
import argparse
from pathlib import Path
import numpy as np
from typing import (
    TypedDict,
)  # This is unavoidable for type safety, must be explicitly imported, and is the recommended approach for 3.12+, no need to change


class ProcessMetrics(TypedDict):
    """Process performance metrics"""

    mean: float
    std: float
    speedup: float
    efficiency: float


class BenchmarkResult(TypedDict):
    """Benchmark test results"""

    method: str
    num_residues: int
    num_waters: int
    process_times: dict[int, ProcessMetrics]
    result_consistency: bool
    serial_time: float
    serial_time_std: float


# Add project path to import modules
sys.path.insert(0, str(Path(__file__).parent))

from core.data_models import AnalysisConfig, MethodType
from io_utils.pdb_loader import PDBLoader
from algorithms.method_factory import MethodFactory


def benchmark_method(
    wet_pdb: str,
    dry_pdb: str,
    method_type: MethodType,
    num_processes_list: list[int],
    num_runs: int = 3,
) -> BenchmarkResult:
    """
    Benchmark testing method

    Args:
        wet_pdb: Hydrated PDB file path
        dry_pdb: Dehydrated PDB file path (currently unused, reserved for future FreeSASA comparison)
        method_type: Method type (centroid or per-atom)
        num_processes_list: List of parallel process counts to test
        num_runs: Number of runs per configuration (average taken)

    Returns:
        Dictionary containing performance results
    """
    print(f"\n{'='*60}")
    print(f"Benchmark test: {method_type.value} method")
    print(f"Test file: {Path(wet_pdb).name}")
    print(f"Process count list: {num_processes_list}")
    print(f"Number of runs: {num_runs}")
    print(f"{'='*60}")

    # Load PDB file
    loader = PDBLoader(quiet=True)
    residues, waters, structure = loader.load(wet_pdb)

    print(
        f"System size: {len(residues)} residues, {len(waters.coords) if waters.coords is not None else 0} water molecules"
    )

    results: BenchmarkResult = {
        "method": method_type.value,
        "num_residues": len(residues),
        "num_waters": len(waters.coords) if waters.coords is not None else 0,
        "process_times": {},
        "result_consistency": True,
        "serial_time": 0.0,
        "serial_time_std": 0.0,
    }

    # First run serial version as baseline and correctness reference
    print(f"\n1. Serial baseline (num_processes=1)")
    serial_config = AnalysisConfig(num_processes=1)
    serial_method = MethodFactory.create_method(method_type, serial_config)

    serial_times = []
    serial_results_list = []

    for run_idx in range(num_runs):
        start_time = time.perf_counter()
        serial_results = serial_method.analyze(residues, waters, structure)
        end_time = time.perf_counter()
        serial_times.append(end_time - start_time)
        serial_results_list.append(serial_results)

        if run_idx == 0:
            # Record serial results for comparison
            serial_accessible = sum(1 for r in serial_results if r.accessible)
            print(
                f"   Run {run_idx+1}: {serial_times[-1]:.3f}s, accessible residues: {serial_accessible}/{len(serial_results)}"
            )

    avg_serial_time = float(np.mean(serial_times))
    std_serial_time = float(np.std(serial_times))
    results["serial_time"] = avg_serial_time
    results["serial_time_std"] = std_serial_time
    print(f"   Average time: {avg_serial_time:.3f}s (±{std_serial_time:.3f}s)")

    # Test different parallel process counts
    for num_proc in num_processes_list:
        print(f"\n2. Parallel test (num_processes={num_proc})")

        # Configure parallel version
        parallel_config = AnalysisConfig(num_processes=num_proc)
        parallel_method = MethodFactory.create_method(method_type, parallel_config)

        parallel_times = []
        all_consistent = True

        for run_idx in range(num_runs):
            start_time = time.perf_counter()
            parallel_results = parallel_method.analyze(residues, waters, structure)
            end_time = time.perf_counter()
            parallel_times.append(end_time - start_time)

            # Verify result consistency (compared with first serial run)
            if run_idx == 0:
                parallel_accessible = sum(1 for r in parallel_results if r.accessible)
                serial_accessible = sum(
                    1 for r in serial_results_list[0] if r.accessible
                )

                # Detailed comparison (optional)
                if len(parallel_results) == len(serial_results_list[0]):
                    mismatches = 0
                    for p_res, s_res in zip(parallel_results, serial_results_list[0]):
                        if p_res.accessible != s_res.accessible:
                            mismatches += 1
                            if mismatches <= 3:  # Only print first few mismatches
                                print(
                                    f"     Warning: residue {p_res.residue.chain}:{p_res.residue.resnum} "
                                    f"accessibility mismatch (parallel: {p_res.accessible}, serial: {s_res.accessible})"
                                )

                    if mismatches > 0:
                        print(
                            f"     Total mismatches: {mismatches}/{len(parallel_results)}"
                        )
                        all_consistent = False
                    else:
                        print(f"     Result consistency: ✓ all matched")
                else:
                    print(
                        f"     Error: result lengths differ (parallel: {len(parallel_results)}, serial: {len(serial_results_list[0])})"
                    )
                    all_consistent = False

                print(
                    f"   Run {run_idx+1}: {parallel_times[-1]:.3f}s, accessible residues: {parallel_accessible}/{len(parallel_results)}"
                )
            else:
                print(f"   Run {run_idx+1}: {parallel_times[-1]:.3f}s")

        avg_parallel_time = float(np.mean(parallel_times))
        std_parallel_time = float(np.std(parallel_times))
        speedup = (
            float(avg_serial_time / avg_parallel_time) if avg_parallel_time > 0 else 0.0
        )

        results["process_times"][num_proc] = {
            "mean": avg_parallel_time,
            "std": std_parallel_time,
            "speedup": speedup,
            "efficiency": float(speedup / num_proc) if num_proc > 0 else 0.0,
        }

        if not all_consistent:
            results["result_consistency"] = False

        print(f"   Average time: {avg_parallel_time:.3f}s (±{std_parallel_time:.3f}s)")
        print(f"   Speedup: {speedup:.2f}x (efficiency: {speedup/num_proc*100:.1f}%)")

    return results


def print_summary_table(results_list: list[BenchmarkResult]) -> None:
    """Print performance summary table"""
    print(f"\n{'='*80}")
    print("Performance comparison summary")
    print(f"{'='*80}")

    for result in results_list:
        method = result["method"]
        print(f"\nMethod: {method.upper()}")
        print(
            f"System: {result['num_residues']} residues, {result['num_waters']} water molecules"
        )
        print(f"Result consistency: {'✓' if result['result_consistency'] else '✗'}")

        print(
            f"\n{'Processes':<8} {'Avg time(s)':<12} {'Std dev':<10} {'Speedup':<10} {'Efficiency(%)':<10}"
        )
        print(f"{'-'*50}")

        # Serial results
        print(
            f"{'1 (serial)':<8} {result['serial_time']:<12.3f} {result['serial_time_std']:<10.3f} {'1.00':<10} {'100.0':<10}"
        )

        # Parallel results
        for num_proc, metrics in sorted(result["process_times"].items()):
            print(
                f"{num_proc:<8} {metrics['mean']:<12.3f} {metrics['std']:<10.3f} "
                f"{metrics['speedup']:<10.2f} {metrics['efficiency']*100:<10.1f}"
            )


def analyze_scalability(results_list: list[BenchmarkResult]) -> None:
    """Analyze scalability characteristics"""
    print(f"\n{'='*80}")
    print("Scalability analysis")
    print(f"{'='*80}")

    for result in results_list:
        method = result["method"]
        print(f"\n{method.upper()} method:")

        if len(result["process_times"]) < 2:
            print("  Insufficient data for scalability analysis")
            continue

        # Calculate ideal speedup (Amdahl's law)
        # Assume parallel fraction p, serial fraction 1-p
        # Simple estimation: use speedup at maximum process count
        max_proc = max(result["process_times"].keys())
        max_speedup = result["process_times"][max_proc]["speedup"]

        if max_speedup > 1:
            # Infer parallel fraction from Amdahl's law
            # S = 1 / ((1-p) + p/N)
            # where S is speedup, N is process count
            # Solve for p = (1/S - 1) / (1/N - 1)
            p = (1 / max_speedup - 1) / (1 / max_proc - 1)
            serial_fraction = 1 - p

            print(f"  Maximum speedup: {max_speedup:.2f}x (with {max_proc} threads)")
            print(f"  Estimated parallel fraction: {p*100:.1f}%")
            print(f"  Serial bottleneck: {serial_fraction*100:.1f}%")

            # Determine scalability type
            efficiency = result["process_times"][max_proc]["efficiency"]
            if efficiency > 0.7:
                print(f"  Scalability: excellent (efficiency > 70%)")
            elif efficiency > 0.5:
                print(f"  Scalability: good (efficiency > 50%)")
            elif efficiency > 0.3:
                print(f"  Scalability: fair (efficiency > 30%)")
            else:
                print(f"  Scalability: poor (efficiency ≤ 30%)")
        else:
            print(f"  No speedup effect")


def main():
    parser = argparse.ArgumentParser(
        description="Parallel performance comparison test script"
    )
    parser.add_argument(
        "--wet-pdb",
        default="./pdb/SUMO1_water.pdb",
        help="Hydrated PDB file path (e.g., ./pdb/SUMO1_water.pdb)",
    )
    parser.add_argument(
        "--dry-pdb",
        default="./pdb/SUMO1.pdb",
        help="Dehydrated PDB file path (default: ./pdb/SUMO1.pdb)",
    )
    parser.add_argument(
        "--processes",
        type=int,
        nargs="+",
        default=[1, 2, 4, 8],
        help="List of parallel process counts to test (default: 1 2 4 8)",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=3,
        help="Number of runs per configuration (default: 3)",
    )
    parser.add_argument(
        "--methods",
        choices=["centroid", "peratom", "both"],
        default="both",
        help="Methods to test (default: both)",
    )

    args = parser.parse_args()

    # Verify file existence
    if not Path(args.wet_pdb).exists():
        print(f"Error: file does not exist: {args.wet_pdb}")
        sys.exit(1)

    if not Path(args.dry_pdb).exists():
        print(
            f"Warning: file does not exist: {args.dry_pdb} (FreeSASA comparison will be skipped)"
        )

    # Determine methods to test
    if args.methods == "both":
        method_types = [MethodType.CENTROID, MethodType.PERATOM]
    elif args.methods == "centroid":
        method_types = [MethodType.CENTROID]
    else:
        method_types = [MethodType.PERATOM]

    # Run benchmark tests
    all_results: list[BenchmarkResult] = []

    for method_type in method_types:
        try:
            result = benchmark_method(
                args.wet_pdb,
                args.dry_pdb,
                method_type,
                args.processes,
                args.runs,
            )
            all_results.append(result)
        except Exception as e:
            print(f"Error testing method {method_type.value}: {e}")
            import traceback

            traceback.print_exc()

    # Output results
    if all_results:
        print_summary_table(all_results)
        analyze_scalability(all_results)

        # Save results to file
        output_file = "benchmark_results.txt"
        with open(output_file, "w", encoding="utf-8") as f:
            import json

            json.dump(all_results, f, indent=2, default=str)
        print(f"\nDetailed results saved to: {output_file}")
    else:
        print("No tests completed successfully")
        sys.exit(1)


if __name__ == "__main__":
    main()
