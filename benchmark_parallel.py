#!/usr/bin/env python3
"""
并行化性能对比脚本

测试KDTree并行查询组件的时/空性能：
1. 比较不同线程数下的运行时间
2. 验证并行结果与串行结果一致性
3. 分析加速比和扩展性
4. 测量内存使用变化

命题验证：每个氨基酸与其最近水分子的距离计算相互独立，
因此在构建好KDTree后可以并行计算。
"""

import sys
import time
import argparse
from pathlib import Path
from typing import List, Dict, Any
import numpy as np

# 添加项目路径以便导入模块
sys.path.insert(0, str(Path(__file__).parent))

from core.data_models import AnalysisConfig, MethodType
from io_utils.pdb_loader import PDBLoader
from algorithms.method_factory import MethodFactory


def benchmark_method(
    wet_pdb: str,
    dry_pdb: str,
    method_type: MethodType,
    num_processes_list: List[int],
    num_runs: int = 3,
) -> Dict[str, Any]:
    """
    基准测试方法

    Args:
        wet_pdb: 含水PDB文件路径
        dry_pdb: 无水PDB文件路径（用于FreeSASA对比）
        method_type: 方法类型（质心法或原子级方法）
        num_processes_list: 要测试的并行进程数列表
        num_runs: 每个配置的运行次数（取平均值）

    Returns:
        包含性能结果的字典
    """
    print(f"\n{'='*60}")
    print(f"基准测试: {method_type.value} 方法")
    print(f"测试文件: {Path(wet_pdb).name}")
    print(f"进程数列表: {num_processes_list}")
    print(f"运行次数: {num_runs}")
    print(f"{'='*60}")

    # 加载PDB文件
    loader = PDBLoader(quiet=True)
    residues, waters, structure = loader.load(wet_pdb)

    print(f"系统规模: {len(residues)} 个残基, {len(waters.coords) if waters.coords is not None else 0} 个水分子")

    results = {
        "method": method_type.value,
        "num_residues": len(residues),
        "num_waters": len(waters.coords) if waters.coords is not None else 0,
        "process_times": {},  # 进程数 -> 时间列表
        "speedups": {},       # 进程数 -> 加速比
        "memory_deltas": {},  # 进程数 -> 内存增量（近似）
        "result_consistency": True,  # 结果一致性标志
    }

    # 首先运行串行版本作为基准和正确性参考
    print(f"\n1. 串行基准 (num_processes=1)")
    serial_config = AnalysisConfig(num_processes=1)
    serial_method = MethodFactory.create_method(method_type, serial_config)

    serial_times = []
    serial_results_list = []

    for run in range(num_runs):
        start_time = time.perf_counter()
        serial_results = serial_method.analyze(residues, waters, structure)
        end_time = time.perf_counter()
        serial_times.append(end_time - start_time)
        serial_results_list.append(serial_results)

        if run == 0:
            # 记录串行结果用于比较
            serial_accessible = sum(1 for r in serial_results if r.accessible)
            print(f"   运行 {run+1}: {serial_times[-1]:.3f}s, 可及残基: {serial_accessible}/{len(serial_results)}")

    avg_serial_time = np.mean(serial_times)
    std_serial_time = np.std(serial_times)
    results["serial_time"] = avg_serial_time
    results["serial_time_std"] = std_serial_time
    print(f"   平均时间: {avg_serial_time:.3f}s (±{std_serial_time:.3f}s)")

    # 测试不同并行进程数
    for num_proc in num_processes_list:
        print(f"\n2. 并行测试 (num_processes={num_proc})")

        # 配置并行版本
        parallel_config = AnalysisConfig(num_processes=num_proc)
        parallel_method = MethodFactory.create_method(method_type, parallel_config)

        parallel_times = []
        all_consistent = True

        for run in range(num_runs):
            start_time = time.perf_counter()
            parallel_results = parallel_method.analyze(residues, waters, structure)
            end_time = time.perf_counter()
            parallel_times.append(end_time - start_time)

            # 验证结果一致性（与第一次串行运行比较）
            if run == 0:
                parallel_accessible = sum(1 for r in parallel_results if r.accessible)
                serial_accessible = sum(1 for r in serial_results_list[0] if r.accessible)

                # 详细比较（可选）
                if len(parallel_results) == len(serial_results_list[0]):
                    mismatches = 0
                    for p_res, s_res in zip(parallel_results, serial_results_list[0]):
                        if p_res.accessible != s_res.accessible:
                            mismatches += 1
                            if mismatches <= 3:  # 只打印前几个不匹配
                                print(f"     警告: 残基 {p_res.residue.chain}:{p_res.residue.resnum} "
                                      f"可及性不匹配 (并行: {p_res.accessible}, 串行: {s_res.accessible})")

                    if mismatches > 0:
                        print(f"     总不匹配数: {mismatches}/{len(parallel_results)}")
                        all_consistent = False
                    else:
                        print(f"     结果一致性: ✓ 全部匹配")
                else:
                    print(f"     错误: 结果长度不同 (并行: {len(parallel_results)}, 串行: {len(serial_results_list[0])})")
                    all_consistent = False

                print(f"   运行 {run+1}: {parallel_times[-1]:.3f}s, 可及残基: {parallel_accessible}/{len(parallel_results)}")
            else:
                print(f"   运行 {run+1}: {parallel_times[-1]:.3f}s")

        avg_parallel_time = np.mean(parallel_times)
        std_parallel_time = np.std(parallel_times)
        speedup = avg_serial_time / avg_parallel_time if avg_parallel_time > 0 else 0

        results["process_times"][num_proc] = {
            "mean": avg_parallel_time,
            "std": std_parallel_time,
            "speedup": speedup,
            "efficiency": speedup / num_proc if num_proc > 0 else 0,
        }

        if not all_consistent:
            results["result_consistency"] = False

        print(f"   平均时间: {avg_parallel_time:.3f}s (±{std_parallel_time:.3f}s)")
        print(f"   加速比: {speedup:.2f}x (效率: {speedup/num_proc*100:.1f}%)")

    return results


def print_summary_table(results_list: List[Dict[str, Any]]) -> None:
    """打印性能汇总表格"""
    print(f"\n{'='*80}")
    print("性能对比汇总")
    print(f"{'='*80}")

    for result in results_list:
        method = result["method"]
        print(f"\n方法: {method.upper()}")
        print(f"系统: {result['num_residues']}残基, {result['num_waters']}水分子")
        print(f"结果一致性: {'✓' if result['result_consistency'] else '✗'}")

        print(f"\n{'进程数':<8} {'平均时间(s)':<12} {'标准差':<10} {'加速比':<10} {'效率(%)':<10}")
        print(f"{'-'*50}")

        # 串行结果
        print(f"{'1 (串行)':<8} {result['serial_time']:<12.3f} {result['serial_time_std']:<10.3f} {'1.00':<10} {'100.0':<10}")

        # 并行结果
        for num_proc, metrics in sorted(result["process_times"].items()):
            print(f"{num_proc:<8} {metrics['mean']:<12.3f} {metrics['std']:<10.3f} "
                  f"{metrics['speedup']:<10.2f} {metrics['efficiency']*100:<10.1f}")


def analyze_scalability(results_list: List[Dict[str, Any]]) -> None:
    """分析扩展性特征"""
    print(f"\n{'='*80}")
    print("扩展性分析")
    print(f"{'='*80}")

    for result in results_list:
        method = result["method"]
        print(f"\n{method.upper()} 方法:")

        if len(result["process_times"]) < 2:
            print("  数据不足进行扩展性分析")
            continue

        # 计算理想加速比（Amdahl定律）
        # 假设并行部分比例p，串行部分比例1-p
        # 这里简单估计：使用最大进程数的加速比
        max_proc = max(result["process_times"].keys())
        max_speedup = result["process_times"][max_proc]["speedup"]

        if max_speedup > 1:
            # 根据Amdahl定律反推并行比例
            # S = 1 / ((1-p) + p/N)
            # 其中S为加速比，N为进程数
            # 解出 p = (1/S - 1) / (1/N - 1)
            p = (1/max_speedup - 1) / (1/max_proc - 1)
            serial_fraction = 1 - p

            print(f"  最大加速比: {max_speedup:.2f}x (在{max_proc}个线程下)")
            print(f"  估计并行比例: {p*100:.1f}%")
            print(f"  串行瓶颈: {serial_fraction*100:.1f}%")

            # 判断扩展性类型
            efficiency = result["process_times"][max_proc]["efficiency"]
            if efficiency > 0.7:
                print(f"  扩展性: 优秀 (效率 > 70%)")
            elif efficiency > 0.5:
                print(f"  扩展性: 良好 (效率 > 50%)")
            elif efficiency > 0.3:
                print(f"  扩展性: 一般 (效率 > 30%)")
            else:
                print(f"  扩展性: 较差 (效率 ≤ 30%)")
        else:
            print(f"  无加速效果")


def main():
    parser = argparse.ArgumentParser(
        description="并行化性能对比测试脚本"
    )
    parser.add_argument(
        "--wet-pdb",
        default="./pdb/SUMO1_water.pdb",
        help="含水PDB文件路径 (默认: ./pdb/SUMO1_water.pdb)"
    )
    parser.add_argument(
        "--dry-pdb",
        default="./pdb/SUMO1.pdb",
        help="无水PDB文件路径 (默认: ./pdb/SUMO1.pdb)"
    )
    parser.add_argument(
        "--processes",
        type=int,
        nargs="+",
        default=[1, 2, 4, 8],
        help="要测试的并行进程数列表 (默认: 1 2 4 8)"
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=3,
        help="每个配置的运行次数 (默认: 3)"
    )
    parser.add_argument(
        "--methods",
        choices=["centroid", "peratom", "both"],
        default="both",
        help="要测试的方法 (默认: both)"
    )

    args = parser.parse_args()

    # 验证文件存在
    if not Path(args.wet_pdb).exists():
        print(f"错误: 文件不存在: {args.wet_pdb}")
        sys.exit(1)

    if not Path(args.dry_pdb).exists():
        print(f"警告: 文件不存在: {args.dry_pdb} (FreeSASA对比将跳过)")

    # 确定要测试的方法
    if args.methods == "both":
        method_types = [MethodType.CENTROID, MethodType.PERATOM]
    elif args.methods == "centroid":
        method_types = [MethodType.CENTROID]
    else:
        method_types = [MethodType.PERATOM]

    # 运行基准测试
    all_results = []

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
            print(f"测试方法 {method_type.value} 时出错: {e}")
            import traceback
            traceback.print_exc()

    # 输出结果
    if all_results:
        print_summary_table(all_results)
        analyze_scalability(all_results)

        # 保存结果到文件
        output_file = "benchmark_results.txt"
        with open(output_file, "w", encoding="utf-8") as f:
            import json
            json.dump(all_results, f, indent=2, default=str)
        print(f"\n详细结果已保存到: {output_file}")
    else:
        print("没有成功完成任何测试")
        sys.exit(1)


if __name__ == "__main__":
    main()