""" Overall command with variable parameters
python cli/main.py \
  --wet-pdb "path/to/wet_structure.pdb" \          # Hydrated PDB file (with water molecules)
  --dry-pdb "path/to/dry_structure.pdb" \          # Raw PDB file (without water, for FreeSASA)
  --method "peratom" \                             # Analysis method: "centroid" or "peratom"
  --threshold 3.5 \                                # Accessibility threshold (Å): distance < threshold = accessible
  --margin 2.0 \                                   # Centroid filter margin (Å): centroid distance > (threshold + margin) = filtered out
  --R 5.0 \                                        # Radius for water counting (Å): count waters within the radius
  --fraction-threshold 0.20 \                      # Fraction threshold for peratom method (0-1): higher than this ratio of atoms within threshold is considered accessible
  --min-hits 1 \                                   # Minimum atom hits for peratom method: at least this count of atoms within threshold is considered accessible
  --small-residue-size 5 \                         # Small residue threshold: residues with ≤ this count atoms apply different rules
  --chunk 5000 \                                   # Chunk size for distance calculations: larger ---> more memory and faster
  --nproc 4 \                                      # Number of parallel processes: use more cores for faster computation
  --output-dir "./results" \                       # Output directory for results files
  --verbose \                                      # Show detailed progress and statistics
  --no-comparison                                  # Skip FreeSASA comparison (dry-pdb still required but not used)
"""




"""
主命令行接口
"""

import argparse
import sys
import os
from pathlib import Path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)

# 添加项目根目录到Python路径
sys.path.insert(0, project_root)

# 现在导入其他模块

from core.data_models import AnalysisConfig, MethodType
from io_utils import PDBLoader, CSVWriter, ResultFormatter
from algorithms import MethodFactory, FreeSASAWrapper


def parse_args(args=None):
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="溶剂可及性分析工具 - 基于水分子接近度"
    )

    # 输入文件
    parser.add_argument(
        "--wet-pdb", required=True, help="水合PDB文件（用于自定义方法分析）"
    )
    parser.add_argument(
        "--dry-pdb", required=True, help="无水PDB文件（用于FreeSASA分析）"
    )

    # 方法选择
    parser.add_argument(
        "--method",
        choices=["centroid", "peratom"],
        default="peratom",
        help="分析方法：centroid（质心法）或 peratom（原子级方法）",
    )

    # 距离参数
    parser.add_argument(
        "--threshold", type=float, default=3.5, help="可及性判断阈值（Å），默认：3.5"
    )
    parser.add_argument(
        "--margin", type=float, default=2.0, help="质心法的额外裕度（Å），默认：2.0"
    )
    parser.add_argument(
        "--R", type=float, default=5.0, help="统计水分子的半径（Å），默认：5.0"
    )

    # 原子级方法参数
    parser.add_argument(
        "--fraction-threshold",
        type=float,
        default=0.20,
        help="原子可及比例阈值（0-1），默认：0.20",
    )
    parser.add_argument(
        "--min-hits", type=int, default=1, help="最小命中原子数，默认：1"
    )
    parser.add_argument(
        "--small-residue-size", type=int, default=5, help="小残基的原子数阈值，默认：5"
    )

    # 计算参数
    parser.add_argument(
        "--chunk", type=int, default=5000, help="分块计算大小，默认：5000"
    )
    parser.add_argument("--nproc", type=int, default=1, help="并行进程数，默认：1")

    # 输出控制
    parser.add_argument(
        "--output-dir", default="./output", help="输出目录，默认：./output"
    )
    parser.add_argument("--verbose", action="store_true", help="显示详细输出")
    parser.add_argument(
        "--no-comparison", action="store_true", help="不进行FreeSASA对比"
    )

    return parser.parse_args(args)


def create_config(args) -> AnalysisConfig:
    """从命令行参数创建配置"""
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
    """运行自定义方法分析"""
    if args.verbose:
        print(f"加载PDB文件: {args.wet_pdb}")

    # 加载PDB
    loader = PDBLoader(quiet=not args.verbose)
    residues, waters, structure = loader.load(args.wet_pdb)

    if args.verbose:
        print(f"  残基数: {len(residues)}")
        print(f"  水分子数: {waters.count}")

    # 创建分析方法
    method = MethodFactory.create_method(args.method, config)

    # 执行分析
    if args.verbose:
        print(f"执行{args.method}分析...")

    results = method.analyze(residues, waters, structure)

    if args.verbose:
        summary = ResultFormatter.format_summary(results)
        print(summary)

    return residues, results


def run_freesasa_analysis(args, config: AnalysisConfig):
    """运行FreeSASA分析"""
    if args.verbose:
        print(f"运行FreeSASA分析: {args.dry_pdb}")

    wrapper = FreeSASAWrapper(config)
    sasa_results = wrapper.compute_residue_sasa(args.dry_pdb)

    if args.verbose:
        accessible = sum(1 for r in sasa_results if r["Accessible"] == "Yes")
        total = len(sasa_results)
        print(f"  FreeSASA结果: {accessible}/{total} 可及")

    return sasa_results


def save_results(args, custom_results, sasa_results=None):
    """保存结果文件"""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 保存自定义方法结果
    prefix = Path(args.wet_pdb).stem
    custom_file = output_dir / f"{prefix}_{args.method}.csv"
    CSVWriter.write_results(str(custom_file), custom_results)

    if args.verbose:
        print(f"保存自定义方法结果: {custom_file}")

    # 保存FreeSASA结果
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
            print(f"保存FreeSASA结果: {sasa_file}")

    return custom_file


def compare_results(custom_results, sasa_results):
    """对比结果 - 使用旧代码逻辑"""
    
    def normalize_chain(c):
        """标准化链标识（与旧代码完全一致）"""
        if isinstance(c, str):
            c = c.strip()
        else:
            c = str(c).strip()
        return c if c else "A"  # 空链转为"A"
    
    # 创建SASA映射（与旧代码一致）
    sasa_map = {}
    for item in sasa_results:
        chain = normalize_chain(item.get("chain", ""))
        resnum = str(item.get("resnum", ""))
        accessible = item.get("Accessible", "No")
        sasa_map[(chain, resnum)] = accessible == "Yes"
    
    # 计算匹配
    match_count = 0
    total = 0
    
    for result in custom_results:
        # 使用相同的标准化函数
        chain = normalize_chain(result.residue.chain)
        resnum = str(result.residue.resnum)
        key = (chain, resnum)
        
        sasa_accessible = sasa_map.get(key, False)
        
        if result.accessible == sasa_accessible:
            match_count += 1
        total += 1
    
    return match_count / total if total > 0 else 0.0


def main(args=None):
    """主函数"""
    if args is None:
        args = parse_args()

    try:
        # 创建配置
        config = create_config(args)

        # 运行自定义方法分析
        residues, custom_results = run_custom_analysis(args, config)

        # 运行FreeSASA分析
        if not args.no_comparison:
            sasa_results = run_freesasa_analysis(args, config)

            # 对比结果
            match_ratio = compare_results(custom_results, sasa_results)

            # 保存对比结果
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
                print(f"\n=== 对比完成 ===")
                print(f"匹配比例: {match_ratio:.4f}")
                print(f"对比结果: {comparison_file}")

        # 保存结果
        custom_file = save_results(
            args, custom_results, sasa_results if not args.no_comparison else None
        )

        if args.verbose:
            print(f"\n分析完成！")
            print(f"结果文件: {custom_file}")
            if not args.no_comparison:
                print(f"对比文件: {Path(args.output_dir) / 'comparison.csv'}")

        return 0

    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
