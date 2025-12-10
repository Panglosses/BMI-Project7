# 基本使用示例

import sys
from pathlib import Path

# 添加父目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent))


def example_python_api():
    """Python API使用示例"""
    print("=== Python API 使用示例 ===")

    from solvent_analysis import ResidueInfo, WaterInfo, AnalysisConfig, MethodType
    from io_utils.pdb_loader import PDBLoader
    from algorithms import MethodFactory
    from io_utils.csv_writer import CSVWriter

    # 1. 创建配置
    config = AnalysisConfig(
        threshold=3.5,  # 可及性阈值（Å）
        radius=5.0,  # 统计半径（Å）
        fraction_threshold=0.20,  # 原子可及比例
        min_hits=1,  # 最小命中原子数
        chunk_size=5000,  # 分块计算大小
    )

    print(f"配置创建成功: threshold={config.threshold}, radius={config.radius}")

    # 2. 加载PDB文件
    pdb_path = "../pdb/SUMO1_water.pdb"
    if not Path(pdb_path).exists():
        print(f"⚠ PDB文件不存在: {pdb_path}")
        print("  请先运行测试或使用其他PDB文件")
        return

    loader = PDBLoader(quiet=True)
    residues, waters, structure = loader.load(pdb_path)

    print(f"PDB加载成功:")
    print(f"  残基数: {len(residues)}")
    print(f"  水分子数: {waters.count}")
    print(f"  示例残基: {residues[0] if residues else '无'}")

    # 3. 创建分析方法
    method = MethodFactory.create_method(MethodType.PERATOM, config)
    print(f"分析方法创建成功: {method.get_method_type()}")

    # 4. 执行分析
    print("执行分析...")
    results = method.analyze(residues, waters, structure)

    # 5. 处理结果
    print(f"分析完成，结果数: {len(results)}")

    # 统计可及残基
    accessible = sum(1 for r in results if r.accessible)
    ratio = accessible / len(results) if results else 0

    print(f"可及残基: {accessible}/{len(results)} ({ratio:.1%})")

    # 显示前几个结果
    print("\n前5个结果:")
    for i, result in enumerate(results[:5]):
        status = "可及" if result.accessible else "不可及"
        print(
            f"  {result.residue.chain}{result.residue.resnum} {result.residue.resname}: "
            f"距离={result.min_distance:.2f}Å, 水分子={result.water_count}, {status}"
        )

    # 6. 保存结果
    output_file = "../output/example_results.csv"
    CSVWriter.write_results(output_file, results)
    print(f"\n结果已保存到: {output_file}")

    return results


def example_configuration():
    """配置使用示例"""
    print("\n=== 配置使用示例 ===")

    from solvent_analysis import AnalysisConfig

    # 默认配置
    default_config = AnalysisConfig()
    print(
        f"默认配置: threshold={default_config.threshold}, radius={default_config.radius}"
    )

    # 自定义配置
    custom_config = AnalysisConfig(
        threshold=4.0,  # 更宽松的阈值
        radius=6.0,  # 更大的统计半径
        fraction_threshold=0.15,  # 更宽松的比例阈值
        min_hits=2,  # 至少2个原子命中
        small_residues=("GLY", "ALA", "SER"),  # 自定义小残基集合
    )

    print(
        f"自定义配置: threshold={custom_config.threshold}, min_hits={custom_config.min_hits}"
    )

    # 验证配置
    try:
        custom_config.validate()
        print("配置验证成功")
    except ValueError as e:
        print(f"配置验证失败: {e}")


def example_advanced_usage():
    """高级使用示例"""
    print("\n=== 高级使用示例 ===")

    from algorithms import FreeSASAWrapper
    from utils.progress import ProgressBar
    from utils.logger import setup_logger
    import logging

    # 1. 设置日志
    logger = setup_logger(level=logging.INFO, console=True)
    logger.info("开始高级示例")

    # 2. 使用FreeSASA
    pdb_path = "../pdb/SUMO1.pdb"
    if Path(pdb_path).exists():
        wrapper = FreeSASAWrapper()
        sasa_results = wrapper.compute_residue_sasa(pdb_path)
        logger.info(f"FreeSASA计算完成: {len(sasa_results)}个残基")
    else:
        logger.warning(f"PDB文件不存在: {pdb_path}")

    # 3. 使用进度条
    print("\n进度条示例:")
    items = list(range(100))
    for item in ProgressBar.iterate(items, prefix="处理中", suffix="完成"):
        # 模拟处理
        pass

    print("高级示例完成")


def main():
    """运行所有示例"""
    print("溶剂可及性分析工具包 - 使用示例")
    print("=" * 50)

    try:
        # 运行示例
        results = example_python_api()
        example_configuration()
        example_advanced_usage()

        print("\n" + "=" * 50)
        print("所有示例运行完成！")
        print("\n更多用法请参考:")
        print("1. solvent_analysis/README.md - 模块说明")
        print("2. MIGRATION_GUIDE.md - 迁移指南")
        print("3. test_refactored.py - 测试脚本")

    except ImportError as e:
        print(f"导入错误: {e}")
        print("请确保在项目根目录运行，且虚拟环境已激活")
    except Exception as e:
        print(f"示例运行错误: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
