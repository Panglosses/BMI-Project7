"""
验证工具
"""

from pathlib import Path

from core.data_models import AnalysisConfig


def validate_pdb_file(filepath: str | Path) -> bool:
    """
    验证PDB文件

    Args:
        filepath: PDB文件路径

    Returns:
        bool: 是否有效
    """
    path = Path(filepath)

    # 检查文件是否存在
    if not path.exists():
        raise FileNotFoundError(f"PDB文件不存在: {filepath}")

    # 检查文件大小
    if path.stat().st_size == 0:
        raise ValueError(f"PDB文件为空: {filepath}")

    # 检查文件扩展名
    if path.suffix.lower() not in [".pdb", ".ent"]:
        raise ValueError(f"文件扩展名不是PDB格式: {filepath}")

    # 简单检查文件内容
    try:
        with open(path, "r") as f:
            first_line = f.readline().strip()
            # PDB文件通常以"HEADER"、"ATOM"或"HETATM"开头
            if not any(
                first_line.startswith(prefix)
                for prefix in ["HEADER", "ATOM", "HETATM", "MODEL"]
            ):
                # 有些PDB文件可能没有HEADER，但第一行应该是有效的PDB记录
                if len(first_line) > 0 and first_line[0] not in [" ", "\t"]:
                    # 检查是否是有效的PDB记录类型
                    record_type = first_line[0:6].strip()
                    if record_type not in [
                        "HEADER",
                        "ATOM",
                        "HETATM",
                        "MODEL",
                        "ENDMDL",
                        "TER",
                        "END",
                    ]:
                        raise ValueError(
                            f"PDB文件格式无效，第一行: {first_line[:50]}..."
                        )
    except UnicodeDecodeError:
        raise ValueError(f"PDB文件不是文本格式: {filepath}")

    return True


def validate_config(config: AnalysisConfig) -> bool:
    """
    验证配置

    Args:
        config: 分析配置

    Returns:
        bool: 是否有效
    """
    config.validate()
    return True


def validate_output_dir(directory: str | Path) -> Path:
    """
    验证输出目录

    Args:
        directory: 目录路径

    Returns:
        Path: 有效的目录路径
    """
    path = Path(directory)

    # 如果目录不存在，尝试创建
    if not path.exists():
        try:
            path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise ValueError(f"无法创建输出目录 {directory}: {e}")

    # 检查是否可写
    if not path.is_dir():
        raise ValueError(f"输出路径不是目录: {directory}")

    # 简单检查是否可写
    test_file = path / ".write_test"
    try:
        test_file.touch()
        test_file.unlink()
    except Exception as e:
        raise ValueError(f"输出目录不可写 {directory}: {e}")

    return path
