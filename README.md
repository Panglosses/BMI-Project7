# 溶剂可及性分析工具包

基于水分子接近度的蛋白质残基溶剂可及性高效计算工具。

## 项目简介

蛋白质折叠由疏水效应驱动，促使疏水残基埋藏于蛋白核心，亲水残基暴露于水环境。传统方法基于几何表面积计算（如FreeSASA），本方法直接评估残基与显式水分子的接近程度，提供更物理真实的溶剂可及性度量。

**核心特性**：
- 双算法支持：质心法（残基质心距离）和原子级方法（原子接触比例）
- 高性能计算：分块计算、KDTree空间索引、向量化操作
- 模块化架构：策略模式设计，易于扩展新算法
- 并行化支持：基于Python 3.14自由线程特性的可拔插并行组件

## 快速开始

### 安装依赖

```bash
# 给脚本添加执行权限（Linux/macOS）
chmod +x setup.sh

# 运行安装脚本
./setup.sh

# 该脚本可自动创建虚拟环境bmi并安装所有依赖
```

### 添加数据

将待处理的pdb文件（如`SUMO1_water.pdb`）放进`./pdb/`内

<!-- ### 命令行使用

```bash
# 原子级方法（默认）
python -m solvent_analysis \
  --wet-pdb ./pdb/SUMO1_water.pdb \
  --dry-pdb ./pdb/SUMO1.pdb \
  --method peratom \
  --verbose
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
print(f"可及残基: {accessible}/{len(results)}")
``` -->

## 详细教程

完整代码讲解、算法原理、性能优化和并行化设计请查看：

**[tutorial.ipynb](tutorial.ipynb)** - Jupyter notebook详细教程，包含：
- 逐步代码解析与算法对比
- 性能基准测试与并行化分析
- 交互式运行和实验

## 项目结构

```
solvent_analysis/          # 主包
├── core/                  # 核心接口与数据模型
├── io_utils/              # 输入输出模块
├── algorithms/            # 算法实现
├── utils/                 # 工具模块
└── cli/                   # 命令行接口
```

**核心模块**：
- `core/data_models.py` - 数据类定义（ResidueInfo, WaterInfo等）
- `core/distance_calculator.py` - 距离计算抽象接口
- `algorithms/centroid_method.py` - 质心法实现
- `algorithms/peratom_method.py` - 原子级方法实现
- `io_utils/pdb_loader.py` - PDB文件加载（BioPython集成）

## 测试与验证

```bash
# 运行单元测试
python test_refactored.py

# 性能基准测试
python benchmark_parallel.py --processes 1 2 4
```

## 性能优化

- **分块计算**：控制内存峰值（`chunk_size`参数）
- **空间索引**：KDTree加速最近邻搜索
- **并行计算**：配置`num_processes`参数利用多核
- **向量化操作**：NumPy广播优化距离计算

## 扩展开发

项目采用策略模式和工厂模式设计，支持轻松扩展：
1. 实现`DistanceCalculator`抽象基类添加新距离算法
2. 实现`AccessibilityEvaluator`抽象基类添加新评估规则
3. 在`MethodFactory`中注册新方法
