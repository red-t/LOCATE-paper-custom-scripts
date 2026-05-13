# TE Insertion Simulation Pipeline

用于生成模拟 TE insertion 数据的通用 pipeline，支持多物种。

## 目录结构

```
simulation/
├── config_template.yaml          # 配置文件模板
├── species_configs/              # 物种特定配置
│   ├── dm6.yaml                  # Drosophila melanogaster 配置示例
│   └── rice.yaml                 # 水稻配置模板
├── scripts/                      # SLURM 脚本
│   ├── simulation_with_slurm.sh
│   ├── merge_with_slurm.sh
│   ├── minimap2_with_slurm.sh
│   ├── label_with_slurm.sh
│   └── classify_germline_somatic.sh
├── generate_filelists.py         # filelist 生成脚本
├── run_simulation.py             # 主控脚本
├── define_population_genome.py   # 定义群体基因组
├── build-population-genome.py    # 构建群体基因组
├── generate_TGS.py               # 生成 TGS reads
├── generate_NGS.py               # 生成 NGS reads
├── define_pgdf_utils.pyx         # Cython 模块
└── label/                        # 标注相关脚本
```

## 环境配置

### 1. 安装依赖

```shell
# 创建 simulation 环境
cd LOCATE-paper-custom-scripts/simulation
mamba create -n simulation
mamba activate simulation

# 安装依赖
mamba install VISOR bcftools samtools scipy pysam
pip install cython pyyaml

# -------------- #

# 创建 TEMP3 环境 (用于 labeling)
cd LOCATE-paper-custom-scripts/TEMP3
mamba create -n TEMP3 python=3.10
mamba activate TEMP3
mamba install -c conda-forge autogluon=1.0.0
mamba install samtools=1.17
mamba install minimap2=2.26
mamba install pysam
pip install cython==3.0.6
python setup.py build_ext -i
```

### 2. 编译 Cython 模块

```shell
cd simulation/
python setup.py build_ext -i && rm -r build && rm *c
```

## 使用方法

### 1. 准备参考文件

每个物种需要准备以下参考文件：

| 文件 | 说明 |要求 |
|-----|------|-----|
| template.fa | 参考基因组模板 | 需要 samtools faidx 索引 |
| transposon.fa | TE consensus 序列 | 需要 samtools faidx 索引 |
| repeat.bed | RepeatMasker 注释 | BED 格式 |
| gap.bed | Gap 区域注释 | BED 格式 |

### 2. 创建配置文件

复制配置模板并修改：

```shell
cp config_template.yaml species_configs/<species>.yaml
# 编辑配置文件，填入参考文件路径和参数
```

### 3. 运行 Pipeline

分步执行模式：

```shell
# 查看所有步骤
python run_simulation.py --config species_configs/dm6.yaml

# Step 0: 定义群体基因组 (本地运行)
python run_simulation.py --config species_configs/dm6.yaml --step 0

# Step 1: 生成 filelists
python run_simulation.py --config species_configs/dm6.yaml --step 1

# Step 2-5: SLURM 任务 (依次提交)
python run_simulation.py --config species_configs/dm6.yaml --step 2
python run_simulation.py --config species_configs/dm6.yaml --step 3
python run_simulation.py --config species_configs/dm6.yaml --step 4
python run_simulation.py --config species_configs/dm6.yaml --step 5

# Step 6: 分类 germline/somatic (本地运行)
python run_simulation.py --config species_configs/dm6.yaml --step 6
```

Dry-run 模式（只显示命令，不执行）：

```shell
python run_simulation.py --config species_configs/dm6.yaml --step 0 --dry-run
```

## Pipeline 步骤说明

| Step | 名称 | 说明 | 运行方式 |
|------|------|------|---------|
| 0 | define_pgdf | 定义群体基因组，生成 pgd 文件 | 本地 |
| 1 | generate_filelists | 生成 SLURM 任务所需的 filelists | 本地 |
| 2 | simulation | 构建 population genome，生成 TGS reads | SLURM |
| 3 | merge | 合并各 chromosome 的 reads | SLURM |
| 4 | minimap2 | 将 reads 比对到参考基因组 | SLURM |
| 5 | label | 检测并标注 TP/FP | SLURM |
| 6 | classify | 分类 germline/somatic insertion | 本地 |

## 输出文件

每个数据集的主要输出：

```
<output_dir>/<dataset_id>_<protocol>/
├── TGS.fasta                     # 模拟的 TGS reads
├── <chrom>/<chrom>.ins.summary   # Insertion summary
├── <chrom>/<chrom>.ins.sequence  # Insertion sequence
├── map-hifi/TGS.bam              # minimap2 比对结果
├── map-pb/TGS.bam
├── map-ont/TGS.bam
└── result_no_secondary/
    ├── TP_clt_G.txt              # Germline TP
    ├── TP_clt_S.txt              # Somatic TP
    ├── TP_clt.txt                # All TP
    └── FP_clt.txt                # FP
```

## 物种支持

### 已有物种特定逻辑

- **human**: 特定的 TE 选择权重 (Alu 83%, LINE1 12%, SVA 4%, HERVK 1%)
- **fly**: 特定的 TE 选择权重 (P_element 58% 等)

### 新物种

新物种使用通用逻辑：
- TE 选择：从 transposon.fa 中均匀随机选择
- Truncation：使用 `trunc_prob` 参数控制
- Frequency：使用均匀分布 (0.1-1)

## 参数说明

### simulation 参数

| 参数 | 说明 | 默认值 |
|-----|------|-------|
| population_size | 群体基因组大小 | 100 |
| sub_population_size | 子群体大小 | 5 |
| germline_count | Germline insertion 数量 | 1000 |
| avg_somatic_count | 平均 somatic insertion 数量 | 50 |
| min_distance | TE insertion 最小距离 | 9000 |
| trunc_prob | Truncation 概率 | 0.1 |
| nest_prob | Nested insertion 概率 | 0.1 |

### sequencing 参数

| 参数 | 说明 | 默认值 |
|-----|------|-------|
| depths |测序深度列表 | [10, 20, 30, 40, 50] |
| protocols | 测序协议 | [ccs, clr, ont] |

## 常见问题

### Q: 如何为新物种准备参考文件？

1. **template.fa**: 从参考基因组中提取主要染色体，移除 alternate contigs
2. **transposon.fa**: 从 RepeatMasker 结果中提取 TE consensus 序列
3. **repeat.bed**: RepeatMasker 输出的 BED 格式
4. **gap.bed**: 从参考基因组中识别 gap 区域

### Q: SLURM 任务失败怎么办？

检查日志文件位于 `<output_dir>/logs/`，常见问题：
- 内存不足：调整 `slurm.mem_*` 参数
- 时间不足：调整 `slurm.time` 参数
- conda 环境问题：检查 `conda_path` 和 `conda_env`

### Q: 如何只生成特定深度的数据集？

修改配置文件中的 `datasets` 部分，例如只生成 30X：

```yaml
datasets:
  30: 20
```