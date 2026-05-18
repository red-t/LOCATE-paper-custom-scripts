# LOCATE — TE Insertion Simulation Pipeline

用于生成模拟 TE（Transposable Element）insertion 数据的通用 pipeline，支持多物种、多测序协议，可生成训练数据和测试数据。

## 仓库结构

```
LOCATE-paper-custom-scripts/
├── simulation/                    # 模拟数据生成 pipeline
│   ├── run_simulation.py          # 主控脚本（7 步 pipeline）
│   ├── generate_filelists.py      # SLURM 任务 filelist 生成
│   ├── define_population_genome.py
│   ├── build-population-genome.py
│   ├── generate_TGS.py            # TGS reads 生成
│   ├── generate_NGS.py            # NGS reads 生成
│   ├── config_template.yaml       # 配置文件模板
│   ├── species_configs/           # 物种特定配置
│   │   ├── dm6.yaml               # D. melanogaster
│   │   ├── rice.yaml              # 水稻（模板）
│   │   └── test.yaml              # 测试用最小配置
│   ├── scripts/                   # SLURM 作业脚本
│   │   ├── simulation_with_slurm.sh
│   │   ├── merge_with_slurm.sh
│   │   ├── minimap2_with_slurm.sh
│   │   └── label_with_slurm.sh
│   ├── label/                     # 标注与分类脚本
│   │   ├── label_protocol.sh
│   │   ├── label_protocol_local.sh
│   │   ├── Filter_TP_and_FP.py
│   │   ├── rename.py
│   │   ├── merge_TP_and_FP.sh
│   │   ├── classify_germline_somatic.py
│   │   └── classify_germline_somatic.sh
│   ├── simulate_training_data/    # 训练数据生成脚本（GRCh38 / dm6）
│   └── simulate_test_data/        # 测试数据生成脚本（GRCh38 / dm6）
├── TEMP3/                         # TE insertion 检测工具
│   ├── TEMP3.py                   # 主入口
│   ├── TEMP3/                     # 核心模块（Cython + C）
│   └── setup.py
└── model_training/                # AutoGluon 分类模型训练
    ├── training.py                # 训练入口（支持 --drop-cols / --pos-label）
    └── pipeline/                  # ML 训练 pipeline（Step 7）
        ├── run_pipeline.py        # 独立编排脚本
        ├── config.py              # YAML 配置解析
        ├── defaults.py            # 列名常量与 mode 默认参数
        ├── utils.py               # 工具函数
        ├── evaluate.py            # 共享评估函数
        ├── stage1_preprocess.py   # 数据预处理（merge / split / blacklist / dedup）
        ├── stage2_shuffle.py      # 按比例采样 P/N 数据
        ├── stage3_train.py        # 遍历比例训练 AutoGluon 模型
        └── stage4_filter_eval.py  # 黑名单过滤 + 模型评估
```

## 环境配置

### 1. 创建 Simulation 环境

```shell
cd simulation/
mamba create -n simulation python=3.12
mamba activate simulation
mamba install VISOR bcftools samtools scipy pysam
pip install cython pyyaml
python setup.py build_ext -i && rm -r build && rm -f *.c
```

### 2. 创建 TEMP3 环境（用于 labeling）

```shell
cd TEMP3/
mamba create -n TEMP3 python=3.10
mamba activate TEMP3
mamba install -c conda-forge autogluon=1.0.0
mamba install samtools=1.17 minimap2=2.26 pysam setuptools=69.5.1
pip install cython==3.0.6
python setup.py build_ext -i && rm -r build && rm -f *.c
```

## 快速开始

### 1. 准备参考文件

每个物种需要以下参考文件：

| 文件 | 说明 | 要求 |
|------|------|------|
| `template.fa` | 参考基因组模板 | 需 samtools faidx 索引 |
| `transposon.fa` | TE consensus 序列 | 需 samtools faidx 索引 |
| `repeat.bed` | RepeatMasker 注释 | BED 格式 |
| `gap.bed` | Gap 区域注释 | BED 格式 |

### 2. 创建配置文件

```shell
cp config_template.yaml species_configs/<species>.yaml
# 编辑配置文件，填入参考文件路径和参数
```

### 3. 运行 Pipeline

```shell
# 查看所有步骤
python run_simulation.py --config species_configs/dm6.yaml

# Step 0: 定义群体基因组（本地运行）
python run_simulation.py --config species_configs/dm6.yaml --step 0

# Step 1: 生成 filelists
python run_simulation.py --config species_configs/dm6.yaml --step 1

# Step 2-5: 依次提交 SLURM 任务
python run_simulation.py --config species_configs/dm6.yaml --step 2
python run_simulation.py --config species_configs/dm6.yaml --step 3
python run_simulation.py --config species_configs/dm6.yaml --step 4
python run_simulation.py --config species_configs/dm6.yaml --step 5

# Step 6: 分类 germline/somatic（本地运行）
python run_simulation.py --config species_configs/dm6.yaml --step 6

# Step 7: ML 模型训练（本地运行）
python run_simulation.py --config species_configs/dm6.yaml --step 7

# 或跳过前置步骤，只执行 Step 7:
python run_simulation.py --config species_configs/dm6.yaml --start 7 --step 7

# Step 7 也可以独立运行（不依赖于 run_simulation.py）:
python model_training/pipeline/run_pipeline.py --config species_configs/dm6.yaml --mode all
```

Dry-run 模式（仅显示命令，不执行）：

```shell
python run_simulation.py --config species_configs/dm6.yaml --step 0 --dry-run
```

## Pipeline 步骤详解

| Step | 名称 | 说明 | 运行方式 | 核心脚本 |
|------|------|------|----------|----------|
| 0 | define_pgdf | 定义群体基因组，生成 PGD 文件 | 本地 | `define_population_genome.py` |
| 1 | generate_filelists | 生成 SLURM 任务所需的 filelists | 本地 | `generate_filelists.py` |
| 2 | simulation | 构建 population genome，生成 TGS reads | SLURM | `scripts/simulation_with_slurm.sh` |
| 3 | merge | 合并各 chromosome 的 reads | SLURM | `scripts/merge_with_slurm.sh` |
| 4 | minimap2 | 将 reads 比对到参考基因组 | SLURM | `scripts/minimap2_with_slurm.sh` |
| 5 | label | 检测 TE insertion 并标注 TP/FP | SLURM | `scripts/label_with_slurm.sh` → TEMP3 + `label/Filter_TP_and_FP.py` |
| 6 | classify | 分类 germline/somatic insertion | 本地 | `label/classify_germline_somatic.py` |
| 7 | ml_training | ML 模型训练（AutoGluon） | 本地 | `model_training/pipeline/run_pipeline.py` |

### Step 0 — define_pgdf

定义群体基因组（Population Genome Definition），为每个 contig 生成 `.pgd` 文件，记录所有 TE insertion 的位置、类型和频率信息。

### Step 1 — generate_filelists

根据配置生成四个 filelist 文件，供 SLURM 任务使用：
- `simulation_filelist` — Step 2 的输入
- `merge_filelist` — Step 3 的输入
- `minimap2_filelist` — Step 4 的输入
- `label_filelist` — Step 5 的输入

### Step 2 — simulation

构建群体基因组并生成 TGS reads。每个 SLURM 任务处理一个 contig/子群体组合，生成对应染色体的 reads。

### Step 3 — merge

将 per-contig 的 TGS reads 合并为单个 `TGS.fasta` 文件。

### Step 4 — minimap2

使用 minimap2 将 reads 比对到参考基因组，根据协议选择不同的 preset（map-hifi / map-pb / map-ont），输出排序后的 BAM 文件。

### Step 5 — label

运行 TEMP3 检测 TE insertion candidates，然后通过 `Filter_TP_and_FP.py` 与 ground truth 比对，分类为 TP（真阳性）和 FP（假阳性）。

### Step 6 — classify

根据 insertion 频率将 TP 分为 germline（高频）和 somatic（低频），输出 `TP_clt_G.txt` 和 `TP_clt_S.txt`。

### Step 7 — ml_training

使用 Step 5-6 输出的标注数据训练 AutoGluon 分类模型。通过 YAML 配置的 `ml_training` 段控制参数。

流程：
1. **preprocess** — 遍历 filelist 读取各样本的 `TP_clt_G.txt`/`TP_clt_S.txt`/`FP_clt.txt`，按 `cltType` 和 `teAlignedFrac` 过滤、插入 sample_id、合并
2. **shuffle** — 按配置的 P/N 比例（如 ORG、1V1、1V30）对正负样本采样
3. **train** — 对每个比例训练 AutoGluon TabularPredictor 模型
4. **evaluate** — 用 BlackList 过滤测试数据、评估模型、输出 summary

支持两种 mode：
- **germline**（高频）：`cltType==0`、`teAlignedFrac>=0.8`、正类来源仅 `TP_clt_G.txt`
- **somatic**（低频）：`cltType>0`、`teAlignedFrac>=1`、正类来源 `TP_clt_G.txt` + `TP_clt_S.txt`

输出目录结构：
```
{output_dir}/For_ML/{species}/
  germline/
    Train/         train_P/N_{ratio}.txt
    Test/          test_P/N_{protocol}_{ratio}.txt
    BlackList*.bed
    {species}_G_{ratio}/        # AutoGluon 模型目录
    {species}_G_summary_Dedup.txt
  somatic/
    Train/
    Test/
    BlackList*.bed
    {species}_S_1V30/
    {species}_S_summary_Dedup.txt
```

## 输出文件

每个数据集的主要输出结构：

```
<output_dir>/<dataset_id>_<protocol>/
├── TGS.fasta                         # 模拟的 TGS reads
├── <chrom>/
│   ├── <chrom>.ins.summary           # Insertion 摘要
│   └── <chrom>.ins.sequence          # Insertion 序列
├── map-hifi/TGS.bam                  # HiFi 比对结果
├── map-pb/TGS.bam                    # PacBio CLR 比对结果
├── map-ont/TGS.bam                   # ONT 比对结果
└── result_no_secondary/
    ├── TP_clt.txt                    # 全部真阳性
    ├── TP_clt_G.txt                  # Germline TP
    ├── TP_clt_S.txt                  # Somatic TP
    └── FP_clt.txt                    # 假阳性
```

### 输出文件字段说明

```
Column  Value               Description

1       chrom               chromosome
2       refStart            cluster start on reference sequence (0-based, included)
3       refEnd              cluster end on reference sequence (0-based, not-included)
4       cltID               cluster ID
5       numSeg              number of segments in the cluster (normalized by bg depth)
6       strand              cluster orientation
7       startIndex          start index in segments array (0-based, include)
8       endIndex            end index in segments array (0-based, not-include)
9       numSeg              number of segments in the cluster (normalized by bg depth)
10      directionFlag       bitwise flag representing cluster direction
                                1: forward
                                2: reverse
                                other: unkown
11      cltType             cluster type
                                0: germline (multiple support reads)
                                1: somatic (1 support read & 1 alignment)
                                2: somatic (1 support read & 2 alignments)
12      locationType        bitwise flag representing cluster location
                                1: inside normal region
                                2: at repeat/gap boundary
                                4: inside repeat/gap
13      numSegType          number of different segment types
14      entropy             entropy based on fraction of different type segments
15      balanceRatio        balance ratio based on number of left- & right-clip segments
16      lowMapQualFrac      fraction of segments with low mapQual (<5)
17      dualClipFrac        fraction of "dual-clip" alignments
18      alnFrac1            fraction of segments with alnLocationType=1
19      alnFrac2            fraction of segments with alnLocationType=2
20      alnFrac4            fraction of segments with alnLocationType=4
21      alnFrac8            fraction of segments with alnLocationType=8
22      alnFrac16           fraction of segments with alnLocationType=16
23      meanMapQual         mean mapQual of cluster
24      meanAlnScore        mean per-base alignment score (based on teAlignments)
25      meanQueryMapFrac    mean query mapped fraction (based on teAlignments)
26      meanDivergence      mean per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
27      bgDiv               background divergence (for normalization)
28      bgDepth             background depth (for normalization)
29      bgReadLen           background read length
30      teAlignedFrac       fraction of TE-aligned segments
31      teTid               majority TE-tid of cluster
32      isInBlacklist       whether cluster intersects with blacklist
33      probability         the probability of the cluster to be a positive insertion
```


## 配置文件说明

### 完整参数列表

```yaml
# === 基本信息 ===
species: "dm6"                # 物种标识
mode: 1                       # 1=训练, 2=测试
share_pgdf: true              # 是否跨测序协议共享 PGD

# === 路径配置 ===
output_dir: "/path/to/output"
reference:
  template: "/path/to/template.fa"
  transposon: "/path/to/transposon.fa"
  repeat_bed: "/path/to/repeat.bed"
  gap_bed: "/path/to/gap.bed"
  blacklist: "/path/to/blacklist.bed"
  germ_model: "/path/to/germline_model"
  soma_model: "/path/to/somatic_model"

# === Simulation 参数 ===
simulation:
  population_size: 100           # 群体大小
  sub_population_size: 5         # 子群体大小
  germline_count: 1000           # Germline insertion 数量
  avg_somatic_count: 50          # 平均 somatic insertion 数量
  min_distance: 9000             # TE insertion 最小间距 (bp)
  divergence_rate: 0.001         # 序列分歧率
  trunc_prob: 0.1                # Truncation 概率
  nest_prob: 0.1                 # Nested insertion 概率

# === 测序参数 ===
sequencing:
  depths: [10, 20, 30, 40, 50]   # 测序深度列表
  protocols: [ccs, clr, ont]     # 测序协议
  ngs:
    read_length: 150
    insert_size: 350
    insert_std: 50

# === 数据集配置 ===
datasets:
  10: 20       # 深度: 数据集数量
  20: 20
  30: 20

# === SLURM 配置 ===
slurm:
  partition: "compute"
  threads:
    simulation: 4
    merge: 1
    minimap2: 16
    label: 16
  mem:
    simulation: "16G"
    merge: "4G"
    minimap2: "32G"
    label: "64G"
  time: "24:00:00"
  array_limit: 100

# === Conda 环境 ===
conda_envs:
  simulation: "simulation"
  merge: "simulation"
  minimap2: "simulation"
  label: "TEMP3"

# === 分类参数 ===
classification:
  somatic_freq: 0.1               # somatic 频率阈值

# === 程序路径 ===
programs:
  simulation_dir: "/path/to/simulation"
  label_dir: "/path/to/simulation/label"
  temp3_dir: "/path/to/TEMP3"
```

## 物种支持

### 已有物种特定逻辑

| 物种 | TE 选择权重 |
|------|------------|
| human | Alu 83%, LINE1 12%, SVA 4%, HERVK 1% |
| fly（dm6） | P_element 58% 等 |

### 添加新物种

新物种使用通用逻辑：
- **TE 选择**：从 `transposon.fa` 中均匀随机选择
- **Truncation**：使用 `trunc_prob` 参数控制
- **Frequency**：使用均匀分布 U(0.1, 1.0)

在 `species_configs/` 下创建新的 YAML 配置文件即可。

## 训练数据 vs 测试数据

`simulate_training_data/` 和 `simulate_test_data/` 目录分别包含 GRCh38 和 dm6 的训练/测试数据生成脚本：

- **训练数据**：用于训练 AutoGluon 分类模型，包含完整的 labeling 和分类流程
- **测试数据**：用于评估模型性能，分为高频和低频 insertion 两类

模型训练脚本位于 `model_training/training.py`，使用 AutoGluon Tabular 进行二分类训练。

## 常见问题

### Q: 如何为新物种准备参考文件？

1. **template.fa**：从参考基因组中提取主要染色体，移除 alternate contigs
2. **transposon.fa**：从 RepeatMasker 结果中提取 TE consensus 序列
3. **repeat.bed**：RepeatMasker 输出的 BED 格式
4. **gap.bed**：从参考基因组中识别 gap 区域

### Q: SLURM 任务失败怎么办？

检查日志文件 `<output_dir>/logs/`，常见问题：
- 内存不足：调整 `slurm.mem_*` 参数
- 时间不足：调整 `slurm.time` 参数
- Conda 环境问题：检查 `conda_envs` 配置是否正确

### Q: 如何只生成特定深度的数据集？

修改配置文件中的 `datasets` 部分，例如只生成 30X：

```yaml
datasets:
  30: 20
```

### Q: 如何跳过某些 step？

使用 `--start` 参数从指定步骤开始，或单独使用 `--step` 执行特定步骤。

### Q: share_pgdf 是什么意思？

启用后，不同测序协议（ccs/clr/ont）共享同一套 population genome definition，避免重复生成 PGD 文件。
