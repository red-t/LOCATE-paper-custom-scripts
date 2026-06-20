# LOCATE — TE Insertion Simulation Pipeline

A general-purpose pipeline for generating simulated TE (Transposable Element) insertion data. Supports multiple species and sequencing protocols, capable of generating both training and testing data.

## Repository Structure

```
LOCATE-paper-custom-scripts/
├── simulation/                    # Simulation data generation pipeline
│   ├── run_simulation.py          # Main orchestration script (7-step pipeline)
│   ├── generate_filelists.py      # SLURM task filelist generation
│   ├── define_population_genome.py
│   ├── build-population-genome.py
│   ├── generate_TGS.py            # TGS reads generation
│   ├── generate_NGS.py            # NGS reads generation
│   ├── config_template.yaml       # Configuration file template
│   ├── species_configs/           # Species-specific configurations
│   │   ├── dm6.yaml               # D. melanogaster
│   │   ├── rice.yaml              # Rice (template)
│   │   └── test.yaml              # Minimal test configuration
│   ├── scripts/                   # SLURM job scripts
│   │   ├── simulation_with_slurm.sh
│   │   ├── merge_with_slurm.sh
│   │   ├── minimap2_with_slurm.sh
│   │   └── label_with_slurm.sh
│   ├── label/                     # Labeling and classification scripts
│   │   ├── label_protocol.sh
│   │   ├── label_protocol_local.sh
│   │   ├── Filter_TP_and_FP.py
│   │   ├── rename.py
│   │   ├── merge_TP_and_FP.sh
│   │   ├── classify_germline_somatic.py
│   │   └── classify_germline_somatic.sh
│   ├── simulate_training_data/    # Training data generation scripts (GRCh38 / dm6)
│   └── simulate_test_data/        # Test data generation scripts (GRCh38 / dm6)
├── TEMP3/                         # TE insertion detection tool
│   ├── TEMP3.py                   # Main entry point
│   ├── TEMP3/                     # Core modules (Cython + C)
│   └── setup.py
└── model_training/                # AutoGluon classification model training
    └── pipeline/                  # ML training pipeline (Step 7)
        ├── run_pipeline.py        # Standalone orchestration script
        ├── config.py              # YAML configuration parser
        ├── defaults.py            # Column name constants and mode defaults
        ├── utils.py               # Utility functions
        ├── evaluate.py            # Shared evaluation functions
        ├── stage1_preprocess.py   # Data preprocessing (merge / split / blacklist / dedup)
        ├── stage2_shuffle.py      # Proportional P/N data sampling
        ├── stage3_train.py        # Train AutoGluon models across ratios
        └── stage4_filter_eval.py  # Blacklist filtering + model evaluation
```

## Environment Setup

### 1. Create simulation environment

```shell
cd simulation/
mamba env create -f simulation.env.yaml
python setup.py build_ext -i && rm -r build && rm -f *.c
```

### 2. Create TEMP3 environment (for labeling)

```shell
cd TEMP3/
mamba env create -f TEMP3.env.yaml
python setup.py build_ext -i && rm -r build && rm -f *.c
```

## Quick Start

### 1. Prepare reference files

Each species requires the following reference files:

| File | Description | Requirements |
|------|------------|--------------|
| `template.fa` | Reference genome template | Requires samtools faidx index |
| `transposon.fa` | TE consensus sequences | Requires samtools faidx index |
| `repeat.bed` | RepeatMasker annotations | BED format |
| `gap.bed` | Gap region annotations | BED format |

### 2. Create configuration file

```shell
cp config_template.yaml species_configs/<species>.yaml
# Edit the configuration file with reference file paths and parameters
```

### 3. Run Pipeline

```shell
# View all steps
python run_simulation.py --config species_configs/dm6.yaml

# Step 0: Define population genome (local execution)
python run_simulation.py --config species_configs/dm6.yaml --step 0

# Step 1: Generate filelists
python run_simulation.py --config species_configs/dm6.yaml --step 1

# Step 2-5: Submit SLURM jobs sequentially
python run_simulation.py --config species_configs/dm6.yaml --step 2
python run_simulation.py --config species_configs/dm6.yaml --step 3
python run_simulation.py --config species_configs/dm6.yaml --step 4
python run_simulation.py --config species_configs/dm6.yaml --step 5

# Step 6: Classify germline/somatic (local execution)
python run_simulation.py --config species_configs/dm6.yaml --step 6

# Step 7: ML model training (local execution)
python run_simulation.py --config species_configs/dm6.yaml --step 7

# Or skip preceding steps and only execute Step 7:
python run_simulation.py --config species_configs/dm6.yaml --start 7 --step 7

# Step 7 can also run independently (not requiring run_simulation.py):
python run_simulation.py --config species_configs/dm6.yaml --step 7 --dry-run
```

Dry-run mode (show commands without executing):

```shell
python run_simulation.py --config species_configs/dm6.yaml --step 0 --dry-run
```

## Pipeline Step Details

| Step | Name | Description | Execution | Core Script |
|------|------|-------------|-----------|-------------|
| 0 | define_pgdf | Define population genome, generate PGD files | Local | `define_population_genome.py` |
| 1 | generate_filelists | Generate filelists required by SLURM tasks | Local | `generate_filelists.py` |
| 2 | simulation | Build population genome, generate TGS reads | SLURM | `scripts/simulation_with_slurm.sh` |
| 3 | merge | Merge per-chromosome reads | SLURM | `scripts/merge_with_slurm.sh` |
| 4 | minimap2 | Align reads to reference genome | SLURM | `scripts/minimap2_with_slurm.sh` |
| 5 | label | Detect TE insertions and classify TP/FP | SLURM | `scripts/label_with_slurm.sh` → TEMP3 + `label/Filter_TP_and_FP.py` |
| 6 | classify | Classify germline/somatic insertions | Local | `label/classify_germline_somatic.py` |
| 7 | ml_training | ML model training (AutoGluon) | Local | `model_training/pipeline/run_pipeline.py` |

### Step 0 — define_pgdf

Defines the Population Genome Definition. Generates `.pgd` files for each contig, recording all TE insertion positions, types, and frequency information.

### Step 1 — generate_filelists

Generates four filelist files for SLURM tasks based on configuration:
- `simulation_filelist` — input for Step 2
- `merge_filelist` — input for Step 3
- `minimap2_filelist` — input for Step 4
- `label_filelist` — input for Step 5

### Step 2 — simulation

Builds the population genome and generates TGS reads. Each SLURM task processes one contig/sub-population combination and generates reads for the corresponding chromosome.

### Step 3 — merge

Merges per-contig TGS reads into a single `TGS.fasta` file.

### Step 4 — minimap2

Aligns reads to the reference genome using minimap2, with protocol-specific presets (map-hifi / map-pb / map-ont). Outputs sorted BAM files.

### Step 5 — label

Runs TEMP3 to detect TE insertion candidates, then compares against ground truth using `Filter_TP_and_FP.py` to classify TP (true positive) and FP (false positive).

### Step 6 — classify

Classifies TP insertions into germline (high-frequency) and somatic (low-frequency) based on insertion frequency. Outputs `TP_clt_G.txt` and `TP_clt_S.txt`.

### Step 7 — ml_training

Trains AutoGluon classification models using labeled data from Steps 5-6. Controlled by the `ml_training` section in the YAML configuration.

Workflow:
1. **preprocess** — iterate through filelists, read per-sample `TP_clt_G.txt`/`TP_clt_S.txt`/`FP_clt.txt`, filter by `cltType` and `teAlignedFrac`, insert sample_id, merge
2. **shuffle** — sample positive and negative data according to configured P/N ratios (e.g., ORG, 1V1, 1V30)
3. **train** — train AutoGluon TabularPredictor models for each ratio
4. **evaluate** — filter test data with BlackList, evaluate models, output summary

Supports two modes:
- **germline** (high-frequency): `cltType==0`, `teAlignedFrac>=0.8`, positive sources from `TP_clt_G.txt` only
- **somatic** (low-frequency): `cltType>0`, `teAlignedFrac>=1`, positive sources from `TP_clt_G.txt` + `TP_clt_S.txt`

Output directory structure:
```
{output_dir}/For_ML/{species}/
  germline/
    Train/         train_P/N_{ratio}.txt
    Test/          test_P/N_{protocol}_{ratio}.txt
    BlackList*.bed
    {species}_G_{ratio}/        # AutoGluon model directory
    {species}_G_summary_Dedup.txt
  somatic/
    Train/
    Test/
    BlackList*.bed
    {species}_S_1V30/
    {species}_S_summary_Dedup.txt
```

## Output Files

Main output structure per dataset:

```
<output_dir>/<dataset_id>_<protocol>/
├── TGS.fasta                         # Simulated TGS reads
├── <chrom>/
│   ├── <chrom>.ins.summary           # Insertion summary
│   └── <chrom>.ins.sequence          # Insertion sequences
├── map-hifi/TGS.bam                  # HiFi alignment results
├── map-pb/TGS.bam                    # PacBio CLR alignment results
├── map-ont/TGS.bam                   # ONT alignment results
└── result_no_secondary/
    ├── TP_clt.txt                    # All true positives
    ├── TP_clt_G.txt                  # Germline TP
    ├── TP_clt_S.txt                  # Somatic TP
    └── FP_clt.txt                    # False positives
```

### Output File Fields

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
                                other: unknown
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

## Configuration File Reference

### Complete Parameter List

```yaml
# === Basic information ===
species: "dm6"                # Species identifier
mode: 1                       # 1=training, 2=testing
share_pgdf: true              # Share PGD across sequencing protocols

# === Path configuration ===
output_dir: "/path/to/output"
reference:
  template: "/path/to/template.fa"
  transposon: "/path/to/transposon.fa"
  repeat_bed: "/path/to/repeat.bed"
  gap_bed: "/path/to/gap.bed"
  blacklist: "/path/to/blacklist.bed"
  germ_model: "/path/to/germline_model"
  soma_model: "/path/to/somatic_model"

# === Simulation parameters ===
simulation:
  population_size: 100           # Population size
  sub_population_size: 5         # Sub-population size
  germline_count: 1000           # Number of germline insertions
  avg_somatic_count: 50          # Average number of somatic insertions
  min_distance: 9000             # Minimum TE insertion distance (bp)
  divergence_rate: 0.001         # Sequence divergence rate
  trunc_prob: 0.1                # Truncation probability
  nest_prob: 0.1                 # Nested insertion probability

# === Sequencing parameters ===
sequencing:
  depths: [10, 20, 30, 40, 50]   # Sequencing depths
  protocols: [ccs, clr, ont]     # Sequencing protocols
  ngs:
    read_length: 150
    insert_size: 350
    insert_std: 50

# === Dataset configuration ===
datasets:
  10: 20       # Depth: number of datasets
  20: 20
  30: 20

# === SLURM configuration ===
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

# === Conda environment ===
conda_envs:
  simulation: "simulation"
  merge: "simulation"
  minimap2: "simulation"
  label: "TEMP3"

# === Classification parameters ===
classification:
  somatic_freq: 0.1               # Somatic frequency threshold

# === Program paths ===
programs:
  simulation_dir: "/path/to/simulation"
  label_dir: "/path/to/simulation/label"
  temp3_dir: "/path/to/TEMP3"
```

## Species Support

### Species-Specific Logic

| Species | TE Selection Weights |
|---------|---------------------|
| human | Alu 83%, LINE1 12%, SVA 4%, HERVK 1% |
| fly (dm6) | P_element 58%, etc. |

### Adding a New Species

New species use general logic:
- **TE selection**: uniform random selection from `transposon.fa`
- **Truncation**: controlled by `trunc_prob` parameter
- **Frequency**: uniform distribution U(0.1, 1.0)

Create a new YAML configuration file in `species_configs/` to add a new species.

## Training Data vs Testing Data

The `simulate_training_data/` and `simulate_test_data/` directories contain data generation scripts for GRCh38 and dm6:

- **Training data**: used to train AutoGluon classification models, includes complete labeling and classification workflow
- **Testing data**: used to evaluate model performance, categorized into high-frequency and low-frequency insertions

Model training scripts are located in `model_training/training.py`, using AutoGluon Tabular for binary classification.

## FAQ

### Q: How do I prepare reference files for a new species?

1. **template.fa**: extract main chromosomes from the reference genome, remove alternate contigs
2. **transposon.fa**: extract TE consensus sequences from RepeatMasker output
3. **repeat.bed**: RepeatMasker output in BED format
4. **gap.bed**: identify gap regions from the reference genome

### Q: What should I do if a SLURM task fails?

Check log files in `<output_dir>/logs/`. Common issues:
- Insufficient memory: adjust `slurm.mem_*` parameters
- Insufficient time: adjust `slurm.time` parameters
- Conda environment issues: verify `conda_envs` configuration

### Q: How do I generate datasets for only specific depths?

Modify the `datasets` section in the configuration file, e.g., for 30X only:

```yaml
datasets:
  30: 20
```

### Q: How do I skip certain steps?

Use the `--start` parameter to begin from a specific step, or `--step` to execute a single step.

### Q: What does share_pgdf mean?

When enabled, different sequencing protocols (ccs/clr/ont) share the same population genome definition, avoiding redundant PGD file generation.
