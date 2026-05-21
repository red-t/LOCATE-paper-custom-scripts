## 1. Generate raw reads

### 1.1 Generate population genome definition files on local machine

```shell
mkdir training_dataset_0 && cd training_dataset_0

bash your_path_to/simulation/simulation_protocol.sh -d ./ -r your_path_to/Dm6_no_alt.fa -t your_path_to/Dm6.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 1000 --avg-somatic-count 50 --min-distance 9000 --depth 50 --species fly --protocol ccs --truncProb 0.1 --nestProb 0.1 --mode 1
```

Repeat this step `N` times with different parameters `--protocol ccs/clr/ont`, to genereate `N` different datasets.

### 1.2 Build population genome and generate raw reads

Prepare the filelist based on the results of `Step 1.1`, for example, `filelist_example/simulation_filelist_dm6`. The required information include:

```shell
path_to_chrom
path_to_transposon_consensus
total_genome_size
chrom_size
total_haplotypes_to_be_simulated # population genome size
haplotypes_to_be_simulated_per_batch
average_coverage
TGS_platform
```

Then, build population genome(s) and generate reads:

```shell
sbatch simulation_with_blades_dm6.sh
```

### 1.3 Merge raw reads from each chromosome

Prepare the required filelist for merging, for example, `filelist_example/merge_filelist_dm6`. The required information include:

```shell
path_to_training_dataset_x
```

Then, merge results:

```shell
sbatch merge_with_blades_dm6.sh

# Main results under training_dataset_x of each datasets:
# 1. TGS.fasta
# 2. chr2L/chr2L.ins.summary
# 3. chr2L/chr2L.ins.sequence
```

In our research, total 300 datasets were generated for Drosophila using different `depth(10-50X)` and `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth  | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ------ | ------------ | -------------------------------- | --------------------------- |
| 100                    | 1000                     | 5000                    | 50                        | 0.1~1                       | 0.01                        | 0.1                            | 0.1                          | 0                         | 6~8        | 10~50X | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

## 2. Labeling

This step aims to find all insertion candidates and label them based on ground truth. In other word, define TP/FP.

### 2.1 Map raw reads to reference genome with minimap2

Prepare the required filelist for mapping, for example, `filelist_example/minimap2_filelist_dm6`. The required information include:

```shell
path_to_reference_genome
path_to_simulated_reads
minimap2_preset # map-hifi for ccs, map-pb for clr, map-ont for ont
```

Then, map reads for each datasets:

```shell
sbtach minimap2_with_blades_dm6.sh
```

### 2.2 Dtecet all candidates on local machine

Prepare the required filelist for detection, for example, `filelist_example/label_filelist_dm6`. The required information include:

```shell
path_to_alignment
path_to_repeat_annotation
path_to_gap_annotation
path_to_transposon_consensus
```

Then, detect candidates for each datasets:

```shell
bash batch_labeling_dm6.sh
```

### 2.3 Label all candidates on blades

Prepare the required filelist for labeling, for example, `filelist_example/label_filelist_dm6`.

Then, label candidates for each datasets:

```shell
sbatch label_with_blades_dm6.sh
```

### 2.4 Classify high-/low- frequency insertion

This step is only required when building the training dataset, which aims to classify insertions with high-/low- frequency

```shell
bash classify_germline_and_smoatic.sh

# Main results under training_dataset_x of each datasets:
# 1. map-hifi/result_no_secondary/TP_clt_G.txt
# 2. map-hifi/result_no_secondary/TP_clt_S.txt
# 3. map-hifi/result_no_secondary/TP_clt.txt
# 4. map-hifi/result_no_secondary/FP_clt.txt
```
