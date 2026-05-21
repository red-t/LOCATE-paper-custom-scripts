## 1. Generate raw reads

### 1.1 Generate population genome definition files on local machine

```shell
# 1. For High-frequency Insertion Testing Datasets
mkdir testing_dataset_0 && cd testing_dataset_0
bash your_path_to/simulation/simulation_protocol.sh -d ./ -r your_path_to/Dm6_no_alt.fa -t your_path_to/Dm6.transposon.fa -N 100 -R 0 --sub-N 5 --germline-count 500 --avg-somatic-count 0 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2

# 2. For Low-frequency Insertion Testing Datasets
mkdir testing_dataset_1 && cd testing_dataset_1
bash your_path_to/simulation/simulation_protocol.sh -d ./ -r your_path_to/Dm6_no_alt.fa -t your_path_to/Dm6.transposon.fa -N 1000 -R 0 --sub-N 50 --germline-count 0 --avg-somatic-count 10 --min-distance 9000 --depth 50 --species human --protocol ccs --truncProb 0.7 --nestProb 0 --mode 2
```

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
path_to_testing_dataset_x
```

Then, merge results:

```shell
sbatch merge_with_blades_dm6.sh

# Main results under testing_dataset_x of each datasets:
# 1. TGS.fasta
# 2. NGS_1.fasta  NGS_2.fasta
# 3. chr2L/chr2L.ins.summary
```

3 datasets was generated for testing datasets using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 100                    | 500                      | 0                       | 0                         | 0.1~1                       | 0                           | 0.1                            | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

3 datasets was generated for testing somatic insertion using different `protocol(ccs, clr, ont)`:

| Population genome size | Germline insertion count | Somatic insertion count | Average somatic insertion | Gemline insertion frequency | Somatic insertion frequency | Insertion truncate probability | Nested insertion probability | Insertion divergence rate | TSD length | Depth | Min distance | TGS read length                  | TGS error rate              |
| ---------------------- | ------------------------ | ----------------------- | ------------------------- | --------------------------- | --------------------------- | ------------------------------ | ---------------------------- | ------------------------- | ---------- | ----- | ------------ | -------------------------------- | --------------------------- |
| 1000                   | 0                        | 10000                   | 10                        | 0                           | 0.001                       | 0.1                            | 0                            | 0                         | 6~20       | 50X   | 9000         | CCS~=13490; CLR~=7896; ONT~=7170 | CCS=0.01; CLR=0.1; ONT=0.07 |

