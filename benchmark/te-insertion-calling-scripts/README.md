# LOCATE benchmark example scripts

This directory contains example scripts used to run LOCATE and comparison tools on one simulation data set. The scripts are intended as reproducible workflow examples for publication. They are not a packaged one-command pipeline, and users should adjust paths, sample names, depths, scheduler options, and tool environments for their own systems.

## Path placeholders

All shell scripts define two placeholders near the top:

```bash
path1="${path1:-path1}"
path2="${path2:-path2}"
```

- `path1`: reference and annotation root, including genome FASTA files, TE libraries, repeat annotation files, gap/blacklist files, and LOCATE model files.
- `path2`: project root, including raw simulation data, results, locally installed third-party tools, conda profile files, and job templates.

You can either edit these two lines in each script or pass them at runtime:

```bash
path1=/path/to/annotations path2=/path/to/project bash run_LOCATE_simulation_hg38.sh
```

The TrEMOLO YAML templates use `template_path1`; the TrEMOLO shell scripts replace it with `path1` after copying the template.

## Files

### LOCATE runs

- `run_LOCATE_simulation_hg38.sh`: LOCATE on the hg38 ONT simulation example.
- `run_LOCATE_simulation_dm6.sh`: LOCATE on the dm6 ONT simulation example.
- `run_LOCATE_simulation_rice.sh`: LOCATE on the rice `IRGSPv1` ONT simulation example.

### Read mapping and downsampling

- `run_map_and_downsmaple_new-ont-hg38.sh`: map and downsample hg38 ONT reads.
- `run_map_and_downsmaple_new-ont-dm6.sh`: map and downsample dm6 ONT reads.
- `run_map_and_downsmaple_new-ont-rice.sh`: map and downsample rice ONT reads.

### Comparison tools

- `run_simulation_hg38_tgs.tldr.sh`: TLDR on hg38 long-read simulations.
- `run_simulation_hg38_tgs.telr.sh`: TELR on hg38 long-read simulations.
- `run_simulation_hg38_tgs.TreMOLO.sh`: TrEMOLO on hg38 long-read simulations.
- `run_simulation_dm6_tgs.sh`: TrEMOLO-focused dm6 long-read simulation script.
- `run_simulation_rice_tgs.sh`: TLDR, TELR, and TrEMOLO on rice long-read simulations.
- `run_GraffiTE_simulation_hg38.sh`, `run_GraffiTE_simulation_dm6.sh`, `run_GraffiTE_simulation_rice.sh`: GraffiTE comparison runs.
- `run_MEHunter_hg38.sh`: MEHunter comparison run for hg38.
- `run_palmer.sh`: PALMER comparison run.
- `run_simulation_rice_NGS_TEMP2.sh`: TEMP2 comparison run for rice short-read simulations.
- `run_simulation_rice_NGS_MELT.sh`: MELT comparison run for rice short-read simulations.
- `run_xTea_giab.sh`: xTea command example for GIAB-style input lists.

### TrEMOLO templates

- `config_template_hg38.yaml`
- `config_template_dm6.yaml`
- `config_template_rice.yaml`
- `config_template_giab.yaml`

These templates contain `template_sample_fa`, `template_outpath`, `template_map_mode`, and `template_path1` placeholders. The corresponding TrEMOLO shell scripts replace them before running `snakemake`.

## Expected input layout

The scripts assume the example data are organized under `path2` with subdirectories such as:

```text
path2/
  rawdata/
    simulation/
      hg38/
      dm6/
      IRGSPv1/
  result/
    simulation/
  tools/
```

The annotation files are expected under `path1`, for example:

```text
path1/
  hg38/
  dm6/
  IRGSPv1/
  from_temp3/
```

These names mirror the sample data used for the benchmark. If your local annotation layout differs, update the variables near the top of the relevant script.

## Software requirements

Install the external tools before running the corresponding scripts:

- LOCATE
- minimap2, samtools, seqtk, bedtools, bwa
- TLDR
- TELR
- TrEMOLO and snakemake
- GraffiTE and nextflow
- MEHunter and cuteSV
- PALMER
- TEMP2
- MELT
- xTea

The scripts activate existing conda environments such as `locate_pub`, `locate_pub2`, `tldr`, `TELR`, `TrEMOLO`, `py2`, `cuteSVenv`, and `MEHunterEnv`. Rename these environment names if your installation uses different names.

## Running an example

Review and edit `path1`, `path2`, sample names, sequencing depths, CPU/thread settings, and tool-specific paths first. Then run the script for the desired genome or tool:

```bash
path1=/path/to/annotations path2=/path/to/project bash run_LOCATE_simulation_hg38.sh
```

Some scripts include cleanup commands for generated intermediate files. Review those commands before running the scripts on directories that already contain results you need to keep.

For publication use, run a static syntax check and a private-path scan before committing:

```bash
bash -n *.sh
rg -n -S "(/[d]ata/|/[z]ata/|/[h]ome/|/[U]sers/|/[V]olumes/|/[m]nt/|[b]oxu|[x]ubo|[z]hongrenhu|[t]users)" .
```
