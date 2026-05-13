#!/usr/bin/env python3
"""
生成 simulation pipeline 所需的 filelists

根据配置文件生成:
- simulation_filelist: 用于 simulation_with_slurm.sh
- merge_filelist: 用于 merge_with_slurm.sh
- minimap2_filelist: 用于 minimap2_with_slurm.sh
- label_filelist: 用于 label_with_slurm.sh

用法:
    python generate_filelists.py --config species_configs/dm6.yaml --output-dir ./filelists
"""

import argparse
import os
import yaml
from pathlib import Path


def load_config(config_path):
    """加载 YAML 配置文件"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def get_contigs(template_path):
    """从 template.fa.fai 获取 contigs 信息"""
    fai_path = template_path + '.fai'
    contigs = []
    genome_size = 0

    with open(fai_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                contig_name = parts[0]
                contig_size = int(parts[1])
                contigs.append((contig_name, contig_size))
                genome_size += contig_size

    return contigs, genome_size


def generate_simulation_filelist(config, output_dir):
    """
    生成 simulation_filelist
    格式: work_dir te_fa genome_size contig_size pop_size sub_pop_size depth protocol
    """
    filelist_path = output_dir / 'simulation_filelist'

    output_base = Path(config['output_dir'])
    template_path = config['reference']['template']
    te_fa = config['reference']['transposon']

    contigs, genome_size = get_contigs(template_path)

    pop_size = config['simulation']['population_size']
    sub_pop_size = config['simulation']['sub_population_size']

    with open(filelist_path, 'w') as f:
        for protocol in config['sequencing']['protocols']:
            for depth, count in config['datasets'].items():
                # 计算 dataset ID 范围
                # 例如 depth=50 的 20 个数据集: ID 1-20
                start_id = 1
                for d, c in config['datasets'].items():
                    if d == depth:
                        break
                    start_id += c

                for dataset_id in range(start_id, start_id + count):
                    dataset_name = f"{dataset_id}_{protocol}"

                    for contig_name, contig_size in contigs:
                        work_dir = output_base / dataset_name / contig_name
                        f.write(f"{work_dir} {te_fa} {genome_size} {contig_size} {pop_size} {sub_pop_size} {depth} {protocol}\n")

    return filelist_path


def generate_merge_filelist(config, output_dir):
    """
    生成 merge_filelist
    格式: dataset_dir
    """
    filelist_path = output_dir / 'merge_filelist'

    output_base = Path(config['output_dir'])

    total_datasets = sum(config['datasets'].values())

    with open(filelist_path, 'w') as f:
        for protocol in config['sequencing']['protocols']:
            for dataset_id in range(1, total_datasets + 1):
                dataset_name = f"{dataset_id}_{protocol}"
                dataset_dir = output_base / dataset_name
                f.write(f"{dataset_dir}\n")

    return filelist_path


def generate_minimap2_filelist(config, output_dir):
    """
    生成 minimap2_filelist
    格式: ref_genome reads_fasta preset

    每个 dataset 根据 protocol 使用对应的 preset:
    - ccs → map-hifi
    - clr → map-pb
    - ont → map-ont
    """
    filelist_path = output_dir / 'minimap2_filelist'

    output_base = Path(config['output_dir'])
    ref_genome = config['reference']['template']

    # minimap2 preset 映射: 每个 protocol 使用对应 preset
    preset_map = {
        'ccs': 'map-hifi',
        'clr': 'map-pb',
        'ont': 'map-ont'
    }

    total_datasets = sum(config['datasets'].values())

    with open(filelist_path, 'w') as f:
        for protocol in config['sequencing']['protocols']:
            for dataset_id in range(1, total_datasets + 1):
                dataset_name = f"{dataset_id}_{protocol}"
                reads_fasta = output_base / dataset_name / 'TGS.fasta'
                preset = preset_map[protocol]
                f.write(f"{ref_genome} {reads_fasta} {preset}\n")

    return filelist_path


def generate_label_filelist(config, output_dir):
    """
    生成 label_filelist
    格式: bam_path repeat_bed gap_bed ref_te

    每个 dataset 根据 protocol 使用对应的 preset:
    - ccs → map-hifi
    - clr → map-pb
    - ont → map-ont
    """
    filelist_path = output_dir / 'label_filelist'

    output_base = Path(config['output_dir'])
    repeat_bed = config['reference']['repeat_bed']
    gap_bed = config['reference']['gap_bed']
    ref_te = config['reference']['transposon']

    # preset 映射: 每个 protocol 使用对应 preset
    preset_map = {
        'ccs': 'map-hifi',
        'clr': 'map-pb',
        'ont': 'map-ont'
    }

    total_datasets = sum(config['datasets'].values())

    with open(filelist_path, 'w') as f:
        for protocol in config['sequencing']['protocols']:
            for dataset_id in range(1, total_datasets + 1):
                dataset_name = f"{dataset_id}_{protocol}"
                preset = preset_map[protocol]
                bam_path = output_base / dataset_name / preset / 'TGS.bam'
                f.write(f"{bam_path} {repeat_bed} {gap_bed} {ref_te}\n")

    return filelist_path


def count_tasks(filelist_path):
    """计算 filelist 中的任务数"""
    with open(filelist_path, 'r') as f:
        return sum(1 for line in f if line.strip())


def main():
    parser = argparse.ArgumentParser(description='生成 simulation pipeline filelists')
    parser.add_argument('--config', required=True, help='配置文件路径 (YAML)')
    parser.add_argument('--output-dir', default='./filelists', help='filelists 输出目录')
    parser.add_argument('--dry-run', action='store_true', help='只显示统计信息，不生成文件')

    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(args.output_dir)

    print(f"物种: {config['species']}")
    print(f"输出目录: {config['output_dir']}")
    print(f"数据集总数: {sum(config['datasets'].values())} x {len(config['sequencing']['protocols'])} protocols")
    print()

    # 获取 contigs 信息
    template_path = config['reference']['template']
    contigs, genome_size = get_contigs(template_path)
    print(f"基因组大小: {genome_size}")
    print(f"Contigs 数量: {len(contigs)}")
    print()

    # 计算各 filelist 的任务数
    print("任务统计:")

    # simulation_filelist
    simulation_tasks = 0
    for protocol in config['sequencing']['protocols']:
        for depth, count in config['datasets'].items():
            simulation_tasks += count * len(contigs)

    print(f"  simulation_filelist: {simulation_tasks} tasks")

    # merge_filelist
    merge_tasks = sum(config['datasets'].values()) * len(config['sequencing']['protocols'])
    print(f"  merge_filelist: {merge_tasks} tasks")

    # minimap2_filelist: 每个 dataset 1 个 preset
    minimap2_tasks = merge_tasks
    print(f"  minimap2_filelist: {minimap2_tasks} tasks")

    # label_filelist: 每个 dataset 1 个 preset
    label_tasks = merge_tasks
    print(f"  label_filelist: {label_tasks} tasks")

    if args.dry_run:
        print("\n[dry-run] 未生成文件")
        return

    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)

    # 生成 filelists
    print("\n生成 filelists:")

    sim_path = generate_simulation_filelist(config, output_dir)
    print(f"  {sim_path}: {count_tasks(sim_path)} lines")

    merge_path = generate_merge_filelist(config, output_dir)
    print(f"  {merge_path}: {count_tasks(merge_path)} lines")

    minimap2_path = generate_minimap2_filelist(config, output_dir)
    print(f"  {minimap2_path}: {count_tasks(minimap2_path)} lines")

    label_path = generate_label_filelist(config, output_dir)
    print(f"  {label_path}: {count_tasks(label_path)} lines")

    print("\n完成!")


if __name__ == '__main__':
    main()