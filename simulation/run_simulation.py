#!/usr/bin/env python3
"""
Simulation Pipeline 主控脚本

分步执行模式:
1. Step 0: 定义群体基因组 (define_pgdf)
2. Step 1: 生成 filelists
3. Step 2: 运行 simulation (SLURM)
4. Step 3: 合并 reads (SLURM)
5. Step 4: minimap2 比对 (SLURM)
6. Step 5: 标注 TP/FP (SLURM)
7. Step 6: 分类 germline/somatic

用法:
    python run_simulation.py --config species_configs/dm6.yaml --step 0-6
    python run_simulation.py --config species_configs/dm6.yaml --dry-run
"""

import argparse
import os
import re
import shutil
import subprocess
import time
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

    with open(fai_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                contig_name = parts[0]
                contigs.append(contig_name)

    return contigs


def get_script_dir():
    """获取脚本目录"""
    script_path = Path(__file__).resolve()
    return script_path.parent


def get_conda_env(config, step, default):
    """获取指定 step 的 conda 环境，支持按 step 覆盖"""
    conda_envs = config.get('conda_envs', {})
    if isinstance(conda_envs, dict):
        return conda_envs.get(step, conda_envs.get('default', default))
    return default


def parse_slurm_job_id(stdout):
    """从 sbatch 输出中提取 job ID"""
    match = re.search(r'Submitted batch job (\d+)', stdout)
    return match.group(1) if match else None


def wait_for_slurm_job(job_id, check_interval=30):
    """轮询等待 SLURM 任务完成"""
    print(f"等待 SLURM 任务 {job_id} 完成...")
    while True:
        result = subprocess.run(['squeue', '-j', str(job_id)],
                                capture_output=True, text=True)
        if result.returncode != 0 or len(result.stdout.strip().split('\n')) <= 1:
            break
        time.sleep(check_interval)
    print(f"SLURM 任务 {job_id} 已完成")


def step0_define_pgdf(config, dry_run=False, **kwargs):
    """
    Step 0: 定义群体基因组

    支持两种模式：
    - share_pgdf=True: 所有 protocol 共享同一套 pgdf 文件（只生成一次）
    - share_pgdf=False: 每个 protocol 独立生成 pgdf
    """
    print("=" * 60)
    print("Step 0: 定义群体基因组 (define_pgdf)")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    template_path = config['reference']['template']
    te_fa = config['reference']['transposon']

    script_dir = get_script_dir()
    simulation_protocol = script_dir / 'simulation_protocol.sh'

    # simulation parameters
    pop_size = config['simulation']['population_size']
    sub_pop_size = config['simulation']['sub_population_size']
    germline_count = config['simulation']['germline_count']
    avg_somatic_count = config['simulation']['avg_somatic_count']
    min_distance = config['simulation']['min_distance']
    divergence_rate = config['simulation']['divergence_rate']
    trunc_prob = config['simulation']['trunc_prob']
    nest_prob = config['simulation']['nest_prob']
    mode = config['mode']

    total_datasets = sum(config['datasets'].values())
    protocols = config['sequencing']['protocols']

    # 是否共享 pgdf
    share_pgdf = config.get('share_pgdf', True)  # 默认共享
    print(f"share_pgdf 模式: {share_pgdf}")

    # 计算每个 depth 的 dataset ID 范围
    depth_ranges = {}
    start_id = 1
    for depth, count in config['datasets'].items():
        depth_ranges[depth] = (start_id, start_id + count - 1)
        start_id += count

    # 构建任务列表
    # 每个任务: (dataset_name, dataset_dir, protocol, cmd_or_None)
    # share_pgdf=True 时，只有第一个 protocol 需要 cmd，其他为 None（复制）
    tasks = []
    first_protocol = protocols[0]

    for depth, (start, end) in depth_ranges.items():
        for dataset_id in range(start, end + 1):
            for protocol in protocols:
                dataset_name = f"{dataset_id}_{protocol}"
                dataset_dir = output_base / dataset_name

                if share_pgdf and protocol != first_protocol:
                    # 共享模式：非首个 protocol 只需复制，不需要生成命令
                    tasks.append((dataset_name, dataset_dir, protocol, None))
                else:
                    # 需要生成 pgdf
                    cmd = [
                        str(simulation_protocol),
                        '-d', str(dataset_dir),
                        '-r', template_path,
                        '-t', te_fa,
                        '-N', str(pop_size),
                        '-R', str(divergence_rate),
                        '--sub-N', str(sub_pop_size),
                        '--germline-count', str(germline_count),
                        '--avg-somatic-count', str(avg_somatic_count),
                        '--min-distance', str(min_distance),
                        '--depth', str(depth),
                        '--species', config['species'],
                        '--protocol', protocol,
                        '--truncProb', str(trunc_prob),
                        '--nestProb', str(nest_prob),
                        '--mode', str(mode)
                    ]
                    tasks.append((dataset_name, dataset_dir, protocol, cmd))

    # 统计
    generate_count = sum(1 for t in tasks if t[3] is not None)
    copy_count = sum(1 for t in tasks if t[3] is None)
    print(f"将生成 {generate_count} 个 pgdf，复制 {copy_count} 个 pgdf")
    print()

    if dry_run:
        print("[dry-run] 命令示例:")
        for i, (name, ddir, protocol, cmd) in enumerate(tasks[:3]):
            if cmd:
                print(f"  {name}: {' '.join(cmd)}")
            else:
                print(f"  {name}: 复制 pgdf (share_pgdf 模式)")
        if len(tasks) > 3:
            print(f"  ... 还有 {len(tasks) - 3} 个数据集")
        return

    # Conda 环境
    conda_path = config.get('conda_path', '/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh')
    conda_env = get_conda_env(config, 'step0', 'simulation')

    # 设置环境变量
    env = os.environ.copy()
    env['CONDA_PATH'] = conda_path
    env['CONDA_ENV'] = conda_env

    # 执行任务
    for i, (name, dataset_dir, protocol, cmd) in enumerate(tasks):
        print(f"[{i+1}/{len(tasks)}] 处理 {name}...")
        dataset_dir.mkdir(parents=True, exist_ok=True)

        log_file = dataset_dir / 'define_pgdf.log'

        if cmd:
            # 生成 pgdf
            with open(log_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
                if result.returncode != 0:
                    print(f"  ERROR: {name} 生成失败 (return code {result.returncode})")
                    continue
        else:
            # 复制 pgdf (share_pgdf 模式)
            # 从第一个 protocol 的数据集复制所有 chr 目录
            source_name = name.replace(f"_{protocol}", f"_{first_protocol}")
            source_dir = output_base / source_name

            if not source_dir.exists():
                print(f"  ERROR: 源目录不存在: {source_dir}")
                continue

            with open(log_file, 'w') as f:
                f.write(f"复制 pgdf from {source_dir}\n")

                # 复制所有染色体目录
                for chr_dir in source_dir.iterdir():
                    if chr_dir.is_dir() and chr_dir.name.startswith('chr'):
                        dest_chr_dir = dataset_dir / chr_dir.name
                        if dest_chr_dir.exists():
                            f.write(f"  {chr_dir.name}: 已存在，跳过\n")
                        else:
                            shutil.copytree(chr_dir, dest_chr_dir)
                            f.write(f"  {chr_dir.name}: 已复制\n")

                # 复制 summary 文件（如果存在）
                for summary_file in source_dir.glob('*.ins.summary'):
                    shutil.copy2(summary_file, dataset_dir / summary_file.name)
                    f.write(f"  {summary_file.name}: 已复制\n")
                    
    print("\nStep 0 完成!")


def step1_generate_filelists(config, dry_run=False, **kwargs):
    """
    Step 1: 生成 filelists

    调用 generate_filelists.py 生成所有 filelists
    """
    print("=" * 60)
    print("Step 1: 生成 filelists")
    print("=" * 60)

    script_dir = get_script_dir()
    generate_script = script_dir / 'generate_filelists.py'
    config_path = Path(config['_config_path'])  # 原始配置文件路径

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'

    cmd = [
        'python', str(generate_script),
        '--config', str(config_path),
        '--output-dir', str(filelists_dir)
    ]

    print(f"输出目录: {filelists_dir}")

    if dry_run:
        print("[dry-run] 命令:")
        print(f"  {' '.join(cmd)}")
        return

    # 检查 filelists 是否已存在
    if filelists_dir.exists():
        existing_files = list(filelists_dir.glob('*'))
        if existing_files:
            print(f"  WARNING: filelists 目录已存在，包含 {len(existing_files)} 个文件")
            print(f"  将覆盖现有文件")

    subprocess.run(cmd)

    print("\nStep 1 完成!")


def step2_simulation(config, dry_run=False, **kwargs):
    """
    Step 2: 运行 simulation (SLURM)

    构建 population genome 并生成 reads
    """
    print("=" * 60)
    print("Step 2: 运行 simulation (SLURM)")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'
    filelist = filelists_dir / 'simulation_filelist'

    if not filelist.exists():
        print(f"ERROR: filelist 不存在: {filelist}")
        print("请先运行 Step 1")
        return None

    # 计算任务数
    with open(filelist, 'r') as f:
        task_count = sum(1 for line in f if line.strip())

    script_dir = get_script_dir()
    slurm_script = script_dir / 'scripts' / 'simulation_with_slurm.sh'

    # SLURM 参数
    slurm_config = config.get('slurm', {})
    partition = slurm_config.get('partition', '4hours')
    threads = slurm_config.get('threads', 20)
    mem = slurm_config.get('mem_simulation', '10G')
    time_limit = slurm_config.get('time', '4:00:00')
    array_limit = slurm_config.get('array_limit', 10)

    # Conda 环境
    conda_path = config.get('conda_path', '/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh')
    conda_env = get_conda_env(config, 'step2', 'simulation')

    # 程序目录
    programs = config.get('programs', {})
    simulation_dir = programs.get('simulation_dir') or str(get_script_dir())

    # 日志目录
    log_dir = output_base / 'logs' / 'simulation'

    cmd = [
        'sbatch',
        f'--partition={partition}',
        f'--time={time_limit}',
        f'--mem={mem}',
        f'-c {threads}',
        f'--array=1-{task_count}%{array_limit}',
        f'--output={log_dir}/simulation-log-%A-%a.out',
    ]
    exclude = slurm_config.get('exclude', '')
    if exclude:
        cmd.append(f'--exclude={exclude}')
    cmd.append(str(slurm_script))

    # 设置环境变量
    env = os.environ.copy()
    env['FILELIST'] = str(filelist)
    env['CONDA_PATH'] = conda_path
    env['CONDA_ENV'] = conda_env
    env['SIMULATION_DIR'] = simulation_dir
    env['SLURM_LOG_DIR'] = str(output_base / 'logs')
    env['PYTHONUNBUFFERED'] = '1'  # Python 输出实时刷新到日志

    print(f"任务数: {task_count}")
    print(f"SLURM 命令:")
    print(f"  {' '.join(cmd)}")
    print(f"环境变量:")
    print(f"  FILELIST={filelist}")
    print(f"  CONDA_PATH={conda_path}")
    print(f"  CONDA_ENV={conda_env}")
    print(f"  SIMULATION_DIR={simulation_dir}")
    print(f"  SLURM_LOG_DIR={output_base / 'logs'}")

    if dry_run:
        print("\n[dry-run] 未提交任务")
        return None

    # 创建日志目录
    log_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    print(result.stdout.rstrip())
    job_id = parse_slurm_job_id(result.stdout)
    print(f"\nSLURM 任务已提交 (job ID: {job_id})")
    return job_id


def step3_merge(config, dry_run=False, **kwargs):
    """
    Step 3: 合并 reads (SLURM)
    """
    print("=" * 60)
    print("Step 3: 合并 reads (SLURM)")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'
    filelist = filelists_dir / 'merge_filelist'

    if not filelist.exists():
        print(f"ERROR: filelist 不存在: {filelist}")
        print("请先运行 Step 1")
        return None

    # 计算任务数
    with open(filelist, 'r') as f:
        task_count = sum(1 for line in f if line.strip())

    script_dir = get_script_dir()
    slurm_script = script_dir / 'scripts' / 'merge_with_slurm.sh'

    # SLURM 参数
    slurm_config = config.get('slurm', {})
    partition = slurm_config.get('partition', '4hours')
    threads = slurm_config.get('threads', 36)
    mem = slurm_config.get('mem_merge', '40G')
    time_limit = slurm_config.get('time', '4:00:00')
    array_limit = slurm_config.get('array_limit', 10)

    # 参考基因组
    ref_template = config['reference']['template']

    # 日志目录
    log_dir = output_base / 'logs' / 'merge'

    dependency = kwargs.get('dependency')
    cmd = [
        'sbatch',
        f'--partition={partition}',
        f'--time={time_limit}',
        f'--mem={mem}',
        f'-c {threads}',
        f'--array=1-{task_count}%{array_limit}',
        f'--output={log_dir}/merge-log-%A-%a.out',
    ]
    exclude = slurm_config.get('exclude', '')
    if exclude:
        cmd.append(f'--exclude={exclude}')
    if dependency and not dry_run:
        cmd.append(f'--dependency=afterok:{dependency}')
    cmd.append(str(slurm_script))

    env = os.environ.copy()
    env['FILELIST'] = str(filelist)
    env['REF_TEMPLATE'] = ref_template
    env['SLURM_LOG_DIR'] = str(output_base / 'logs')

    print(f"任务数: {task_count}")
    print(f"SLURM 命令:")
    print(f"  {' '.join(cmd)}")

    if dry_run:
        print("\n[dry-run] 未提交任务")
        return None

    log_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    print(result.stdout.rstrip())
    job_id = parse_slurm_job_id(result.stdout)
    print(f"\nSLURM 任务已提交 (job ID: {job_id})")
    return job_id


def step4_minimap2(config, dry_run=False, **kwargs):
    """
    Step 4: minimap2 比对 (SLURM)
    """
    print("=" * 60)
    print("Step 4: minimap2 比对 (SLURM)")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'
    filelist = filelists_dir / 'minimap2_filelist'

    if not filelist.exists():
        print(f"ERROR: filelist 不存在: {filelist}")
        print("请先运行 Step 1")
        return None

    # 计算任务数
    with open(filelist, 'r') as f:
        task_count = sum(1 for line in f if line.strip())

    script_dir = get_script_dir()
    slurm_script = script_dir / 'scripts' / 'minimap2_with_slurm.sh'

    # SLURM 参数
    slurm_config = config.get('slurm', {})
    partition = slurm_config.get('partition', '4hours')
    threads = slurm_config.get('threads', 36)
    mem = slurm_config.get('mem_minimap2', '60G')
    time_limit = slurm_config.get('time', '4:00:00')
    array_limit = slurm_config.get('array_limit', 10)

    # Conda 环境
    conda_path = config.get('conda_path', '/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh')
    conda_env = get_conda_env(config, 'step4', 'TEMP3')

    # 日志目录
    log_dir = output_base / 'logs' / 'minimap2'

    dependency = kwargs.get('dependency')
    cmd = [
        'sbatch',
        f'--partition={partition}',
        f'--time={time_limit}',
        f'--mem={mem}',
        f'-c {threads}',
        f'--array=1-{task_count}%{array_limit}',
        f'--output={log_dir}/minimap2-log-%A-%a.out',
    ]
    exclude = slurm_config.get('exclude', '')
    if exclude:
        cmd.append(f'--exclude={exclude}')
    if dependency and not dry_run:
        cmd.append(f'--dependency=afterok:{dependency}')
    cmd.append(str(slurm_script))

    env = os.environ.copy()
    env['FILELIST'] = str(filelist)
    env['CONDA_PATH'] = conda_path
    env['CONDA_ENV'] = conda_env
    env['THREADS'] = str(threads)
    env['SLURM_LOG_DIR'] = str(output_base / 'logs')

    print(f"任务数: {task_count}")
    print(f"SLURM 命令:")
    print(f"  {' '.join(cmd)}")

    if dry_run:
        print("\n[dry-run] 未提交任务")
        return None

    log_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    print(result.stdout.rstrip())
    job_id = parse_slurm_job_id(result.stdout)
    print(f"\nSLURM 任务已提交 (job ID: {job_id})")
    return job_id


def step5_label(config, dry_run=False, **kwargs):
    """
    Step 5: 标注 TP/FP (SLURM)
    """
    print("=" * 60)
    print("Step 5: 标注 TP/FP (SLURM)")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'
    filelist = filelists_dir / 'label_filelist'

    if not filelist.exists():
        print(f"ERROR: filelist 不存在: {filelist}")
        print("请先运行 Step 1")
        return None

    # 计算任务数
    with open(filelist, 'r') as f:
        task_count = sum(1 for line in f if line.strip())

    script_dir = get_script_dir()
    slurm_script = script_dir / 'scripts' / 'label_with_slurm.sh'

    # SLURM 参数
    slurm_config = config.get('slurm', {})
    partition = slurm_config.get('partition', '4hours')
    threads = slurm_config.get('threads', 24)
    mem = slurm_config.get('mem_label', '30G')
    time_limit = slurm_config.get('time', '4:00:00')
    array_limit = slurm_config.get('array_limit', 10)

    # Conda 环境
    conda_path = config.get('conda_path', '/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh')
    conda_env = get_conda_env(config, 'step5', 'TEMP3')

    # 参考基因组
    ref_template = config['reference']['template']

    # 程序目录
    programs = config.get('programs', {})
    script_dir = get_script_dir()
    simulation_dir = programs.get('simulation_dir') or str(script_dir)
    label_dir = programs.get('label_dir') or str(Path(simulation_dir) / 'label')

    # 日志目录
    log_dir = output_base / 'logs' / 'label'

    dependency = kwargs.get('dependency')
    cmd = [
        'sbatch',
        f'--partition={partition}',
        f'--time={time_limit}',
        f'--mem={mem}',
        f'-c {threads}',
        f'--array=1-{task_count}%{array_limit}',
        f'--output={log_dir}/label-log-%A-%a.out',
    ]
    exclude = slurm_config.get('exclude', '')
    if exclude:
        cmd.append(f'--exclude={exclude}')
    if dependency and not dry_run:
        cmd.append(f'--dependency=afterok:{dependency}')
    cmd.append(str(slurm_script))

    env = os.environ.copy()
    env['FILELIST'] = str(filelist)
    env['CONDA_PATH'] = conda_path
    env['CONDA_ENV'] = conda_env
    env['REF_TEMPLATE'] = ref_template
    env['LABEL_DIR'] = label_dir
    env['TEMP3_DIR'] = programs.get('temp3_dir') or str(Path(simulation_dir).parent / 'TEMP3')

    ref = config.get('reference', {})
    env['BLACKLIST'] = ref.get('blacklist', '')
    env['GERM_MODEL'] = ref.get('germ_model', '')
    env['SOMA_MODEL'] = ref.get('soma_model', '')
    env['SUBSIZE'] = str(config['simulation']['sub_population_size'])
    env['SLURM_LOG_DIR'] = str(output_base / 'logs')

    print(f"任务数: {task_count}")
    print(f"SLURM 命令:")
    print(f"  {' '.join(cmd)}")

    if dry_run:
        print("\n[dry-run] 未提交任务")
        return None

    log_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    print(result.stdout.rstrip())
    job_id = parse_slurm_job_id(result.stdout)
    print(f"\nSLURM 任务已提交 (job ID: {job_id})")
    return job_id


def step6_classify(config, dry_run=False, **kwargs):
    """
    Step 6: 分类 germline/somatic

    本地运行 classify_germline_and_somatic.sh
    """
    print("=" * 60)
    print("Step 6: 分类 germline/somatic")
    print("=" * 60)

    output_base = Path(config['output_dir'])
    filelists_dir = output_base / 'filelists'
    filelist = filelists_dir / 'label_filelist'

    if not filelist.exists():
        print(f"ERROR: filelist 不存在: {filelist}")
        print("请先运行 Step 1")
        return

    # 计算任务数
    with open(filelist, 'r') as f:
        task_count = sum(1 for line in f if line.strip())

    script_dir = get_script_dir()
    classify_script = script_dir / 'label' / 'classify_germline_somatic.sh'

    classification_config = config.get('classification', {})
    somatic_freq = classification_config.get('somatic_freq', 0.01)

    print(f"处理 {task_count} 个数据集")

    if dry_run:
        print("[dry-run] 未执行分类")
        return

    wait_job_id = kwargs.get('wait_job_id')
    if wait_job_id:
        wait_for_slurm_job(wait_job_id)

    # 执行分类脚本
    env = os.environ.copy()
    env['FILELIST'] = str(filelist)
    env['SOMATIC_FREQ'] = str(somatic_freq)

    subprocess.run(['bash', str(classify_script)], env=env)
    print("\nStep 6 完成!")


def step7_ml_training(config, dry_run=False, **kwargs):
    """
    Step 7: ML 模型训练

    使用 Step 6 输出的标注数据（TP_clt_G.txt, TP_clt_S.txt, FP_clt.txt）
    训练 AutoGluon 分类模型（germline 和/或 somatic）。

    配置参数来自 config 中的 ml_training 段。
    调用 model_training/pipeline/run_pipeline.py 执行。
    """
    print("=" * 60)
    print("Step 7: ML 模型训练")
    print("=" * 60)

    ml_config = config.get('ml_training', {})
    if not ml_config.get('enabled', True):
        print("ML training 已禁用 (ml_training.enabled = false)，跳过")
        return

    script_dir = get_script_dir()
    pipeline_script = script_dir.parent / 'model_training' / 'pipeline' / 'run_pipeline.py'

    if not pipeline_script.exists():
        print(f"ERROR: pipeline 脚本不存在: {pipeline_script}")
        return

    # 确定 mode
    modes = ml_config.get('modes', {})
    has_germline = modes.get('germline', {}).get('enabled', True)
    has_somatic = modes.get('somatic', {}).get('enabled', True)

    if has_germline and has_somatic:
        mode_arg = 'all'
    elif has_germline:
        mode_arg = 'germline'
    elif has_somatic:
        mode_arg = 'somatic'
    else:
        print("WARNING: 没有启用的 mode")
        return

    # 设置 conda 环境 —— 使用目标环境的 Python 解释器完整路径
    conda_path = config.get('conda_path', '/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh')
    conda_env = get_conda_env(config, 'step7', 'TEMP3')

    # 从 conda.sh 路径推导目标环境的 Python 解释器路径
    #   conda.sh -> <conda_root>/etc/profile.d/conda.sh
    #   python    -> <conda_root>/envs/<env>/bin/python
    conda_root = Path(conda_path).parent.parent.parent
    python_exe = str(conda_root / 'envs' / conda_env / 'bin' / 'python')

    cmd = [
        python_exe, str(pipeline_script),
        '--config', str(config['_config_path']),
        '--mode', mode_arg,
    ]

    ml_stage = kwargs.get('ml_stage')
    if ml_stage:
        cmd += ['--stage', ml_stage]

    print(f"执行命令:")
    print(f"  {' '.join(cmd)}")

    if dry_run:
        print("\n[dry-run] 未执行 ML 训练")
        return

    # 设置 PATH 使子进程能找到目标环境中的工具（如 bedtools）
    env = os.environ.copy()
    env_bin = str(conda_root / 'envs' / conda_env / 'bin')
    env['PATH'] = f"{env_bin}:{env.get('PATH', '')}"

    subprocess.run(cmd, check=True, env=env)
    print("\nStep 7 完成!")


def main():
    parser = argparse.ArgumentParser(description='Simulation Pipeline 主控脚本')
    parser.add_argument('--config', required=True, help='配置文件路径 (YAML)')
    parser.add_argument('--step', type=int, default=None,
                        help='执行到哪一步 (0-7)。不指定则显示所有步骤信息。配合 --start 可以只执行中间某几步')
    parser.add_argument('--start', type=int, default=0,
                        help='从哪一步开始执行 (默认 0)')
    parser.add_argument('--dry-run', action='store_true',
                        help='只显示命令，不实际执行')
    parser.add_argument('--ml-stage', default=None,
                        help='Step 7 的 ML pipeline stage 范围 (e.g., "1", "1-3", "1,3")。不指定则跑全部 stage')

    args = parser.parse_args()

    config = load_config(args.config)
    config['_config_path'] = args.config  #保存原始配置文件路径

    print(f"物种: {config['species']}")
    print(f"输出目录: {config['output_dir']}")
    print(f"数据集总数: {sum(config['datasets'].values())} x {len(config['sequencing']['protocols'])} protocols")
    print()

    steps = [
        ("Step 0", "定义群体基因组 (define_pgdf)", step0_define_pgdf),
        ("Step 1", "生成 filelists", step1_generate_filelists),
        ("Step 2", "运行 simulation (SLURM)", step2_simulation),
        ("Step 3", "合并 reads (SLURM)", step3_merge),
        ("Step 4", "minimap2 比对 (SLURM)", step4_minimap2),
        ("Step 5", "标注 TP/FP (SLURM)", step5_label),
        ("Step 6", "分类 germline/somatic", step6_classify),
        ("Step 7", "ML 模型训练", step7_ml_training),
    ]

    if args.step is None:
        print("可用步骤:")
        for i, (name, desc, _) in enumerate(steps):
            print(f"  {i}: {desc}")
        print()
        print("用法: python run_simulation.py --config <config.yaml> --step <0-7>")
        print("      python run_simulation.py --config <config.yaml> --dry-run")
        return

    if args.step < 0 or args.step > 7:
        print(f"ERROR: step 必须在 0-7 之间")
        return

    if args.start > args.step:
        print(f"ERROR: --start ({args.start}) 不能大于 --step ({args.step})")
        return

    # 执行指定步骤，自动传递 SLURM 作业依赖
    last_job_id = None
    for i in range(args.start, args.step + 1):
        name, desc, func = steps[i]

        kwargs = {'dry_run': args.dry_run}
        if i == 7 and args.ml_stage:
            kwargs['ml_stage'] = args.ml_stage
        if i in (3, 4, 5) and last_job_id:
            kwargs['dependency'] = last_job_id
        elif i == 6 and last_job_id:
            kwargs['wait_job_id'] = last_job_id

        job_id = func(config, **kwargs)
        if i in (2, 3, 4, 5):
            last_job_id = job_id


if __name__ == '__main__':
    main()