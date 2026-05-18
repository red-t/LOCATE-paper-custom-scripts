"""
Stage 2: Generate different P/N ratio datasets by shuffling.

Reads train_P.txt / train_N.txt from Stage 1, auto-detects counts,
then samples according to configured ratios. Also handles test data.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from ..pipeline.utils import extract_protocol
from ..pipeline.config import load_ml_config


def _compute_sample_count(total_pos, total_neg, ratio_def):
    """Determine positive/negative sample count for a ratio.

    If ratio_def has explicit counts, use them.
    For ORG: use all available data.
    For N_V_M (e.g. 1V1, 1V30): bidirectional downsampling to maintain target ratio.

    The target ratio P:N = p_ratio:n_ratio. We compare the actual data ratio
    against the target to decide which side to downsample:

      total_neg * p_ratio >= total_pos * n_ratio  →  negatives are abundant
        → keep all positives, downsample negatives to total_pos * n_ratio / p_ratio
      total_neg * p_ratio < total_pos * n_ratio   →  positives are abundant
        → keep all negatives, downsample positives to total_neg * p_ratio / n_ratio
    """
    if ratio_def.positive is not None:
        n_pos = min(ratio_def.positive, total_pos)
    elif ratio_def.name.upper() == 'ORG':
        n_pos = total_pos
    else:
        # Try to parse ratio from name like "1V1", "1V30"
        try:
            parts = ratio_def.name.upper().split('V')
            if len(parts) == 2:
                p_ratio = int(parts[0])
                n_ratio = int(parts[1])
                if total_neg * p_ratio >= total_pos * n_ratio:
                    n_pos = total_pos
                else:
                    n_pos = total_neg * p_ratio // n_ratio
            else:
                n_pos = total_pos
        except ValueError:
            n_pos = total_pos

    if ratio_def.negative is not None:
        n_neg = min(ratio_def.negative, total_neg)
    elif ratio_def.name.upper() == 'ORG':
        n_neg = total_neg
    else:
        try:
            parts = ratio_def.name.upper().split('V')
            if len(parts) == 2:
                p_ratio = int(parts[0])
                n_ratio = int(parts[1])
                if total_neg * p_ratio >= total_pos * n_ratio:
                    n_neg = total_pos * n_ratio // p_ratio
                else:
                    n_neg = total_neg
            else:
                n_neg = total_neg
        except ValueError:
            n_neg = total_neg

    # Ensure we don't exceed available counts
    n_pos = min(n_pos, total_pos)
    n_neg = min(n_neg, total_neg)
    return n_pos, n_neg


def _shuffle_and_write(df, n, replace, output_path):
    """Sample n rows from df (with or without replacement) and write."""
    if n <= 0 or df.empty:
        pd.DataFrame().to_csv(output_path, sep='\t', header=False, index=False)
        return
    sampled = df.sample(n=n, replace=replace)
    sampled.to_csv(output_path, sep='\t', header=False, index=False)


def run(config, mode_config, work_dir):
    """Run Stage 2: shuffle and generate ratio datasets."""
    print(f"[Stage 2] Generating ratio datasets for mode='{mode_config.mode_name}'")
    work_dir = Path(work_dir)
    train_dir = work_dir / 'Train'
    test_dir = work_dir / 'Test'

    # === Load training data ===
    train_p_path = train_dir / 'train_P.txt'
    train_n_path = train_dir / 'train_N.txt'

    if not train_p_path.exists():
        print(f"  ERROR: {train_p_path} not found. Run Stage 1 first.")
        sys.exit(1)

    df_train_p = pd.read_table(train_p_path, header=None) if train_p_path.stat().st_size > 0 else pd.DataFrame()
    df_train_n = pd.read_table(train_n_path, header=None) if train_n_path.exists() and train_n_path.stat().st_size > 0 else pd.DataFrame()

    total_pos = len(df_train_p)
    total_neg = len(df_train_n)
    print(f"  Available: {total_pos} positive, {total_neg} negative")

    # === Generate training ratios ===
    train_out = work_dir / 'Train'
    train_out.mkdir(parents=True, exist_ok=True)

    for ratio_def in mode_config.train_ratios:
        n_pos, n_neg = _compute_sample_count(total_pos, total_neg, ratio_def)
        replace = mode_config.negative_sampling == 'with_replacement'

        print(f"    Ratio '{ratio_def.name}': sampling P={n_pos}/{total_pos}, N={n_neg}/{total_neg}")
        _shuffle_and_write(df_train_p, n_pos, replace=False,
                           output_path=train_out / f'train_P_{ratio_def.name}.txt')
        _shuffle_and_write(df_train_n, n_neg, replace=replace,
                           output_path=train_out / f'train_N_{ratio_def.name}.txt')

    # === Generate test ratios ===
    test_out = work_dir / 'Test'
    test_out.mkdir(parents=True, exist_ok=True)

    # Group test files by protocol
    from collections import defaultdict
    test_files = defaultdict(lambda: {'p': None, 'n': None})
    for fpath in test_dir.glob('test_*.txt'):
        name = fpath.stem  # e.g. test_P_ccs_map-hifi
        parts = name.split('_', 3)  # ['test', 'P', 'ccs', 'map-hifi']
        if len(parts) >= 4:
            pn = parts[1]   # P or N
            protocol = parts[2]
            preset = parts[3]
            key = f'{protocol}_{preset}'
            if pn == 'P':
                test_files[key]['p'] = fpath
            else:
                test_files[key]['n'] = fpath

    for key, files in test_files.items():
        df_p = pd.read_table(files['p'], header=None) if files['p'] and files['p'].stat().st_size > 0 else pd.DataFrame()
        df_n = pd.read_table(files['n'], header=None) if files['n'] and files['n'].stat().st_size > 0 else pd.DataFrame()
        total_p = len(df_p)
        total_n = len(df_n)

        for ratio_def in mode_config.test_ratios:
            n_pos, n_neg = _compute_sample_count(total_p, total_n, ratio_def)
            replace = mode_config.negative_sampling == 'with_replacement'
            ratio_name = ratio_def.name

            print(f"    Test '{key}' ratio '{ratio_name}': P={n_pos}, N={n_neg}")
            _shuffle_and_write(df_p, n_pos, replace=False,
                               output_path=test_out / f'test_P_{key}_{ratio_name}.txt')
            _shuffle_and_write(df_n, n_neg, replace=replace,
                               output_path=test_out / f'test_N_{key}_{ratio_name}.txt')

    print(f"[Stage 2] Complete. Output in: {work_dir / 'Train'}, {work_dir / 'Test'}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage 2: Generate ratio datasets')
    parser.add_argument('--config', required=True, help='Pipeline config YAML')
    parser.add_argument('--mode', required=True, choices=['germline', 'somatic'])
    parser.add_argument('--work-dir', default=None)
    args = parser.parse_args()

    cfg = load_ml_config(args.config)
    mode_cfg = cfg.germline if args.mode == 'germline' else cfg.somatic
    work_dir = args.work_dir or cfg.work_dir / args.mode

    run(cfg, mode_cfg, work_dir)
