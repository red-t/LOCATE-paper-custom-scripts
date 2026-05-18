"""
Stage 3: Train AutoGluon models for each ratio.

Scans the Train/ directory for train_P_*.txt / train_N_*.txt pairs,
invokes training.py for each ratio to train a separate model.
"""

import argparse
import subprocess
import sys
from pathlib import Path

from ..pipeline.config import load_ml_config


def run(config, mode_config, work_dir):
    """Run Stage 3: train models for each ratio."""
    print(f"[Stage 3] Training models for mode='{mode_config.mode_name}'")
    work_dir = Path(work_dir)
    train_dir = work_dir / 'Train'

    # Discover ratio files (train_P_*.txt)
    ratio_names = set()
    for fpath in train_dir.glob('train_P_*.txt'):
        name = fpath.stem  # train_P_ORG, train_P_1V1, etc.
        parts = name.split('_', 2)  # ['train', 'P', 'ORG']
        if len(parts) >= 3:
            ratio_names.add(parts[2])

    if not ratio_names:
        # Fall back to configured ratios
        ratio_names = set(r.name for r in mode_config.train_ratios)

    if not ratio_names:
        print("  ERROR: No ratio files found in Train/ directory.")
        sys.exit(1)

    # Find training.py
    script_dir = Path(__file__).resolve().parent.parent  # model_training/
    training_py = script_dir / 'training.py'

    model_dir = work_dir / 'models'
    model_dir.mkdir(parents=True, exist_ok=True)

    drop_cols_str = ','.join(mode_config.drop_columns)

    for ratio_name in sorted(ratio_names):
        pos_file = train_dir / f'train_P_{ratio_name}.txt'
        neg_file = train_dir / f'train_N_{ratio_name}.txt'

        if not pos_file.exists() or pos_file.stat().st_size == 0:
            print(f"    WARNING: {pos_file} is empty. Skipping ratio '{ratio_name}'.")
            continue
        if not neg_file.exists() or neg_file.stat().st_size == 0:
            print(f"    WARNING: {neg_file} is empty. Skipping ratio '{ratio_name}'.")
            continue

        # Model name: {species}_{mode}_{ratio}/
        species = config.species
        mode_short = 'G' if mode_config.mode_name == 'germline' else 'S'
        model_name = f'{species}_{mode_short}_{ratio_name}'
        model_path = work_dir / model_name

        if model_path.exists():
            print(f"    Model '{model_name}' already exists. Skipping.")
            continue

        cmd = [
            'python', str(training_py),
            '--positive', str(pos_file),
            '--negative', str(neg_file),
            '--model_name', str(model_path),
            '--drop-cols', drop_cols_str,
            '--pos-label', str(mode_config.pos_label),
        ]

        print(f"    Training '{model_name}' (P:{pos_file.stat().st_size}B, N:{neg_file.stat().st_size}B)...")
        print(f"      Drop cols: {drop_cols_str}")
        print(f"      Pos label: {mode_config.pos_label}")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"      ERROR: Training failed for '{model_name}'")
            print(f"      {result.stderr[:500]}")
        else:
            print(f"      Model '{model_name}' trained successfully.")

    print(f"[Stage 3] Complete. Models in: {work_dir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage 3: Train models')
    parser.add_argument('--config', required=True, help='Pipeline config YAML')
    parser.add_argument('--mode', required=True, choices=['germline', 'somatic'])
    parser.add_argument('--work-dir', default=None)
    args = parser.parse_args()

    cfg = load_ml_config(args.config)
    mode_cfg = cfg.germline if args.mode == 'germline' else cfg.somatic
    work_dir = args.work_dir or cfg.work_dir / args.mode

    run(cfg, mode_cfg, work_dir)
