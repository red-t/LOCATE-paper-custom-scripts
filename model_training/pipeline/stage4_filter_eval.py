"""
Stage 4: Filter test data with blacklist and evaluate all models.

For each test ratio + protocol combination:
  1. Intersect test data with BlackList.bed to separate filtered/unfiltered
  2. Load each trained model
  3. Predict on unfiltered + filtered data combined
  4. Compute metrics (precision/recall/f1/macro-f1)
  5. Write summary
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

from ..pipeline.config import load_ml_config
from ..pipeline.utils import run_bedtools_intersect, clamp_coordinates_in_file
from ..pipeline.evaluate import load_data, compute_metrics
from ..pipeline.defaults import COLUMN_NAMES


def run(config, mode_config, work_dir):
    """Run Stage 4: filter test data and evaluate models."""
    print(f"[Stage 4] Evaluating models for mode='{mode_config.mode_name}'")
    work_dir = Path(work_dir)
    test_dir = work_dir / 'Test'
    blacklist = work_dir / 'BlackList.bed'

    # Discover test ratio files
    # Files are named: test_P_{protocol}_{preset}_{ratio}.txt
    test_ratios = set()
    for fpath in test_dir.glob('test_P_*.txt'):
        parts = fpath.stem.split('_')  # ['test', 'P', 'ccs', 'map-hifi', 'ORG']
        if len(parts) >= 5:
            test_ratios.add(parts[-1])  # ratio name

    if not test_ratios:
        print("  WARNING: No test ratio files found.")
        return

    # Discover trained models
    model_dirs = []
    mode_prefix = 'G' if mode_config.mode_name == 'germline' else 'S'
    for d in work_dir.glob(f'{config.species}_{mode_prefix}_*'):
        if d.is_dir() and (d / 'predictor.pkl').exists():
            model_dirs.append(d)

    if not model_dirs:
        print("  WARNING: No trained models found. Run Stage 3 first.")
        return

    print(f"  Found {len(model_dirs)} models, {len(test_ratios)} test ratios")

    # Filter test data with blacklist
    if blacklist.exists() and blacklist.stat().st_size > 0:
        print("  Filtering test data with blacklist...")
        for fpath in list(test_dir.glob('test_*.txt')):
            clamped = test_dir / f'.tmp_clamped_{fpath.name}'
            filtered_out = test_dir / f'filtered_{fpath.name}'  # -wa (overlap)
            unfiltered_out = test_dir / f'unfilt_{fpath.name}'  # -v (no overlap)

            # Clamp coordinates
            clamp_coordinates_in_file(fpath, clamped)

            # Get filtered records (intersect with blacklist)
            run_bedtools_intersect(clamped, blacklist, filtered_out, invert=False)
            # Get unfiltered records (no blacklist overlap)
            run_bedtools_intersect(clamped, blacklist, unfiltered_out, invert=True)

            clamped.unlink(missing_ok=True)

    # Track if we wrote to summary
    summary_written = False

    # Evaluate each model
    mode_short = 'G' if mode_config.mode_name == 'germline' else 'S'
    summary_path = work_dir / f'{config.species}_{mode_short}_summary_Dedup.txt'
    fout = open(summary_path, 'a')

    for model_dir in sorted(model_dirs):
        # Extract the ratio name from model dir name
        model_name = model_dir.name  # e.g., dm6_G_ORG
        ratio_name = model_name.split('_')[-1]  # e.g., ORG

        # Load model
        try:
            from autogluon.tabular import TabularPredictor
        except ImportError:
            print("  ERROR: autogluon not installed.")
            sys.exit(1)

        print(f"  Loading model: {model_dir.name}")
        predictor = TabularPredictor.load(str(model_dir))

        # Process each test ratio × protocol combination
        for test_ratio in sorted(test_ratios):
            for fpath in sorted(test_dir.glob(f'test_P_*_{test_ratio}.txt')):
                # Extract protocol and preset from the filename
                parts = fpath.stem.split('_')
                # test_P_{protocol}_{preset}_{ratio}
                protocol = parts[2]
                preset = parts[3]

                # Load unfiltered data
                p_unfilt = test_dir / f'unfilt_test_P_{protocol}_{preset}_{test_ratio}.txt'
                n_unfilt = test_dir / f'unfilt_test_N_{protocol}_{preset}_{test_ratio}.txt'
                # Load filtered data
                p_filt = test_dir / f'filtered_test_P_{protocol}_{preset}_{test_ratio}.txt'
                n_filt = test_dir / f'filtered_test_N_{protocol}_{preset}_{test_ratio}.txt'

                # If no blacklist, use original files as unfiltered
                if not p_unfilt.exists():
                    p_unfilt = fpath
                    n_unfilt = fpath.parent / f'test_N_{protocol}_{preset}_{test_ratio}.txt'
                    # No filtered data
                    p_filt = None

                # Load data
                df_test = load_data(p_unfilt, n_unfilt)
                if df_test.empty:
                    continue

                # Load filtered data (if any) — mark as class 0 (caught by blacklist)
                df_filtered = pd.DataFrame()
                if p_filt and p_filt.exists() and p_filt.stat().st_size > 0:
                    df_filtered = load_data(p_filt, n_filt)
                    if not df_filtered.empty:
                        df_filtered['class'] = 0
                        df_filtered['filtered'] = 1
                if not df_filtered.empty:
                    df_test['filtered'] = 0
                    df_combined = pd.concat([df_test, df_filtered], ignore_index=True)
                else:
                    df_combined = df_test

                # Predict
                y_pred = np.asarray(predictor.predict(df_combined).astype(int))
                y_true = np.asarray(df_combined['class'].astype(int))

                # Compute metrics
                metrics = compute_metrics(y_true, y_pred, pos_label=mode_config.pos_label)

                # Write summary
                platform = preset  # e.g., map-hifi, map-ont
                line = f"{ratio_name}\t{test_ratio}\t{protocol}\t{platform}\t"
                line += f"{metrics['precision']:.6f}\t{metrics['recall']:.6f}\t"
                line += f"{metrics['f1']:.6f}\t{metrics['macro_f1']:.6f}"
                print(f"    {line}")
                print(line, file=fout)
                summary_written = True

    fout.close()

    # Clean up temp files
    for f in test_dir.glob('unfilt_*'):
        f.unlink(missing_ok=True)
    for f in test_dir.glob('filtered_*'):
        f.unlink(missing_ok=True)

    if summary_written:
        print(f"  Summary written to: {summary_path}")
    print(f"[Stage 4] Complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage 4: Filter and evaluate')
    parser.add_argument('--config', required=True, help='Pipeline config YAML')
    parser.add_argument('--mode', required=True, choices=['germline', 'somatic'])
    parser.add_argument('--work-dir', default=None)
    args = parser.parse_args()

    cfg = load_ml_config(args.config)
    mode_cfg = cfg.germline if args.mode == 'germline' else cfg.somatic
    work_dir = args.work_dir or cfg.work_dir / args.mode

    run(cfg, mode_cfg, work_dir)
