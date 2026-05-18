"""
Stage 1: Data Preprocessing (merge → split → train/test → blacklist → dedup).

Reads labeled TP/FP data from each sample's result_no_secondary/ directory,
filters by mode-specific criteria, and prepares training/test datasets.
"""

import argparse
import sys
import shutil
from pathlib import Path

import pandas as pd

from ..pipeline.defaults import COLUMN_NAMES
from ..pipeline.utils import (
    parse_filelist,
    run_bedtools_intersect,
    run_bedtools_merge,
    clamp_coordinates_in_file,
    stratified_train_test_split,
)
from ..pipeline.config import load_ml_config


def run(config, mode_config, work_dir):
    """Run Stage 1 preprocessing.

    Steps:
      1. Merge — read all samples, filter, insert sample_id, concatenate
      2. Split — per-sample files (not strictly needed, kept for compatibility)
      3. Train/test split — stratified by (depth, protocol)
      4. Blacklist construction — from FP regions across all samples
      5. Deduplication — remove training records in blacklist regions
    """
    print(f"[Stage 1] Preprocessing data for mode='{mode_config.mode_name}'")
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Resolve cltType and teAlignedFrac column indices (0-based) in raw 33-col input
    # Current format: cltType is col 11 → index 10, teAlignedFrac is col 30 → index 29
    clt_type_idx = 10     # $11 in 1-based
    score_idx = 29        # $30 in 1-based

    # Parse filter expression
    filter_op = mode_config.positive_filter.cltType_expr  # e.g. "==0", ">0"
    score_threshold = mode_config.positive_filter.score_threshold

    def _filter_clt_type(val):
        if filter_op == "==0":
            return val == 0
        elif filter_op == ">0":
            return val > 0
        elif filter_op == ">=0":
            return val >= 0
        elif filter_op == "==1":
            return val == 1
        elif filter_op == ">1":
            return val > 1
        # Add more as needed
        return True

    # === Step 1: Merge all samples ===
    print("  Merging all samples...")
    entries = parse_filelist(config.filelist_path)

    all_positive = []
    all_negative = []
    all_fp_blacklist = []

    for entry in entries:
        sample_id = entry['sample_id']
        input_dir = entry['input_dir']

        # Read positive sources
        for src in mode_config.positive_sources:
            src_path = input_dir / src
            if not src_path.exists():
                continue
            df = pd.read_table(src_path, header=None, sep='\s+')
            if df.shape[1] < 33:
                print(f"    WARNING: {src_path} has {df.shape[1]} columns, expected 33. Skipping.")
                continue
            # Filter by cltType and score
            if len(df) > 0:
                mask = df.apply(
                    lambda row: _filter_clt_type(row.iloc[clt_type_idx])
                                and float(row.iloc[score_idx]) >= score_threshold,
                    axis=1,
                )
                df_filtered = df[mask].copy()
                # Insert sample_id at position 6 (0-based, becomes col 7)
                df_filtered.insert(6, 'sample_id', sample_id)
                all_positive.append(df_filtered)

        # Read negative data (FP)
        fp_path = input_dir / 'FP_clt.txt'
        if fp_path.exists():
            df_fp = pd.read_table(fp_path, header=None, sep='\s+')
            if df_fp.shape[1] >= 33:
                df_fp.insert(6, 'sample_id', sample_id)
                all_negative.append(df_fp)

                # For blacklist: collect FP with coordinate info
                # Use the same filter as positive to get relevant FP
                if len(df_fp) > 0:
                    fp_bl = df_fp.iloc[:, [0, 1, 2, 4, 5]].copy()  # chrom, st, ed, freq, strand
                    fp_bl.columns = range(5)
                    # Add sample_id to last column (for merge counting)
                    fp_bl[5] = sample_id
                    all_fp_blacklist.append(fp_bl)

    if not all_positive:
        print("  ERROR: No positive data found. Check your data paths and filters.")
        sys.exit(1)

    # Concatenate all positive and negative
    df_pos = pd.concat(all_positive, ignore_index=True)
    df_pos.columns = range(34)  # 33 original + 1 inserted id
    df_neg = pd.concat(all_negative, ignore_index=True) if all_negative else pd.DataFrame()
    if not df_neg.empty:
        df_neg.columns = range(34)

    # Assign standard column names
    df_pos.columns = COLUMN_NAMES
    if not df_neg.empty:
        df_neg.columns = COLUMN_NAMES

    print(f"    Positive records: {len(df_pos)}")
    print(f"    Negative records: {len(df_neg) if not df_neg.empty else 0}")

    # Ensure refStart >= 0
    df_pos.loc[df_pos['refStart'] < 0, 'refStart'] = 0
    if not df_neg.empty:
        df_neg.loc[df_neg['refStart'] < 0, 'refStart'] = 0

    # === Step 2: Train/Test split ===
    print("  Splitting train/test (90/10 stratified)...")
    train_entries, test_entries = stratified_train_test_split(
        entries, test_frac=config.test_frac, random_seed=config.random_seed,
    )
    train_ids = set(e['sample_id'] for e in train_entries)
    test_ids = set(e['sample_id'] for e in test_entries)

    print(f"    Training samples: {len(train_ids)}")
    print(f"    Testing samples: {len(test_ids)}")

    # Split positive and negative data
    df_pos_train = df_pos[df_pos['sample_id'].isin(train_ids)]
    df_pos_test = df_pos[df_pos['sample_id'].isin(test_ids)]
    df_neg_train = df_neg[df_neg['sample_id'].isin(train_ids)] if not df_neg.empty else pd.DataFrame()
    df_neg_test = df_neg[df_neg['sample_id'].isin(test_ids)] if not df_neg.empty else pd.DataFrame()

    # Write train files
    train_dir = work_dir / 'Train'
    train_dir.mkdir(parents=True, exist_ok=True)

    df_pos_train.to_csv(train_dir / 'train_P.txt', sep='\t', header=False, index=False)
    if not df_neg_train.empty:
        df_neg_train.to_csv(train_dir / 'train_N.txt', sep='\t', header=False, index=False)
    else:
        # Create empty file
        (train_dir / 'train_N.txt').touch()

    print(f"    Train positive: {len(df_pos_train)}")
    print(f"    Train negative: {len(df_neg_train) if not df_neg_train.empty else 0}")

    # Write test files — grouped by protocol
    test_dir = work_dir / 'Test'
    test_dir.mkdir(parents=True, exist_ok=True)

    for entry in test_entries:
        sid = entry['sample_id']
        protocol = entry['protocol']
        preset = entry['preset']
        suffix = f"{protocol}_{preset}"

        p_test = df_pos_test[df_pos_test['sample_id'] == sid]
        n_test = df_neg_test[df_neg_test['sample_id'] == sid] if not df_neg_test.empty else pd.DataFrame()

        if not p_test.empty:
            p_test.to_csv(test_dir / f'test_P_{suffix}.txt', sep='\t', header=False, index=False, mode='a')
        if not n_test.empty:
            n_test.to_csv(test_dir / f'test_N_{suffix}.txt', sep='\t', header=False, index=False, mode='a')

    # === Step 3: Blacklist construction ===
    print("  Constructing blacklist...")
    if all_fp_blacklist:
        # Concatenate all FP records for blacklist
        df_fp_all = pd.concat(all_fp_blacklist, ignore_index=True)
        df_fp_all.columns = ['chrom', 'st', 'ed', 'freq', 'strand', 'sample_id']
        df_fp_all.loc[df_fp_all['st'] < 0, 'st'] = 0

        # Write temporary files for bedtools
        tmp_dir = work_dir / '.tmp'
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Split into G (cltType==0) and S (cltType>0) FP
        # Note: Here we use positive_examples' cltType to determine if a FP
        # belongs to germline or somatic
        # Actually, looking at the original preprocess.sh, BOTH G and S blacklists
        # are built from FP_clt.txt, differentiated by the same cltType filter
        # that's applied to positives. Let me re-read the original code...

        # Original preprocess.sh builds blacklist by iterating filelist again:
        # awk -v id ... '($11==0 && $32>=0.8){print ...}' FP_clt.txt >> tmp_black_1_G.bed
        # awk -v id ... '($11>0 && $32>=1){print ...}' FP_clt.txt >> tmp_black_1_S.bed
        # So it applies the SAME filter on FP data.
        # For our refactored code, we need separate passes for G and S blacklist.

        # Since we only build one blacklist at a time (one mode per pass),
        # we use the mode's own filter on FP data.
        # Actually the original code always builds BOTH G and S blacklist regardless
        # of which mode we're in. Let me build only the relevant one for now,
        # and the second pass can reuse the first one's output.

        # Re-read FP data with original cltType column to split
        bl_fp_records = []
        for entry in entries:
            input_dir = entry['input_dir']
            fp_path = input_dir / 'FP_clt.txt'
            if not fp_path.exists():
                continue
            df = pd.read_table(fp_path, header=None, sep='\s+')
            if df.shape[1] < 33:
                continue
            sample_id = entry['sample_id']
            for _, row in df.iterrows():
                try:
                    clt_val = int(row.iloc[clt_type_idx])
                    score_val = float(row.iloc[score_idx])
                except (ValueError, TypeError):
                    continue
                st = max(0, int(row.iloc[1]))
                ed = int(row.iloc[2])
                bl_fp_records.append({
                    'chrom': row.iloc[0],
                    'st': st,
                    'ed': ed,
                    'sample_id': sample_id,
                    'cltType': clt_val,
                    'score': score_val,
                })

        if bl_fp_records:
            df_bl = pd.DataFrame(bl_fp_records)

            # Always build both G and S blacklists regardless of current mode.
            # G filter: cltType == 0 AND teAlignedFrac >= 0.8
            # S filter: cltType > 0 AND teAlignedFrac >= 1
            # These thresholds match the original pipeline behavior.
            df_bl_G = df_bl[df_bl.apply(lambda r: r['cltType'] == 0 and r['score'] >= 0.8, axis=1)]
            df_bl_S = df_bl[df_bl.apply(lambda r: r['cltType'] > 0 and r['score'] >= 1, axis=1)]

            # Write BED files
            for bl_name, bl_df in [('G', df_bl_G), ('S', df_bl_S)]:
                if bl_df.empty:
                    continue
                bed_raw = tmp_dir / f'tmp_black_1_{bl_name}.bed'
                bl_df[['chrom', 'st', 'ed', 'sample_id']].to_csv(
                    bed_raw, sep='\t', header=False, index=False)

                bed_sorted = tmp_dir / f'tmp_black_2_{bl_name}.bed'
                bed_sorted_fname = str(bed_sorted)
                bed_merged = tmp_dir / f'tmp_black_3_{bl_name}.bed'
                bed_filtered = tmp_dir / f'tmp_black_4_{bl_name}.bed'

                # Sort
                df_sorted = pd.read_table(bed_raw, header=None)
                df_sorted = df_sorted.sort_values([0, 1])
                df_sorted.to_csv(bed_sorted, sep='\t', header=False, index=False)

                # Merge
                run_bedtools_merge(bed_sorted, bed_merged, col=4, op='count_distinct')

                # Filter: count >= min_overlap AND length < max_length
                df_merged = pd.read_table(bed_merged, header=None)
                if len(df_merged) > 0:
                    df_merged['length'] = df_merged[2] - df_merged[1]
                    df_filtered = df_merged[
                        (df_merged[3] >= config.blacklist_min_overlap)
                        & (df_merged['length'] < config.blacklist_max_length)
                    ]
                    df_filtered = df_filtered.drop(columns=['length'])

                    if not df_filtered.empty:
                        df_filtered.to_csv(bed_filtered, sep='\t', header=False, index=False)
                        bl_out = work_dir / f'BlackList_{bl_name}.bed'
                        shutil.copy2(bed_filtered, bl_out)

            # Combine G + S into BlackList.bed
            bl_g = work_dir / 'BlackList_G.bed'
            bl_s = work_dir / 'BlackList_S.bed'
            bl_all = work_dir / 'BlackList.bed'

            dfs_to_combine = []
            if bl_g.exists() and bl_g.stat().st_size > 0:
                dfs_to_combine.append(pd.read_table(bl_g, header=None))
            if bl_s.exists() and bl_s.stat().st_size > 0:
                dfs_to_combine.append(pd.read_table(bl_s, header=None))

            if dfs_to_combine:
                df_combined = pd.concat(dfs_to_combine, ignore_index=True)
                df_combined = df_combined.sort_values([0, 1])
                df_combined.to_csv(bl_all, sep='\t', header=False, index=False)

            # Clean up tmp
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # === Step 4: Deduplication ===
    print("  Deduplicating training data...")
    dedup_blacklist = work_dir / mode_config.dedup_blacklist
    train_p = train_dir / 'train_P.txt'
    train_n = train_dir / 'train_N.txt'

    if dedup_blacklist.exists() and dedup_blacklist.stat().st_size > 0:
        for fpath in [train_p, train_n]:
            if not fpath.exists() or fpath.stat().st_size == 0:
                continue
            # Clamp coordinates
            clamped = train_dir / '.tmp_clamped'
            clamp_coordinates_in_file(fpath, clamped)
            # Remove blacklist overlaps
            tmp_out = train_dir / '.tmp_dedup'
            run_bedtools_intersect(clamped, dedup_blacklist, tmp_out, invert=True)
            shutil.move(tmp_out, fpath)
            clamped.unlink(missing_ok=True)
    elif not dedup_blacklist.exists():
        print(f"    WARNING: Blacklist {dedup_blacklist} not found. Skipping dedup.")

    # Print final counts
    n_pos = len(pd.read_table(train_p, header=None)) if train_p.exists() and train_p.stat().st_size > 0 else 0
    n_neg = len(pd.read_table(train_n, header=None)) if train_n.exists() and train_n.stat().st_size > 0 else 0
    print(f"    After dedup - Positive: {n_pos}, Negative: {n_neg}")

    print(f"[Stage 1] Complete. Output in: {work_dir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage 1: Preprocess')
    parser.add_argument('--config', required=True, help='Pipeline config YAML')
    parser.add_argument('--mode', required=True, choices=['germline', 'somatic'],
                        help='Mode: germline or somatic')
    parser.add_argument('--work-dir', default=None, help='Override working directory')
    args = parser.parse_args()

    cfg = load_ml_config(args.config)
    mode_cfg = cfg.germline if args.mode == 'germline' else cfg.somatic
    work_dir = args.work_dir or cfg.work_dir / args.mode

    run(cfg, mode_cfg, work_dir)
