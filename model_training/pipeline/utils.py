"""
Shared utilities for the ML training pipeline.
"""

import subprocess
import random
from pathlib import Path
from collections import defaultdict


def parse_filelist(filelist_path):
    """Parse a label_filelist.

    Each line format: bam_path repeat_bed gap_bed ref_te

    Returns list of dicts with resolved paths and extracted metadata.
    """
    entries = []
    with open(filelist_path) as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            bam = parts[0]
            workdir = Path(bam).parent
            input_dir = workdir / 'result_no_secondary'
            preset = workdir.name  # e.g., map-hifi, map-ont
            id1 = workdir.parent.name  # e.g., 1_ccs, 42_clr

            entry = {
                'index': i,
                'bam': bam,
                'workdir': workdir,
                'input_dir': input_dir,
                'preset': preset,
                'id1': id1,
                'sample_id': f'{id1}_{preset}',
            }
            entries.append(entry)
    return entries


def extract_protocol(sample_id):
    """Extract sequencing protocol from a sample_id like '1_ccs_map-hifi'."""
    parts = sample_id.split('_')
    if len(parts) >= 2:
        return parts[1]  # ccs, clr, ont
    return 'unknown'


def stratified_train_test_split(entries, test_frac=0.1, random_seed=42):
    """Split entries into train/test by stratified sampling on (depth, protocol).

    Groups entries by (depth, protocol), then selects test_frac from each group.
    Returns (train_entries, test_entries).
    """
    rng = random.Random(random_seed)

    # Infer depth and protocol from each entry
    for entry in entries:
        entry['protocol'] = extract_protocol(entry['sample_id'])
        # Extract depth from id1 (e.g., '1_ccs' -> depth comes from
        # the dataset numbering scheme; we use id1 directly)
        entry['depth'] = entry['id1']

    # Group by protocol for stratified sampling
    by_protocol = defaultdict(list)
    for entry in entries:
        by_protocol[entry['protocol']].append(entry)

    train = []
    test = []
    for protocol, group in by_protocol.items():
        # Sort for reproducible splits
        group.sort(key=lambda e: e['index'])
        n_test = max(1, round(len(group) * test_frac))
        rng.shuffle(group)
        test.extend(group[:n_test])
        train.extend(group[n_test:])

    train.sort(key=lambda e: e['index'])
    test.sort(key=lambda e: e['index'])
    return train, test


def run_bedtools_intersect(input_bed, blacklist_bed, output_bed, invert=False):
    """Run bedtools intersect.

    Args:
        invert: If True, use -v flag (exclude overlaps).
    """
    cmd = ['bedtools', 'intersect', '-a', str(input_bed), '-b', str(blacklist_bed)]
    if invert:
        cmd.append('-v')
    else:
        cmd.append('-wa')

    with open(output_bed, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)


def run_bedtools_merge(input_bed, output_bed, col=4, op='count_distinct'):
    """Run bedtools merge with -c and -o options."""
    cmd = [
        'bedtools', 'merge',
        '-i', str(input_bed),
        '-d', '0',
        '-c', str(col),
        '-o', op,
    ]
    with open(output_bed, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)


def clamp_coordinates_in_file(input_path, output_path):
    """Ensure refStart (column 2) >= 0 in a tab-delimited file."""
    with open(input_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    val = int(parts[1])
                    if val < 0:
                        parts[1] = '0'
                except ValueError:
                    pass
            fout.write('\t'.join(parts) + '\n')
