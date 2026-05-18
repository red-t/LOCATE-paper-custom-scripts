"""
ML Training Pipeline — Standalone orchestrator.

Runs a subset or all of the ML training stages for a given mode.

Usage:
    # Run all stages for germline mode
    python run_pipeline.py --config species_configs/dm6.yaml --mode germline

    # Run all modes
    python run_pipeline.py --config species_configs/dm6.yaml --mode all

    # Run specific stages
    python run_pipeline.py --config species_configs/dm6.yaml --mode germline --stage 1-2

    # Override working directory
    python run_pipeline.py --config species_configs/dm6.yaml --mode all --work-dir /path/to/output
"""

import argparse
import sys
from pathlib import Path

# Add project root to Python path
_project_root = Path(__file__).resolve().parent.parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

from model_training.pipeline.config import load_ml_config
from model_training.pipeline import stage1_preprocess
from model_training.pipeline import stage2_shuffle
from model_training.pipeline import stage3_train
from model_training.pipeline import stage4_filter_eval


# Stage registry: (name, module.run function)
STAGES = [
    ('stage1_preprocess', stage1_preprocess.run),
    ('stage2_shuffle', stage2_shuffle.run),
    ('stage3_train', stage3_train.run),
    ('stage4_filter_eval', stage4_filter_eval.run),
]


def parse_stage_range(stage_arg):
    """Parse --stage argument (e.g., '1', '1-3', '1,3'). Returns list of 0-based indices."""
    if stage_arg is None:
        return list(range(len(STAGES)))

    indices = set()
    for part in stage_arg.split(','):
        part = part.strip()
        if '-' in part:
            a, b = part.split('-', 1)
            start = int(a.strip()) - 1
            end = int(b.strip()) - 1
            indices.update(range(start, end + 1))
        else:
            indices.add(int(part.strip()) - 1)

    return sorted(i for i in indices if 0 <= i < len(STAGES))


def main():
    parser = argparse.ArgumentParser(description='ML Training Pipeline Orchestrator')
    parser.add_argument('--config', required=True, help='Path to species YAML config')
    parser.add_argument('--mode', default='all', choices=['germline', 'somatic', 'all'],
                        help='Which mode(s) to run')
    parser.add_argument('--stage', default=None,
                        help='Stage range (e.g., "1", "1-3", "1,3"). Default: all stages')
    parser.add_argument('--work-dir', default=None,
                        help='Override working directory')
    args = parser.parse_args()

    # Load config
    config = load_ml_config(args.config)
    if not config.enabled:
        print("ML training is disabled in config.")
        return

    # Determine which modes to run
    if args.mode == 'all':
        modes_to_run = []
        if config.germline.enabled:
            modes_to_run.append(('germline', config.germline))
        if config.somatic.enabled:
            modes_to_run.append(('somatic', config.somatic))
    else:
        mode_cfg = config.germline if args.mode == 'germline' else config.somatic
        if not mode_cfg.enabled:
            print(f"Mode '{args.mode}' is disabled in config.")
            return
        modes_to_run = [(args.mode, mode_cfg)]

    if not modes_to_run:
        print("No enabled modes to run.")
        return

    # Parse stage range
    stage_indices = parse_stage_range(args.stage)
    stage_names = [STAGES[i][0] for i in stage_indices]

    print("=" * 60)
    print(f"ML Training Pipeline")
    print(f"  Config: {args.config}")
    print(f"  Species: {config.species}")
    print(f"  Modes: {[m[0] for m in modes_to_run]}")
    print(f"  Stages: {stage_names}")
    print(f"  Work dirs: config.work_dir / {{mode}}")
    print("=" * 60)

    # Run each mode
    for mode_name, mode_cfg in modes_to_run:
        print(f"\n{'=' * 60}")
        print(f"Mode: {mode_name}")
        print(f"{'=' * 60}")

        work_dir = args.work_dir or str(config.work_dir / mode_name)

        for stage_idx in stage_indices:
            stage_name, stage_func = STAGES[stage_idx]
            print(f"\n--- Running {stage_name} for {mode_name} ---")
            stage_func(config, mode_cfg, work_dir)

    print("\n=== Pipeline complete! ===")


if __name__ == '__main__':
    main()
