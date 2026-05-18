"""
Config parsing for the ML training pipeline (Step 7).

Reads the ml_training section from a species YAML config and
builds MLConfig / ModeConfig dataclasses with validated defaults.
"""

import yaml
from dataclasses import dataclass, field
from typing import List, Optional
from pathlib import Path

from .defaults import DEFAULT_DROP_COLUMNS, COLUMN_NAMES


@dataclass
class PositiveFilter:
    """Filter criteria for positive examples."""
    cltType_expr: str = "==0"       # e.g. "==0" for germline, ">0" for somatic
    score_threshold: float = 0.8    # teAlignedFrac threshold


@dataclass
class RatioDef:
    """A named ratio definition for sampling."""
    name: str                       # e.g. "ORG", "1V1", "1V30"
    positive: Optional[int] = None  # None = auto-detect
    negative: Optional[int] = None  # None = auto-detect


@dataclass
class ModeConfig:
    """Configuration for a single mode (germline or somatic)."""
    enabled: bool = True
    mode_name: str = "germline"
    pos_label: int = 1
    positive_sources: List[str] = field(default_factory=lambda: ['TP_clt_G.txt'])
    positive_filter: PositiveFilter = field(default_factory=PositiveFilter)
    dedup_blacklist: str = "BlackList_G.bed"
    drop_columns: List[str] = field(default_factory=list)
    train_ratios: List[RatioDef] = field(default_factory=list)
    test_ratios: List[RatioDef] = field(default_factory=list)
    negative_sampling: str = "without_replacement"  # or "with_replacement"


@dataclass
class MLConfig:
    """Top-level config for the ML training pipeline."""
    species: str = "dm6"
    enabled: bool = True
    output_dir: Optional[Path] = None
    work_dir: Optional[Path] = None
    filelist_path: Optional[Path] = None
    test_frac: float = 0.1
    random_seed: int = 42

    # Blacklist
    blacklist_min_overlap: int = 10
    blacklist_max_length: int = 50

    # Modes
    germline: ModeConfig = field(default_factory=ModeConfig)
    somatic: ModeConfig = field(default_factory=ModeConfig)

    # Column positions in raw input (34-col format)
    column_schema_version: str = "v4"  # for future format compatibility


def _parse_ratio_list(ratio_list_raw):
    """Parse a list of ratio dicts from YAML into RatioDef objects."""
    ratios = []
    for item in ratio_list_raw:
        if isinstance(item, str):
            ratios.append(RatioDef(name=item))
        elif isinstance(item, dict):
            ratios.append(RatioDef(
                name=item.get('name', item.get('ratio', 'unknown')),
                positive=item.get('positive'),
                negative=item.get('negative'),
            ))
    return ratios


def build_mode_config(mode_name, mode_raw, default_drop_cols):
    """Build a ModeConfig from raw YAML dict + defaults."""
    filter_raw = mode_raw.get('positive_filter', {})

    # Default ratios per mode
    default_train_ratios = {
        'germline': [{'name': 'ORG'}, {'name': '1V1'}],
        'somatic': [{'name': '1V30'}],
    }
    default_test_ratios = {
        'germline': [{'name': 'ORG'}, {'name': '1V1'}],
        'somatic': [{'name': '1V30'}],
    }
    default_neg_sampling = {
        'germline': 'without_replacement',
        'somatic': 'with_replacement',
    }

    ratios_raw = mode_raw.get('train_ratios', default_train_ratios.get(mode_name, []))
    test_ratios_raw = mode_raw.get('test_ratios', default_test_ratios.get(mode_name, []))

    return ModeConfig(
        enabled=mode_raw.get('enabled', True),
        mode_name=mode_name,
        pos_label=mode_raw.get('pos_label', 1 if mode_name == 'germline' else 2),
        positive_sources=mode_raw.get('positive_sources',
                                       ['TP_clt_G.txt'] if mode_name == 'germline'
                                       else ['TP_clt_G.txt', 'TP_clt_S.txt']),
        positive_filter=PositiveFilter(
            cltType_expr=filter_raw.get('cltType_expr', "==0" if mode_name == 'germline' else ">0"),
            score_threshold=filter_raw.get('score_threshold', 0.8 if mode_name == 'germline' else 1),
        ),
        dedup_blacklist=mode_raw.get('dedup_blacklist',
                                      'BlackList_G.bed' if mode_name == 'germline'
                                      else 'BlackList_S.bed'),
        drop_columns=mode_raw.get('drop_columns', default_drop_cols),
        train_ratios=_parse_ratio_list(ratios_raw),
        test_ratios=_parse_ratio_list(test_ratios_raw),
        negative_sampling=mode_raw.get('negative_sampling', default_neg_sampling.get(mode_name, 'without_replacement')),
    )


def load_ml_config(config_path):
    """Load the full pipeline config, extract and validate ml_training section.

    Args:
        config_path: Path to the species YAML config file.

    Returns:
        An MLConfig dataclass with all validated settings.
    """
    with open(config_path) as f:
        raw = yaml.safe_load(f)

    species = raw.get('species', 'unknown')
    output_base = Path(raw.get('output_dir', '.'))
    filelists_dir = output_base / 'filelists'

    ml_raw = raw.get('ml_training', {})

    # Resolve work_dir
    work_dir_raw = ml_raw.get('output_dir') or ml_raw.get('work_dir')
    if work_dir_raw:
        work_dir = Path(work_dir_raw)
    else:
        work_dir = output_base / 'For_ML' / species

    # Resolve filelist
    filelist_path = ml_raw.get('filelist_path')
    if not filelist_path:
        filelist_path = filelists_dir / 'label_filelist'

    config = MLConfig(
        species=species,
        enabled=ml_raw.get('enabled', True),
        output_dir=work_dir,
        work_dir=work_dir,
        filelist_path=Path(filelist_path) if filelist_path else None,
        test_frac=ml_raw.get('test_frac', 0.1),
        random_seed=ml_raw.get('random_seed', 42),
        blacklist_min_overlap=ml_raw.get('blacklist', {}).get('min_sample_overlap', 10),
        blacklist_max_length=ml_raw.get('blacklist', {}).get('max_length', 50),
    )

    # Build mode configs
    modes_raw = ml_raw.get('modes', {})

    germ_raw = modes_raw.get('germline', {})
    config.germline = build_mode_config(
        'germline', germ_raw,
        DEFAULT_DROP_COLUMNS['germline'],
    )

    soma_raw = modes_raw.get('somatic', {})
    config.somatic = build_mode_config(
        'somatic', soma_raw,
        DEFAULT_DROP_COLUMNS['somatic'],
    )

    return config
