"""
Standard column schema and mode-specific default parameters.

Current simulation pipeline output has 33 columns (the "v4" format).
After Stage 1 (preprocess) inserts sample_id at position 7 (0-based index 6),
the training files have 34 columns. This schema defines the standard names
for the 34-column format (post-insertion).

NOTE: Column 9 in raw (col 10 after insertion) is a duplicate of column 5's
numSeg. It is named 'numSeg2' to avoid pandas duplicate column name issues.
"""

# Standard 34-column names (after sample_id insertion at position 7)
COLUMN_NAMES = [
    'chrom',           # 1
    'refStart',        # 2
    'refEnd',          # 3
    'cid',             # 4  (cltID / cluster ID)
    'numSeg',          # 5  (number of segments, normalized)
    'strand',          # 6  (cluster orientation)
    'sample_id',       # 7  (inserted by preprocess)
    'startIndex',      # 8  (0-based, include)
    'endIndex',        # 9  (0-based, not-include)
    'numSeg2',         # 10 (duplicate of col 5)
    'directionFlag',   # 11 (1=forward, 2=reverse)
    'cltType',         # 12 (0=germline, >0=somatic)
    'locationType',    # 13 (bitwise location flag)
    'numSegType',      # 14 (number of segment types)
    'entropy',         # 15
    'balanceRatio',    # 16
    'lowMapQualFrac',  # 17
    'dualClipFrac',    # 18
    'alnFrac1',        # 19
    'alnFrac2',        # 20
    'alnFrac4',        # 21
    'alnFrac8',        # 22
    'alnFrac16',       # 23
    'meanMapQual',     # 24
    'meanAlnScore',    # 25
    'meanQueryMapFrac',# 26
    'meanDivergence',  # 27
    'bgDiv',           # 28
    'bgDepth',         # 29
    'bgReadLen',       # 30
    'teAlignedFrac',   # 31
    'teTid',           # 32
    'isInBlacklist',   # 33
    'probability',     # 34
]

# Default drop columns per mode
# These are columns excluded from AutoGluon training features.
# Germline: keeps numSeg as a feature, drops directionFlag
# Somatic: drops numSeg (both copies), directionFlag, and bgDepth
DEFAULT_DROP_COLUMNS = {
    'germline': [
        'chrom', 'refStart', 'refEnd', 'cid',
        'strand', 'sample_id', 'startIndex', 'endIndex', 'numSeg2',
        'directionFlag', 'cltType',
        'teAlignedFrac', 'teTid', 'isInBlacklist', 'probability',
    ],
    'somatic': [
        'chrom', 'refStart', 'refEnd', 'cid',
        'strand', 'sample_id', 'startIndex', 'endIndex',
        'numSeg', 'numSeg2', 'directionFlag',
        'bgDepth', 'teAlignedFrac', 'teTid', 'isInBlacklist', 'probability',
    ],
}

# Column names that should be categorical type
CATEGORICAL_COLUMNS = ['cltType', 'locationType', 'numSegType']
