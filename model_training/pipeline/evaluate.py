"""
Shared evaluation utilities for the ML training pipeline.

Provides load_data() and compute_metrics() for evaluating
AutoGluon models against labeled test data.
"""

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score

from .defaults import COLUMN_NAMES, CATEGORICAL_COLUMNS


def load_data(positive_fn, negative_fn):
    """Load test data from positive and negative files.

    Handles missing/empty files gracefully (returns empty DataFrame).
    Assigns standard column names and categorical types.
    """
    try:
        tp_data = pd.read_table(positive_fn, header=None)
        if len(tp_data) > 0:
            tp_data.columns = COLUMN_NAMES
            tp_data['class'] = 1
            np_check = True
        else:
            tp_data = pd.DataFrame()
            np_check = False
    except Exception:
        tp_data = pd.DataFrame()
        np_check = False

    try:
        fp_data = pd.read_table(negative_fn, header=None)
        if len(fp_data) > 0:
            fp_data.columns = COLUMN_NAMES
            fp_data['class'] = 0
            nn_check = True
        else:
            fp_data = pd.DataFrame()
            nn_check = False
    except Exception:
        fp_data = pd.DataFrame()
        nn_check = False

    if not np_check and not nn_check:
        return pd.DataFrame()

    if not np_check:
        ret_data = fp_data
    elif not nn_check:
        ret_data = tp_data
    else:
        ret_data = pd.concat([tp_data, fp_data], ignore_index=True)

    for col in CATEGORICAL_COLUMNS:
        if col in ret_data.columns:
            ret_data[col] = ret_data[col].astype('category')

    if 'numSeg' in ret_data.columns:
        ret_data['numSeg'] = ret_data['numSeg'].astype('category')
    if 'numSeg2' in ret_data.columns:
        ret_data['numSeg2'] = ret_data['numSeg2'].astype('category')

    return ret_data


def compute_metrics(y_true, y_pred, pos_label=1):
    """Compute precision, recall, F1, and macro-F1 scores.

    Returns dict with keys: precision, recall, f1, macro_f1.
    """
    return {
        'precision': precision_score(y_true, y_pred, pos_label=pos_label),
        'recall': recall_score(y_true, y_pred, pos_label=pos_label),
        'f1': f1_score(y_true, y_pred, pos_label=pos_label),
        'macro_f1': f1_score(y_true, y_pred, average='macro'),
    }
