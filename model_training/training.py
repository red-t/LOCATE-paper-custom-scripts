import sys
from pathlib import Path

# 将项目根目录加入 sys.path，确保子进程能找到 model_training 模块
_project_root = Path(__file__).resolve().parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

from autogluon.tabular import TabularDataset, TabularPredictor
import pandas as pd
import os
import argparse

from model_training.pipeline.defaults import COLUMN_NAMES, CATEGORICAL_COLUMNS


def load_train_data(tp_fn, fp_fn, pos_label=1):
    """Load positive and negative training data.

    The input files are tab-delimited with 34 columns (after Stage 1 processing).
    Column names are assigned based on the standard COLUMN_NAMES schema.
    """
    tp_data = pd.read_table(tp_fn, header=None)
    tp_data.columns = COLUMN_NAMES
    tp_data['class'] = pos_label

    fp_data = pd.read_table(fp_fn, header=None)
    fp_data.columns = COLUMN_NAMES
    fp_data['class'] = 0

    ret_data = pd.concat([tp_data, fp_data])
    for col in CATEGORICAL_COLUMNS:
        if col in ret_data.columns:
            ret_data[col] = ret_data[col].astype('category')
    ret_data['class'] = ret_data['class'].astype('category')
    return TabularDataset(ret_data)


if __name__ == '__main__':
    default_drop = ('chrom,refStart,refEnd,cid,strand,sample_id,'
                    'startIndex,endIndex,numSeg2,directionFlag,cltType,'
                    'teAlignedFrac,teTid,isInBlacklist,probability')

    parser = argparse.ArgumentParser(description="Train a binary classification model using AutoGluon.")
    parser.add_argument('--positive', type=str, required=True, help="Path to the positive dataset file.")
    parser.add_argument('--negative', type=str, required=True, help="Path to the negative dataset file.")
    parser.add_argument('--model_name', type=str, required=True, help="Name of the model to save.")
    parser.add_argument('--drop-cols', type=str, default=default_drop,
                        help="Comma-separated column names to drop")
    parser.add_argument('--pos-label', type=int, default=1,
                        help="Positive class label value (default: 1)")

    args = parser.parse_args()

    label = 'class'
    drop_cols = [c.strip() for c in args.drop_cols.split(',')]

    raw_data = load_train_data(args.positive, args.negative, pos_label=args.pos_label)
    train_data = raw_data.drop(columns=drop_cols)
    save_path = args.model_name
    predictor = TabularPredictor(label=label, path=save_path, problem_type='binary').fit(train_data, presets='best_quality')