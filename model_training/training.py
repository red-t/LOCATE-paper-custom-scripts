from autogluon.tabular import TabularDataset, TabularPredictor
import pandas as pd
import os
import argparse

def load_train_data(tp_fn, fp_fn, n=0):
    tp_data = pd.read_table(tp_fn, header=None)
    tp_data.columns = ['chrom', 'refStart', 'refEnd', 'cid', 'freq', 'strand', 'sample_id', 'startIndex', 'endIndex',
                       'numSeg', 'directionFlag', 'cltType', 'locationType', 'numSegType', 'entropy', 'balanceRatio',
                       'lowMapQualFrac', 'dualClipFrac', 'alnFrac1', 'alnFrac2', 'alnFrac4', 'alnFrac8', 'alnFrac16',
                       'meanMapQual', 'meanAlnScore', 'meanQueryMapFrac', 'meanDivergence', 'avg_de', 'bgDiv', 'back_de',
                       'bgDepth', 'bgReadLen', 'teAlignedFrac', 'flag', 'teTid']
    if n:
        tp_data = tp_data.sample(n, replace=False)
    tp_data['class'] = 1

    fp_data = pd.read_table(fp_fn, header=None)
    fp_data.columns = ['chrom', 'refStart', 'refEnd', 'cid', 'freq', 'strand', 'sample_id', 'startIndex', 'endIndex',
                       'numSeg', 'directionFlag', 'cltType', 'locationType', 'numSegType', 'entropy', 'balanceRatio',
                       'lowMapQualFrac', 'dualClipFrac', 'alnFrac1', 'alnFrac2', 'alnFrac4', 'alnFrac8', 'alnFrac16',
                       'meanMapQual', 'meanAlnScore', 'meanQueryMapFrac', 'meanDivergence', 'avg_de', 'bgDiv', 'back_de',
                       'bgDepth', 'bgReadLen', 'teAlignedFrac', 'flag', 'teTid']
    if n:
        fp_data = fp_data.sample(n, replace=False)
    fp_data['class'] = 0

    ret_data = pd.concat([tp_data, fp_data])
    ret_data['cltType'] = ret_data['cltType'].astype('category')
    ret_data['locationType'] = ret_data['locationType'].astype('category')
    ret_data['numSegType'] = ret_data['numSegType'].astype('category')
    ret_data['class'] = ret_data['class'].astype('category')
    return TabularDataset(ret_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Train a binary classification model using AutoGluon.")
    parser.add_argument('--positive', type=str, required=True, help="Path to the positive dataset file.")
    parser.add_argument('--negative', type=str, required=True, help="Path to the negative dataset file.")
    parser.add_argument('--model_name', type=str, required=True, help="Name of the model to save.")

    args = parser.parse_args()

    label = 'class'
    raw_data = load_train_data(args.positive, args.negative)
    train_data = raw_data.drop(['chrom', 'refStart', 'refEnd', 'cid', 'freq', 'strand', 'sample_id', 'startIndex', 'endIndex', 'directionFlag', 'cltType', 'avg_de', 'back_de', 'teAlignedFrac', 'flag', 'teTid'], axis=1)
    save_path = args.model_name
    predictor = TabularPredictor(label=label, path=save_path, problem_type='binary').fit(train_data, presets='best_quality')