import chromatics
import json
import os
import pandas as pd

def add_enhancer_distance_to_promoter(df, bin_count = None, bins = None):
    df['window_start'] = df[['promoter_end', 'enhancer_end']].min(axis = 1) + 1
    df['window_end'] = df[['promoter_start', 'enhancer_start']].max(axis = 1) - 1
    df['enhancer_distance_to_promoter'] = 0

    non_overlapping_mask = df.eval('promoter_end < enhancer_start or enhancer_end < promoter_start')
    df.loc[non_overlapping_mask, 'enhancer_distance_to_promoter'] = df[non_overlapping_mask].eval('window_end - window_start')

    # bin_count is set for positives
    # bins is returned from positives and re-used for negatives
    # neither is used for genome-wide predictions
    if bin_count is not None:
        df['bin'], bins = pd.qcut(
            df['enhancer_distance_to_promoter'],
            bin_count,
            precision = 1,
            retbins = True)
        return bins
    elif bins is not None:
        df['bin'] = pd.cut(
            df['enhancer_distance_to_promoter'],
            bins,
            precision = 1,
            retbins = False,
            include_lowest = True)

def parse_config(config_fn):
    config = json.load(open(config_fn))
    if 'base_config_fn' in config:
        base_config_fn = os.path.join(os.path.dirname(config_fn), config['base_config_fn'])
        base_config = json.load(open(base_config_fn))
        base_config.update(config)
        config = base_config
    if 'working_dir' in config:
        config['working_dir'] = os.path.join(os.path.dirname(config_fn), config['working_dir'])
    else:
        config['working_dir'] = os.path.dirname(config_fn)
    return config

# pipeline parameters
min_enhancer_distance_to_promoter = 10000
max_enhancer_distance_to_promoter = 2000000
cell_lines = ['K562', 'GM12878', 'HeLa-S3', 'HUVEC', 'IMR90', 'NHEK', 'combined']
