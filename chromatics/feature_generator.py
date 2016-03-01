#!/usr/bin/env python

import chromatics
import numpy as np
import pandas as pd
import sklearn.externals.joblib as joblib

def generate_average_signal_features(chunk_df, region, dataset):
    assert (chunk_df[region + '_end'] > chunk_df[region + '_start']).all()

    region_bed_columns = ['{}_{}'.format(region, _) for _ in chromatics.generic_bed_columns]
    signal_df = chromatics.bedtools('intersect -wa -wb', chunk_df[region_bed_columns].drop_duplicates(region + '_name'), dataset, right_names = chromatics.signal_bed_columns)

    group_columns = ['{}_{}'.format(region, _) for _ in ['name', 'start', 'end']] + ['dataset']
    average_signal_df = signal_df.groupby(group_columns, sort = False, as_index = False).aggregate({'signal_value': sum})
    average_signal_df['signal_value'] /= average_signal_df[region + '_end'] - average_signal_df[region + '_start']
    average_signal_df['dataset'] += ' ({})'.format(region)

    return average_signal_df.pivot_table(index = region + '_name', columns = 'dataset', values = 'signal_value')

def generate_chunk_features(pairs_df, regions, generators, chunk_size, chunk_number, max_chunks):
    print(chunk_number, max_chunks - 1)

    chunk_lower_bound = chunk_number * chunk_size
    chunk_upper_bound = chunk_lower_bound + chunk_size
    chunk_df = pairs_df.iloc[chunk_lower_bound:chunk_upper_bound]
    assert 0 < len(chunk_df) <= chunk_size

    index_columns = ['{}_name'.format(region) for region in regions]
    features_df = chunk_df[index_columns]
    for region in regions:
        region_features = [generator(chunk_df, region, dataset) for generator, dataset in generators]
        region_features_df = pd.concat(region_features, axis = 1)
        features_df = pd.merge(features_df, region_features_df, left_on = '{}_name'.format(region), right_index = True, how = 'left')
    return features_df.set_index(index_columns)

def generate_training(pairs_df, regions, generators, chunk_size = 2**16, n_jobs = -1):
    for region in regions:
        region_bed_columns = {'{}_{}'.format(region, _) for _ in chromatics.generic_bed_columns}
        assert region_bed_columns.issubset(pairs_df.columns)

    max_chunks = int(np.ceil(len(pairs_df) / chunk_size))
    results = joblib.Parallel(n_jobs)(
        joblib.delayed(generate_chunk_features)(pairs_df, regions, generators, chunk_size, chunk_number, max_chunks)
        for chunk_number in range(max_chunks))

    features_df = pd.concat(results).fillna(0)
    training_df = pd.merge(pairs_df, features_df, left_on = ['{}_name'.format(region) for region in regions], right_index = True)
    assert training_df.index.is_unique
    assert training_df.columns.is_unique
    return training_df

def get_random_pairs(pair_count, region_a_prefix = 'r1', region_b_prefix = 'r2', random_state = 0):
    random_state = np.random.RandomState(random_state)
    f1_start = random_state.randint(0, 1e6, pair_count)
    f2_start = random_state.randint(0, 1e6, pair_count)
    pair_coordinates = [
        random_state.choice(chromatics.chroms, pair_count),
        f1_start,
        f1_start + random_state.randint(1, 50000, pair_count),
        ['{}_{}'.format(region_a_prefix, _) for _ in range(pair_count)],
        random_state.choice(chromatics.chroms, pair_count),
        f2_start,
        f2_start + random_state.randint(1, 50000, pair_count),
        ['{}_{}'.format(region_b_prefix, _) for _ in range(pair_count)]
        ]
    pair_columns = ['{}_{}'.format(region_a_prefix, _) for _ in chromatics.generic_bed_columns] + \
        ['{}_{}'.format(region_b_prefix, _) for _ in chromatics.generic_bed_columns]
    return pd.DataFrame(dict(zip(pair_columns, pair_coordinates)), columns = pair_columns)

def test_generate_average_signal_features():
    enhancers_df = chromatics.read_bed('enhancers.bed', names = chromatics.enhancer_bed_columns)
    average_signal_df = generate_average_signal_features(enhancers_df, 'enhancer', 'peaks.bed')
    assert average_signal_df.loc['enhancer1', 'RAD21 (enhancer)'] == 2.1
    assert average_signal_df.loc['enhancer1', 'CTCF (enhancer)'] == 0.502
    print(average_signal_df)

def test_generate_training():
    regions = ['enhancer', 'promoter']
    pairs_df = get_random_pairs(100, regions[0], regions[1])
    print(pairs_df.head())

    signal_df = chromatics.read_bed('wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak.gz', names = chromatics.narrowpeak_bed_columns, usecols = ['chrom', 'start', 'end', 'signal_value'])
    signal_df['dataset'] = 'DNase'
    signal_df = signal_df[chromatics.signal_bed_columns]
    generators = [(generate_average_signal_features, signal_df)]
    training_df = generate_training(pairs_df, regions, generators, chunk_size = len(pairs_df) // 2)
    print(training_df.head(), '\n')

if __name__ == '__main__':
    test_generate_average_signal_features()
    test_generate_training()
