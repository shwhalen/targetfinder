#!/usr/bin/env python

import common
import chromatics
import os
import pandas as pd
import sys

from glob import glob

config_fn = sys.argv[1]
cell_line = config_fn.split('/')[0]
config = common.parse_config(config_fn)
os.chdir(os.path.expanduser(config['working_dir']))

distance_bin_count = 5
negatives_per_bin = 20
left_fragment_columns = ['f1_' + _ for _ in chromatics.generic_bed_columns]
right_fragment_columns = ['f2_' + _ for _ in chromatics.generic_bed_columns]

# parse hi-c interactions
interactions_fn = glob('../hi-c/*looplist.txt.gz')[0]
hic_interactions_df = pd.read_csv(interactions_fn, sep = '\t', usecols = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'fdr_h'])
hic_interactions_df.columns = left_fragment_columns[:-1] + right_fragment_columns[:-1] + ['qvalue']

# fix chromosome names
hic_interactions_df['f1_chrom'] = 'chr' + hic_interactions_df['f1_chrom']
hic_interactions_df['f2_chrom'] = 'chr' + hic_interactions_df['f2_chrom']

# assign names for quick fragment searching
chromatics.add_names(hic_interactions_df, left_fragment_columns, cell_line)
chromatics.add_names(hic_interactions_df, right_fragment_columns, cell_line)
hic_interactions_df['interaction_id'] = hic_interactions_df['f1_name'] + '.' + hic_interactions_df['f2_name']

# find interactions with at least 1 promoter fragment and 1 enhancer fragment
enhancers_df = chromatics.read_bed('enhancers.bed', names = chromatics.enhancer_bed_columns)
promoters_df = chromatics.read_bed('promoters.bed', names = chromatics.promoter_bed_columns)
positives_df = chromatics.get_interaction_elements(
    hic_interactions_df,
    'interaction_id',
    left_fragment_columns,
    right_fragment_columns,
    enhancers_df,
    promoters_df)

# re-arrange columns so enhancer is always on the left
positives_df = positives_df[chromatics.enhancer_bed_columns + chromatics.promoter_bed_columns + ['interaction_id']]

# cap interaction distance
common.add_enhancer_distance_to_promoter(positives_df)
positives_df = positives_df.query('@common.min_enhancer_distance_to_promoter < enhancer_distance_to_promoter < @common.max_enhancer_distance_to_promoter')

# add distance bins
positive_bins = common.add_enhancer_distance_to_promoter(positives_df, bin_count = distance_bin_count)
positives_df['label'] = 1

# generate negative candidates from all pairs of enhancers and active promoters
negative_candidates_df = pd.merge(enhancers_df, promoters_df, left_on = 'enhancer_chrom', right_on = 'promoter_chrom')
negative_candidates_df = negative_candidates_df.query('not (enhancer_name in @positives_df.enhancer_name and promoter_name in @positives_df.promoter_name)')
print('enhancers: {} active promoters: {} negative candidate pairs: {}'.format(enhancers_df.shape[0], promoters_df.shape[0], negative_candidates_df.shape[0]))

# cap interaction distance: distances outside positive bins will be assigned NA, so drop NAs
common.add_enhancer_distance_to_promoter(negative_candidates_df, bins = positive_bins)
negative_candidates_df.dropna(inplace = True)

print('\npositive distance bins:')
print(positives_df['bin'].value_counts())

print('\nnegative candidate distance bins:')
print(negative_candidates_df['bin'].value_counts())

# distance match negatives to positives
fewest_binned_positives = positives_df['bin'].value_counts().min()
negatives_df = negative_candidates_df.groupby('bin', as_index = False).apply(lambda x: x.sample(fewest_binned_positives * negatives_per_bin, random_state = 0))
negatives_df['label'] = 0

# combine negatives with positives and remove potential overlap
pairs_df = pd.concat([positives_df, negatives_df], ignore_index = True)
pairs_df.drop_duplicates(['enhancer_name', 'promoter_name'], keep = 'first', inplace = True)

# print stats and save
print('\nenhancer lengths:')
print(pairs_df.eval('enhancer_end - enhancer_start').describe())
assert pairs_df.eval('enhancer_end - enhancer_start').max() <= enhancers_df.eval('enhancer_end - enhancer_start').max()

print('\npromoter lengths:')
print(pairs_df.eval('promoter_end - promoter_start').describe())
assert pairs_df.eval('promoter_end - promoter_start').max() <= promoters_df.eval('promoter_end - promoter_start').max()

print('\nwindow lengths:')
print(pairs_df['enhancer_distance_to_promoter'].describe())

print('\ndistance bins:')
print(pairs_df.groupby('label')['bin'].value_counts())

print('\nclasses:')
print(pairs_df['label'].value_counts())

print('\nenhancers per promoter (positives only):')
print(pairs_df.query('label == 1').groupby('promoter_name')['enhancer_name'].nunique().describe())

print('\npromoters per enhancer (positives only):')
print(pairs_df.query('label == 1').groupby('enhancer_name')['promoter_name'].nunique().describe(), '\n')

pairs_df['window_chrom'] = pairs_df['enhancer_chrom']
chromatics.add_names(pairs_df, chromatics.window_bed_columns, cell_line)

# add a few useful features here -- positive interactions already in the window
interactions_in_window_df = chromatics.bedtools('coverage -counts -F 1.0'.format(cell_line), pairs_df[chromatics.window_bed_columns], pairs_df.query('label == 1')[chromatics.window_bed_columns], left_names = chromatics.window_bed_columns, right_names = ['interactions_in_window']).iloc[:, -2:]
assert len(pairs_df) == len(interactions_in_window_df)
pairs_df = pd.merge(pairs_df, interactions_in_window_df, on = 'window_name')

# active genes skipped over by the loop
active_promoters_in_window_df = chromatics.bedtools('coverage -counts'.format(cell_line), pairs_df[chromatics.window_bed_columns], promoters_df, left_names = chromatics.window_bed_columns, right_names = ['active_promoters_in_window']).iloc[:, -2:]
assert len(pairs_df) == len(active_promoters_in_window_df)
pairs_df = pd.merge(pairs_df, active_promoters_in_window_df, on = 'window_name').drop('interaction_id', axis = 1)

# save
assert pairs_df.duplicated().sum() == 0
pairs_df.to_csv('pairs.csv', index = False)
