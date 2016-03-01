#!/usr/bin/env python

import chromatics
import common
import os
import pandas as pd
import sys

from glob import glob

config_fn = sys.argv[1]
cell_line = config_fn.split('/')[0]
config = common.parse_config(config_fn)
os.chdir(os.path.expanduser(config['working_dir']))

peaks_fn = 'peaks.bed.gz'
methylation_fn = 'methylation.bed.gz'
cage_fn = 'cage.bed.gz'
generators = []

# preprocess peaks
if os.path.exists('../peaks'):
    assays = []
    for name, filename, source, accession in pd.read_csv('../peaks/filenames.csv').itertuples(index = False):
        columns = chromatics.narrowpeak_bed_columns if filename.endswith('narrowPeak') else chromatics.broadpeak_bed_columns
        assay_df = chromatics.read_bed('../peaks/{}.gz'.format(filename), names = columns, usecols = chromatics.generic_bed_columns + ['signal_value'])
        assay_df['name'] = name
        assays.append(assay_df)
    peaks_df = pd.concat(assays, ignore_index = True)
    chromatics.write_bed(peaks_df, peaks_fn, compression = 'gzip')
    generators.append((chromatics.generate_average_signal_features, peaks_fn))

# preprocess methylation
if os.path.exists('../methylation'):
    assays = [chromatics.read_bed(_, names = chromatics.methylation_bed_columns, usecols = chromatics.generic_bed_columns + ['mapped_reads', 'percent_methylated']) for _ in glob('../methylation/*.bed.gz')]
    methylation_df = pd.concat(assays, ignore_index = True).query('mapped_reads >= 10 and percent_methylated > 0')
    methylation_df['name'] = 'Methylation'
    del methylation_df['mapped_reads']
    chromatics.write_bed(methylation_df, methylation_fn, compression = 'gzip')
    generators.append((chromatics.generate_average_signal_features, methylation_fn))

# preprocess cage
if os.path.exists('../cage'):
    cage_df = chromatics.read_bed(glob('../cage/*.bed.gz')[0], names = chromatics.cage_bed_columns, usecols = chromatics.cage_bed_columns[:5])
    cage_df['name'] = 'CAGE'
    chromatics.write_bed(cage_df, cage_fn, compression = 'gzip')
    generators.append((chromatics.generate_average_signal_features, cage_fn))

# generate features
pairs_df = pd.read_csv('pairs.csv')
assert pairs_df.duplicated().sum() == 0
training_df = chromatics.generate_training(pairs_df, config['regions'], generators, chunk_size = 2**14, n_jobs = 1)

# save
training_df.to_hdf('training.h5', 'training', mode = 'w', complevel = 1, complib = 'zlib')
