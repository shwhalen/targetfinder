#!/usr/bin/env python

import gzip
import pandas as pd
import sys

cell_line = sys.argv[1]
region = sys.argv[2]

training_df = pd.read_hdf('{}/output-{}/training.h5'.format(cell_line, region), 'training')
training_df.to_csv(gzip.open('{}/output-{}/training.csv.gz'.format(cell_line, region), 'wt'), index = False)
