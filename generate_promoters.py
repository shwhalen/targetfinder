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

expression_cutoff = 0.3
idr_cutoff = 0.1

segmentation_df = chromatics.read_bed(glob('../segmentation/*.bed.gz')[0], usecols = range(4), names = chromatics.generic_bed_columns)
promoter_states = {'TSS', '1_TssA'}
chromhmm_promoters_df = segmentation_df.query('name in @promoter_states').copy()
chromhmm_promoters_df.columns = chromatics.promoter_bed_columns
chromatics.add_names(chromhmm_promoters_df, cell_line = cell_line)

# take a single gencode tss per gene
tss_df = pd.read_csv('../../expression/gencode.v19.TSS.notlow.gff.gz', sep = '\t', header = None, usecols = [0, 3, 6, 8])
tss_df.columns = ['gene_chrom', 'gene_tss', 'strand', 'attributes']
tss_df['gene_id'] = tss_df['attributes'].str.extract('gene_id ([^ ]+)')
pos_strand_tss_df = tss_df.query('strand == "+"').groupby('gene_id')[['gene_chrom', 'gene_tss']].min()
neg_strand_tss_df = tss_df.query('strand == "-"').groupby('gene_id')[['gene_chrom', 'gene_tss']].max()
final_tss_df = pd.concat([pos_strand_tss_df, neg_strand_tss_df]) # don't ignore index

# gencode expression
expression_df = pd.read_csv('../../expression/gencodev19_genes_with_RPKM_and_npIDR_oct2014.txt.gz', sep = ' ')

# pre-filter columns for large speedup
cell_line_mask = expression_df.columns.to_series().str.contains(cell_line)
cell_line_mask['gene_id'] = True
expression_df = pd.melt(expression_df.loc[:, cell_line_mask], id_vars = 'gene_id').set_index('gene_id')

# cannot split on commas since a few localizations also use commas, don't really need lab ids so this is ok
expression_variables_df = expression_df['variable'].str.split('[:.]', expand = True)
expression_variables_df.columns = ['lab_ids', 'rna_extract', 'cell_line', 'localization']

expression_values_df = expression_df['value'].str.split('[:]', expand = True)
expression_values_df.columns = ['rpkm1', 'rpkm2', 'idr']

# extract relevant expression values
expression_df = pd.concat([expression_df.drop(['variable', 'value'], axis = 1), expression_variables_df, expression_values_df], axis = 1)

# grab polyA+ genes in the cell since cytosol doesn't have replicates for all cell lines
expression_df = expression_df.query('rna_extract == "longPolyA" and localization == "cell"')
expression_df['rpkm1'] = pd.to_numeric(expression_df['rpkm1'])
expression_df['rpkm2'] = pd.to_numeric(expression_df['rpkm2'])
expression_df['idr'] = pd.to_numeric(expression_df['idr'], errors = 'coerce')

# drop inconsistently expressed genes and genes with low expression using cutoff from Ramskold et al., "An Abundance of Ubiquitously Expressed Genes Revealed by Tissue Transcriptome Sequence Data", PLoS Comp Bio 2009
print('{:.2%} of genes exceed IDR cutoff'.format(expression_df.eval('idr > 0.1').sum() / len(expression_df)))
print('expression cutoff: {} rpkm'.format(expression_cutoff))
expression_df = expression_df.query('idr <= @idr_cutoff and ((rpkm1 + rpkm2) / 2) > @expression_cutoff')

# combine tss and expression data
active_promoters_df = pd.concat([expression_df, final_tss_df], axis = 1, join = 'inner').reset_index()
active_promoters_df['gene_tss_dup'] = active_promoters_df['gene_tss']
active_promoters_df = active_promoters_df[['gene_chrom', 'gene_tss', 'gene_tss_dup', 'gene_id']]
chromatics.write_bed(active_promoters_df, 'tss.bed')

# find active chromhmm promoters
promoters_df = chromatics.bedtools('intersect -wa -u', chromhmm_promoters_df, active_promoters_df)

# optionally expand enhancer coordinates
promoters_df['promoter_start'] -= config['promoter_extension_size'] if 'promoter_extension_size' in config else 0
promoters_df['promoter_end'] += config['promoter_extension_size'] if 'promoter_extension_size' in config else 0

# save
assert promoters_df.duplicated().sum() == 0
chromatics.write_bed(promoters_df, 'promoters.bed')
print(promoters_df.eval('promoter_end - promoter_start').describe(), '\n')
