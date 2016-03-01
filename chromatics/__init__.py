from .bedtools import *
from .feature_generator import *
from .interactions import *
from .samtools import *

chroms = ['chr{}'.format(_) for _ in list(range(1, 22 + 1)) + ['X', 'Y']]

# http://genome.ucsc.edu/FAQ/FAQformat.html#format1
generic_bed_columns = ['chrom', 'start', 'end', 'name']
bed6_columns = generic_bed_columns + ['score', 'strand']
bed9_columns = bed6_columns + ['thick_start', 'thick_end', 'item_rgb']

# http://genome.ucsc.edu/FAQ/FAQformat.html#format13
broadpeak_bed_columns = bed6_columns + ['signal_value', 'p_value', 'q_value']

# http://genome.ucsc.edu/FAQ/FAQformat.html#format12
narrowpeak_bed_columns = broadpeak_bed_columns + ['peak']

# http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeHaibMethylRrbs
methylation_bed_columns = bed9_columns + ['mapped_reads', 'percent_methylated']

# http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeRikenCage
cage_bed_columns = bed6_columns + ['rpkm1', 'rpkm2', 'idr']

# custom
enhancer_bed_columns = ['enhancer_' + _ for _ in generic_bed_columns]
promoter_bed_columns = ['promoter_' + _ for _ in generic_bed_columns]
window_bed_columns = ['window_' + _ for _ in generic_bed_columns]
signal_bed_columns = generic_bed_columns[:3] + ['dataset', 'signal_value']

# GOTHiC ouptut (non-bed)
gothic_p_columns = ['chrom', 'start', 'end', 'symbol', 'ensembl_id', 'expression_quartile']
gothic_po_columns = ['f1_' + _ for _ in gothic_p_columns] + ['f2_' + _ for _ in generic_bed_columns[:-1]] + ['count', 'log_ratio']
gothic_pp_columns = ['f1_' + _ for _ in gothic_p_columns] + ['f2_' + _ for _ in gothic_p_columns] + ['count', 'log_ratio']
