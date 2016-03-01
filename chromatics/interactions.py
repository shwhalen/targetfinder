import chromatics
import numpy as np
import pandas as pd
import scipy.stats as stats

def get_enrichment(a, b, c):
    a_with_c = chromatics.bedtools('intersect -sorted -u -f 1.0', a, c)
    b_with_c = chromatics.bedtools('intersect -sorted -u -f 1.0', b, c)

    ct = np.zeros([2, 2])
    ct[0, 0] = len(b) - len(b_with_c)
    ct[1, 0] = len(a) - len(a_with_c)
    ct[0, 1] = len(b_with_c)
    ct[1, 1] = len(a_with_c)

    return stats.fisher_exact(ct)[1]

def get_labeled_enhancers(enhancers_fn, pairs_fn):
    enhancers_df = chromatics.read_bed(enhancers_fn, names = chromatics.enhancer_bed_columns)
    pairs_df = pd.read_csv(pairs_fn)

    interacting_enhancers_df = pairs_df.query('label == 1').drop_duplicates('enhancer_name')[chromatics.enhancer_bed_columns].sort_values(chromatics.enhancer_bed_columns[:3])
    noninteracting_enhancers_df = enhancers_df.query('enhancer_name not in @interacting_enhancers_df.enhancer_name').drop_duplicates('enhancer_name').sort_values(chromatics.enhancer_bed_columns[:3])

    return interacting_enhancers_df, noninteracting_enhancers_df

def get_labeled_promoters(promoters_fn, pairs_fn):
    promoters_df = chromatics.read_bed(promoters_fn, names = chromatics.promoter_bed_columns)
    pairs_df = pd.read_csv(pairs_fn)

    interacting_promoters_df = pairs_df.query('label == 1').drop_duplicates('promoter_name')[chromatics.promoter_bed_columns].sort_values(chromatics.promoter_bed_columns[:3])
    noninteracting_promoters_df = promoters_df.query('promoter_name not in @interacting_promoters_df.promoter_name').drop_duplicates('promoter_name').sort_values(chromatics.promoter_bed_columns[:3])

    return interacting_promoters_df, noninteracting_promoters_df

def add_names(df, columns = None, cell_line = None):
    if columns is None:
        columns = df.columns.tolist()

    df[columns[-1]] = df[columns[0]].astype(str) + ':' + df[columns[1]].astype(str) + '-' + df[columns[2]].astype(str)

    if cell_line is not None:
        df[columns[-1]] = cell_line + '|' + df[columns[-1]]

def add_windows(df, left_columns, right_columns, cell_line = None):
    left_chrom, left_start, left_end = left_columns[:3]
    right_chrom, right_start, right_end = right_columns[:3]

    df['window_chrom'] = df[left_chrom]
    df['window_start'] = df[[left_end, right_end]].min(axis = 1) + 1
    df['window_end'] = df[[left_start, right_start]].max(axis = 1) - 1
    chromatics.add_names(df, chromatics.window_bed_columns, cell_line)

def correct_fragment_order(df, left_columns, right_columns, flipped_mask):
    unaltered_columns = list(set(df.columns) - set(left_columns) - set(right_columns))
    corrected_df = df[left_columns + right_columns + unaltered_columns].copy()
    corrected_df[flipped_mask] = corrected_df.loc[flipped_mask, right_columns + left_columns + unaltered_columns].values
    return corrected_df

def get_interaction_types(interactions_df, interaction_id_column, left_fragment_columns, right_fragment_columns, elements):
    def get_fragment_elements(fragments_df, elements_df):
        return set(chromatics.bedtools('intersect -wa -u', fragments_df, elements_df).iloc[:, -1])

    interaction_types_df = interactions_df.copy()
    left_fragments_df = interaction_types_df[left_fragment_columns]
    right_fragments_df = interaction_types_df[right_fragment_columns]

    for element_label, elements_df in sorted(elements.items()):
        left_fragment_elements = get_fragment_elements(left_fragments_df, elements_df)
        interaction_types_df.eval('left_fragment_{} = {} in @left_fragment_elements'.format(element_label, left_fragments_df.columns[-1]))

    for element_label, elements_df in sorted(elements.items()):
        right_fragment_elements = get_fragment_elements(right_fragments_df, elements_df)
        interaction_types_df.eval('right_fragment_{} = {} in @right_fragment_elements'.format(element_label, right_fragments_df.columns[-1]))

    return interaction_types_df

def get_interaction_elements(interactions_df, interaction_id_column, left_fragment_columns, right_fragment_columns, left_elements_df, right_elements_df):
    left_element_name_column = left_elements_df.columns[-1]
    right_element_name_column = right_elements_df.columns[-1]
    assert left_element_name_column.endswith('_name') and right_element_name_column.endswith('_name')

    # include interaction ids for re-merging pairs
    left_fragments_df = interactions_df[left_fragment_columns + [interaction_id_column]]
    right_fragments_df = interactions_df[right_fragment_columns + [interaction_id_column]]

    left_elements_left_fragments_df = chromatics.bedtools('intersect -wa -wb', left_elements_df, left_fragments_df)
    right_elements_right_fragments_df = chromatics.bedtools('intersect -wa -wb', right_elements_df, right_fragments_df)
    left_first_pairs_df = pd.merge(left_elements_left_fragments_df, right_elements_right_fragments_df, on = interaction_id_column)

    right_elements_left_fragments_df = chromatics.bedtools('intersect -wa -wb', right_elements_df, left_fragments_df)
    left_elements_right_fragments_df = chromatics.bedtools('intersect -wa -wb', left_elements_df, right_fragments_df)
    right_first_pairs_df = pd.merge(right_elements_left_fragments_df, left_elements_right_fragments_df, on = interaction_id_column)

    interaction_elements_df = pd.concat([left_first_pairs_df, right_first_pairs_df], ignore_index = True)
    interaction_elements_df.drop_duplicates([left_element_name_column, right_element_name_column], inplace = True)
    return interaction_elements_df

def get_ordered_interaction_elements(interactions_df, interaction_id_column, left_fragment_columns, right_fragment_columns, left_elements_df, right_elements_df):
    left_element_name_column = left_elements_df.columns[-1]
    right_element_name_column = right_elements_df.columns[-1]
    assert left_element_name_column.endswith('_name') and right_element_name_column.endswith('_name')

    left_fragments_df = interactions_df[left_fragment_columns + [interaction_id_column]]
    left_elements_df = chromatics.bedtools('intersect -wa -wb', left_elements_df, left_fragments_df)

    right_fragments_df = interactions_df[right_fragment_columns + [interaction_id_column]]
    right_elements_df = chromatics.bedtools('intersect -wa -wb', right_elements_df, right_fragments_df)

    return pd.merge(left_elements_df, right_elements_df, on = interaction_id_column)

def test_correct_fragment_order():
    left_fragment_bed_columns = ['f1_' + _ for _ in chromatics.generic_bed_columns]
    right_fragment_bed_columns = ['f2_' + _ for _ in chromatics.generic_bed_columns]
    interactions_df = chromatics.read_bed('flipped.bed', names = left_fragment_bed_columns + ['interaction_id'] + right_fragment_bed_columns)

    corrected_df = correct_fragment_order(interactions_df, left_fragment_bed_columns, right_fragment_bed_columns, interactions_df.eval('f1_start > f2_start'))
    print(corrected_df)

def test_get_interaction_types():
    left_fragment_bed_columns = ['f1_' + _ for _ in chromatics.generic_bed_columns]
    right_fragment_bed_columns = ['f2_' + _ for _ in chromatics.generic_bed_columns]
    enhancers_df = chromatics.read_bed('enhancers.bed', names = chromatics.enhancer_bed_columns)
    promoters_df = chromatics.read_bed('promoters.bed', names = chromatics.promoter_bed_columns)
    interactions_df = chromatics.read_bed('interactions.bed', names = left_fragment_bed_columns + ['interaction_id'] + right_fragment_bed_columns)

    print(enhancers_df, '\n')
    print(promoters_df, '\n')
    print(interactions_df, '\n')
    typed_interactions_df = get_interaction_types(
        interactions_df,
        'interaction_id',
        left_fragment_bed_columns,
        right_fragment_bed_columns,
        {'enhancer': enhancers_df, 'promoter': promoters_df}
        )
    print(typed_interactions_df)

def test_get_interaction_elements():
    left_fragment_bed_columns = ['f1_' + _ for _ in chromatics.generic_bed_columns]
    right_fragment_bed_columns = ['f2_' + _ for _ in chromatics.generic_bed_columns]
    enhancers_df = chromatics.read_bed('enhancers.bed', names = chromatics.enhancer_bed_columns)
    promoters_df = chromatics.read_bed('promoters.bed', names = chromatics.promoter_bed_columns)
    interactions_df = chromatics.read_bed('interactions.bed', names = left_fragment_bed_columns + ['interaction_id'] + right_fragment_bed_columns)

    print(enhancers_df, '\n')
    print(promoters_df, '\n')
    print(interactions_df, '\n')
    interaction_elements_df = get_interaction_elements(
        interactions_df,
        'interaction_id',
        left_fragment_bed_columns,
        right_fragment_bed_columns,
        enhancers_df,
        promoters_df
        )
    print(interaction_elements_df)
    assert set(interaction_elements_df['interaction_id']) == {'chr1:250-350.chr1:5000-5050', 'chr1:2000-3000.chr1:210-290'}

if __name__ == '__main__':
    test_correct_fragment_order()
    test_get_interaction_types()
    test_get_interaction_elements()
