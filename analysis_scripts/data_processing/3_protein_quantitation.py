""" 
Calculate abundance ratio for individual proteins according to the control sample, then filter to proteins quantified in at least 3/4 replicates. Determine significantly regulated proteins using one-sample t-test and apply standard cutoffs (fold change in abundance > 2, p-value < 0.05).
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp, ttest_ind
from collections import Counter

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_folder = f'python_results/normalisation/'
output_folder = f'python_results/protein_quantitation/'
sample_name = 'AGG'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def one_sample_ttest(compiled, sample_cols, popmean=1, index='Sequence'):
    df = compiled.copy()
    ttest_results = []
    for index_name, df in df.groupby(index):
        results = []
        for col in sample_cols:
            test_vals = df[col].values
            if len(test_vals) > 1:
                results.append(
                    list(ttest_1samp(test_vals, popmean=popmean, nan_policy='omit')) + [np.mean(test_vals)])
            else:
                results.append(tuple([np.nan, np.nan, np.nan]))
        results = pd.DataFrame(results)
        results.columns = ['t-stat', 'p-val', 'mean_val']
        results['sample'] = sample_cols
        results[index] = index_name
        ttest_results.append(results)
    ttest_results = pd.concat(ttest_results)
    ttest_results[['t-stat', 'p-val']
                  ] = ttest_results[['t-stat', 'p-val']].astype(float)

    return ttest_results  # 5% of the points detected as significantly different


# Generate treated/control ratio following method described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3751587/
# After normalization, the relative quantification for treated/control is derived from the ratio of the individual treated reporter channels over the average of the control reporter channels.

# read in normalised protein data, drop non-sample labels
proteins = pd.read_excel(f'{input_folder}{sample_name}_IRS_normalised.xlsx', sheet_name='Proteins')
proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
proteins.dropna(inplace=True, subset=['sample_name'])

# create ratiometric measure with 8A mutant
## Calculate average abundance for 8ala proteins as control
control_ratios = proteins[proteins['sample_name'].str.contains('WT_')].copy()
control_ratios = control_ratios.groupby('Proteins').mean().reset_index()
control_ratios = dict(control_ratios[['Proteins', 'med_norm']].values)

## apply control abundance to all proteins 
proteins['control_abundance'] = proteins['Proteins'].map(control_ratios)

## Calculate abundance ratio
proteins['control_ratio'] = proteins['med_norm'] / \
    proteins['control_abundance']

# Filter for proteins identified in > 3 replicates of each sample group
# group column by protein and sample (mutant type). make sure creating list for proteins for individual sample types 
#quantified_proteins = proteins.dropna(
#    subset=['control_ratio']).groupby('Proteins').count().reset_index()
#quantified_proteins = quantified_proteins[quantified_proteins['ref_norm'] >= 3]['Proteins'].tolist(
#)


sample_prot = proteins[['sample_type', 'Proteins', 'control_ratio']].copy()
quantified_proteins = sample_prot.groupby(['sample_type', 'Proteins']).count()
quantified_proteins = quantified_proteins[quantified_proteins['control_ratio'] >= 3].reset_index()
quantified_proteins = quantified_proteins.drop_duplicates('Proteins')
quantified_proteins = quantified_proteins['Proteins'].tolist()


# Quantified = [(Protein_name, 8ala), (Protein_name, WT)]

# proteins['quant'] = [1 if (protein, sample) in quantified_proteins else np.nan for protein, sample in proteins[['Proteins', 'sample_type']].values]

proteins = proteins[proteins['Proteins'].isin(quantified_proteins)]
proteins[['sample_type', 'replicate']] = proteins['sample_name'].str.split('_', expand=True)

# one sample ttest
#  Determine significantly-different proteins via one-sample t-test against pop-mean of 1.0
# Apply one-sample ttest
compiled = pd.pivot_table(proteins.copy(), index=['Proteins', 'replicate'], columns=['sample_type'], values=['control_ratio'])
compiled.columns = compiled.columns.get_level_values(1)
ttest_results = one_sample_ttest(compiled.reset_index(), sample_cols=['8ala', 'L14A', 'Ex4', '3S', '9S', 'FY', '4Y'], popmean=1, index='Proteins')
ttest_results['-log10(pval)'] = - np.log10(ttest_results['p-val'])
ttest_results['log2(abundance)'] = np.log2(ttest_results['mean_val'])

significant_proteins = ttest_results.copy()
significant_proteins = significant_proteins[significant_proteins['-log10(pval)'] >= 1.3]
significant_proteins['|log2(abundance)|'] = abs(significant_proteins['log2(abundance)'])
significant_proteins = significant_proteins[significant_proteins['|log2(abundance)|'] >= 1]

for sample, df in ttest_results.groupby('sample'):
    sns.scatterplot(data=df, x='log2(abundance)', y='-log10(pval)')
    plt.title(sample)
    plt.xlim(-5,5)
    plt.ylim(0,6)
    plt.axhline(1.3, color='gray', ls='--')
    plt.axvline(-1, color='gray', ls='--')
    plt.axvline(1, color='gray', ls='--')
    sns.scatterplot(data=significant_proteins[significant_proteins['sample'] == sample], x='log2(abundance)', y='-log10(pval)', color='red')
    #x_val = list(significant_proteins[significant_proteins['sample'] == sample]['log2(abundance)'])
    #y_val = list(significant_proteins[significant_proteins['sample'] == sample]['-log10(pval)'])
    #text = significant_proteins[significant_proteins['sample'] == sample]['Proteins']
    #for i, txt in enumerate(text):
    #    plt.annotate(txt, (x_val[i]-1, y_val[i]))
    plt.savefig(f'{output_folder}{sample}_volcano_plot.png')
    plt.show()

# Save to excel
FileHandling.df_to_excel(
    output_path=f'{output_folder}{sample_name}_quantitation.xlsx',
    sheetnames=['Proteins', 'sign_results'],
    data_frames=[proteins, ttest_results]
)

