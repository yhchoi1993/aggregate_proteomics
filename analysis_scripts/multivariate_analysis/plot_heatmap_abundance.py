## script for plotting the log abundance of proteins for proteomics
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
from GEN_Utils import FileHandling
from seaborn.palettes import color_palette
from seaborn.widgets import choose_colorbrewer_palette
from lmfit import Model
from scipy.optimize import curve_fit

logger.info('Import OK')

input_path = f'exp_data/heatmap_proteins_100_2.xlsx'
output_folder = f'python_results/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# read in csv file for plotting
abundance = pd.read_excel(input_path, sheet_name='All')
abundance = abundance.drop(['total_abundance', 'logAb'], axis=1)
abundance = abundance.set_index('Proteins')
logabundance = np.log10(abundance)

# plot scatter plot for phase diagram
fig, ax = plt.subplots(figsize=(20,20))
sns.heatmap(data=logabundance, cmap='coolwarm')
plt.savefig(f'{output_folder}log_ab_heatmap.pdf', transparent=True)
plt.show()
