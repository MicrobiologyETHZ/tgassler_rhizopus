# -*- coding: utf-8 -*-
"""
Created on Wed May 21 18:15:15 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams

# --- Set figure size in millimeters ---
fig_width_mm = 100  # Increased width to give space for long labels
fig_height_mm = 80
figsize_inches = (fig_width_mm / 25.4, fig_height_mm / 25.4)

# --- Set font and fontsize globally ---
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 6

# --- Create figure ---
fig, ax = plt.subplots(figsize=figsize_inches, dpi=1200)

# --- Load and prepare data ---
def prepare_data(file_path, N=5):
    df = pd.read_csv(file_path, sep='\t')
    df_sorted = df.sort_values(by='strength', ascending=False)
    return df_sorted.head(N)

file_paths = ["Both_Up.tsv", "Both_Down.tsv"]
dfs = []

for file_path in file_paths:
    df = prepare_data(file_path)
    sample_name = file_path.split('_')[1].split('.')[0]
    df['sample'] = sample_name
    dfs.append(df)

combined_df = pd.concat(dfs)
samples = combined_df['sample'].unique()[:2]
go_terms = combined_df['term description'].unique()

# --- Color mapping ---
fdr_values = combined_df['false discovery rate']
norm = Normalize(vmin=fdr_values.min(), vmax=fdr_values.max())
cmap = plt.get_cmap('viridis')

# --- Bubble plot ---
for i, sample in enumerate(samples):
    sample_df = combined_df[combined_df['sample'] == sample]
    x_offset = 0.3 if i == 0 else -0.55
    x = [i + x_offset] * len(sample_df)
    y = [go_terms.tolist().index(go) for go in sample_df['term description']]
    bubble_sizes = sample_df['strength'] * 75
    bubble_colors = cmap(norm(sample_df['false discovery rate']))
    ax.scatter(x, y, s=bubble_sizes, c=bubble_colors, label=sample, alpha=0.8)

# --- Axis customization ---
ax.set_xticks([0.25, 0.5])
ax.set_xticklabels(['Up', 'Down'])

# Move y-axis to the right and align labels properly
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_yticks(range(len(go_terms)))
ax.set_yticklabels(go_terms, ha='left')  # Align label text left
ax.tick_params(axis='y', labelrotation=0, pad=5)  # No rotation; add padding

# --- Color bar ---
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.1)
cbar.set_label('FDR')

# --- Layout adjustments ---
plt.subplots_adjust(left=0.1, right=0.92, top=0.95, bottom=0.2)
plt.tight_layout()
#plt.show()
plt.savefig("Supplementary_Figure_10.png", bbox_inches='tight', dpi=300)