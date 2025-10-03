# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:32:59 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp
import matplotlib as mpl
import numpy as np
from scipy import stats

# ========== USER SETTINGS FOR STYLING ==========
width_mm = 40  # figure width in mm
height_mm = 60  # figure height in mm
font_size = 8   # desired font size
font_family = "Arial"

# Convert mm to inches
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Set global font properties
mpl.rcParams['font.family'] = font_family
plt.rcParams.update({'font.size': font_size})
# ===============================================

# 1) Read the combined data
df_all = pd.read_excel("Figure_3_b.xlsx", sheet_name="Figure_3_b", engine="openpyxl")

# 2) Build a dictionary to hold each condition’s DataFrame
growth_rate_dfs = {}
doubling_times = []
conditions = []

unique_conditions = df_all["Condition"].unique()
for cond in unique_conditions:
    # Subset for each condition
    sub_df = df_all[df_all["Condition"] == cond].copy()
    growth_rate_dfs[cond] = sub_df
    
    # For Kruskal–Wallis and Dunn’s tests
    doubling_times.append(sub_df["doubling_time_min"].to_numpy())
    conditions.extend([cond] * len(sub_df))

# 3) Kruskal–Wallis test
stat, p = kruskal(*doubling_times)
print("Kruskal-Wallis test:")
print(f"  Statistic: {stat}")
print(f"  P-value:   {p}")

# 4) Dunn’s post-hoc test (Bonferroni correction)
dunn_result = sp.posthoc_dunn(doubling_times, p_adjust='bonferroni')
print("\nDunn's post-hoc test results:")
print(dunn_result)

# 5) Print which conditions were loaded
print("\nLoaded growth rate details for conditions:")
for c in growth_rate_dfs.keys():
    print(" ", c)

# 6) Create the same scatter+mean figure
fig = plt.figure(figsize=(width_in, height_in), dpi=1200)

# Colors for each condition
colors = plt.cm.tab20.colors

# Plot data
for idx, (condition, df_cond) in enumerate(growth_rate_dfs.items()):
    x_vals = np.random.normal(idx, 0.1, size=len(df_cond['doubling_time_min']))
    y_vals = df_cond['doubling_time_min']
    
    # Scatter plot of individual points
    plt.scatter(
        x_vals, y_vals,
        facecolors=colors[idx % len(colors)],  # solid fill color
        edgecolors='black',                    # black outline
        linewidths=1,                          # thin edge line
        alpha=0.5,                             # partial transparency
        s=30
    )
    
    # Mean line for each condition
    if len(y_vals) > 0:
        mean_y = y_vals.mean()
        plt.plot([idx - 0.2, idx + 0.2], [mean_y, mean_y],
                 color='black', linewidth=2)

# 7) Adjust x-axis
plt.xticks(range(len(growth_rate_dfs)), growth_rate_dfs.keys(), rotation=45, ha='right')

# 8) Y-axis label and range
plt.ylabel('Fungal doubling time (min)')
plt.ylim(70, 275)

plt.tight_layout()
#plt.show()
fig.savefig("Figure_3_b.png", bbox_inches='tight', dpi=300)
