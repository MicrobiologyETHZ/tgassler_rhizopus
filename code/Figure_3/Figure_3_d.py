# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:34:21 2025

@author: tgassler
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import matplotlib as mpl
import matplotlib.ticker as ticker

# --- USER SETTINGS: size in mm, font size ---
width_mm = 80   # e.g., 80 mm wide
height_mm = 45  # e.g., 45 mm tall
font_size = 8
font_family = "Arial"

# Convert mm -> inches
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Set font properties
mpl.rcParams['font.family'] = font_family
plt.rcParams.update({'font.size': font_size})

# ========== 1) READ SOURCE DATA FROM EXCEL ==========
df_all_spores = pd.read_excel("Figure_3_d.xlsx", sheet_name="Figure_3_d", engine="openpyxl")

# ========== 2) MANN–WHITNEY TESTS (Yes vs. No FOR EACH CONDITION) ==========
p_values = {}
unique_conditions = df_all_spores['Condition'].unique()

for condition in unique_conditions:
    df_cond = df_all_spores[df_all_spores['Condition'] == condition]
    growing_vals = df_cond[df_cond['Growth Status'] == 'Yes']['Initial Intensity']
    nongrowing_vals = df_cond[df_cond['Growth Status'] == 'No']['Initial Intensity']
    
    # Ensure both sets are non-empty before testing
    if len(growing_vals) > 0 and len(nongrowing_vals) > 0:
        # SciPy's default is two-sided in recent versions; this matches your original call
        _, p_val = mannwhitneyu(growing_vals, nongrowing_vals)
        p_values[condition] = p_val
    else:
        p_values[condition] = None

# ---- Benjamini–Hochberg FDR adjustment across the within-condition tests ----
def benjamini_hochberg(pdict):
    items = [(cond, p) for cond, p in pdict.items() if p is not None]
    m = len(items)
    if m == 0:
        return {cond: None for cond in pdict.keys()}
    # sort by raw p ascending
    items_sorted = sorted(items, key=lambda x: x[1])
    q_vals = {}
    prev = 1.0
    # work from largest rank to smallest to enforce monotonicity
    for i in range(m - 1, -1, -1):
        cond, p = items_sorted[i]
        rank = i + 1
        q = p * m / rank
        prev = min(prev, q)
        q_vals[cond] = min(1.0, prev)
    # include Nones for conditions without a p-value
    for cond in pdict.keys():
        if cond not in q_vals:
            q_vals[cond] = None
    return q_vals

q_bh = benjamini_hochberg(p_values)

print("Mann–Whitney U test results (Condition: raw p-value):")
for cond in unique_conditions:
    print(f"  {cond}: p = {p_values[cond]}")

print("\nBenjamini–Hochberg FDR-adjusted p-values (q-values):")
for cond in unique_conditions:
    q = q_bh.get(cond, None)
    print(f"  {cond}: q = {q:.4g}" if q is not None else f"  {cond}: q = NA")

# ========== 3) CREATE SUBPLOTS (ONE PER CONDITION) ==========
fig, axes = plt.subplots(
    1,
    len(unique_conditions),
    figsize=(width_in, height_in),
    dpi=1200,
    sharey=True
)

# If there's only 1 condition, wrap the single Axes in a list for consistent handling
if len(unique_conditions) == 1:
    axes = [axes]

# ========== 4) PLOT BOX + STRIP FOR EACH CONDITION ==========
for i, condition in enumerate(unique_conditions):
    df_cond = df_all_spores[df_all_spores['Condition'] == condition]
    
    # Boxplot
    sns.boxplot(
        x='Growth Status',
        y='Initial Intensity',
        data=df_cond,
        ax=axes[i],
        palette="Set2",
        hue='Growth Status',
        legend=False,
        showfliers=False
    )
    
    # Stripplot (individual data points)
    sns.stripplot(
        x='Growth Status',
        y='Initial Intensity',
        data=df_cond,
        ax=axes[i],
        color='black',
        alpha=0.5,
        size=3,
        jitter=True,
        dodge=True
    )
    
    # Centered subplot title with some extra space above
    axes[i].set_title(condition, y=1.12, ha='center')
    
    # Remove x-axis label (already labeled 'Growth Status')
    axes[i].set_xlabel('')

# Shared y-axis label on the first subplot
axes[0].set_ylabel('Bacterial Load (a.u.)')

# ========== 5) SCIENTIFIC NOTATION & LAYOUT ADJUSTMENTS ==========
# Force scientific notation on Y-axis for all subplots
for ax in axes:
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# Switch to MathText so '1e5' appears as '1×10^5'
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))  # always use scientific notation
axes[-1].yaxis.set_major_formatter(formatter)  # e.g., apply to the last subplot or any you prefer

# Increase top margin so titles aren’t clipped
plt.tight_layout(rect=[0, 0, 1, 1])
#plt.show()
fig.savefig("Figure_3_d.png", bbox_inches='tight', dpi=300)

print("\nNumber of 'growing' and 'not growing' spores per condition:")
for condition in unique_conditions:
    df_cond = df_all_spores[df_all_spores['Condition'] == condition]
    n_growing = (df_cond['Growth Status'] == 'Yes').sum()
    n_not_growing = (df_cond['Growth Status'] == 'No').sum()
    print(f"  {condition}: grew = {n_growing}, not grew = {n_not_growing}")

