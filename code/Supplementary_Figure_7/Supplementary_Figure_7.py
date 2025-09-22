# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 11:55:17 2025

@author: tgassler
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.formula.api import ols
import statsmodels.api as sm

# ==== USER SETTINGS ====
fig_width_mm = 60
fig_height_mm = 120
font_size = 10
font_family = "Arial"
# =======================

# Convert mm to inches
fig_size_inches = (fig_width_mm / 25.4, fig_height_mm / 25.4)
plt.rcParams.update({"font.size": font_size, "font.family": font_family})

# Load Excel data
file_path = r"sourcedata_SF7g.xlsx"
df = pd.read_excel(file_path)

# Log-transform the load values
df["log_load"] = np.log10(df["load"])

# ANOVA
model = ols('log_load ~ C(Sample)', data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

# Tukey test
mc = MultiComparison(df['log_load'], df['Sample'])
tukey = mc.tukeyhsd(alpha=0.05)

# Print results in the terminal
print("ANOVA results:\n", anova_table, "\n")
print("Tukey HSD results:\n", tukey)

# Prepare data for plotting
grouped = df.groupby("Sample")["load"]
means = grouped.mean()
stds = grouped.std()
x_pos = np.arange(len(means))

# Plot
plt.figure(figsize=fig_size_inches, dpi=1200)

# Colors by sample
colors = []
for sample in means.index:
    if "R19" in sample:
        colors.append('lightpink')
    elif "Evolved" in sample:
        colors.append('seagreen')
    elif "Gfp" in sample:
        colors.append('darkseagreen')
    else:
        colors.append('seagreen')

# Bars and points
plt.bar(x_pos, means, yerr=stds, capsize=5, color=colors, width=0.8)
for i, sample in enumerate(means.index):
    y = grouped.get_group(sample)
    jitter = np.random.normal(0, 0.05, size=len(y))
    plt.scatter(np.full_like(y, x_pos[i]) + jitter, y, color='black', s=40, zorder=10)

# Add dashed line at 1
plt.axhline(y=1, color='gray', linestyle='--', linewidth=1)

# Labels and grid
plt.yscale("log")
plt.xticks(x_pos, means.index, rotation=20, ha="right")
plt.ylabel("copy number ratio fungus:bacteria (log10)")
plt.grid(axis='y')
plt.gca().set_axisbelow(True)
plt.tight_layout()
#plt.show()

plt.savefig("Supplementary_Figure_7.png", bbox_inches='tight', dpi=300)