# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 13:32:06 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# ========== USER SETTINGS FOR STYLING ==========
width_mm = 50     # figure width in mm
height_mm = 50    # figure height in mm
font_size = 8     # desired font size
font_family = "Arial"
dpi_val = 1200

# Convert mm to inches
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Set global font properties
mpl.rcParams['font.family'] = font_family
mpl.rcParams['font.size'] = font_size

#############
# Panel 2 e #
#############

# Load the data
file_path = r'sourcedata_Figure_2_e_f_data.xlsx'
df = pd.read_excel(file_path)

# Clean up columns and convert columns with '%' to numeric
df.columns = df.columns.str.strip()

# Extract day numbers and sort
df['Day_num'] = df['Day'].str.extract(r'(\d+)').astype(int)
df = df.sort_values(by='Day_num')

# Filter data
df_r02_pos = df[df['Sample'].str.contains('R02') & df['Sample'].str.contains('_pos')]
df_r02_neg = df[df['Sample'].str.contains('R02') & df['Sample'].str.contains('_neg')]
df_r20_pos = df[df['Sample'].str.contains('R20') & df['Sample'].str.contains('_pos')]
df_r20_neg = df[df['Sample'].str.contains('R20') & df['Sample'].str.contains('_neg')]

# Group means for Germination Rate
r02_pos_means = df_r02_pos.groupby('Day_num')['Germination Rate Plates'].mean()
r02_neg_means = df_r02_neg.groupby('Day_num')['Germination Rate Plates'].mean()
r20_pos_means = df_r20_pos.groupby('Day_num')['Germination Rate Plates'].mean()
r20_neg_means = df_r20_neg.groupby('Day_num')['Germination Rate Plates'].mean()

# All unique days
all_days = sorted(df['Day_num'].unique())

# Create main figure without legend
fig, ax1 = plt.subplots(figsize=(width_in, height_in), dpi=dpi_val)

width = 0.2
x_positions = range(len(all_days))

# Bar plots without error bars
ax1.bar([x - width * 1.5 for x in x_positions], r02_pos_means.reindex(all_days, fill_value=0).values * 100,
        width, label='R2 B$_{pos}$', color="#3D296D", alpha=0.7)
ax1.bar([x - width / 2 for x in x_positions], r02_neg_means.reindex(all_days, fill_value=0).values * 100,
        width, label='R2 B$_{neg}$', color='#9E94B6', alpha=0.7)
ax1.bar([x + width / 2 for x in x_positions], r20_pos_means.reindex(all_days, fill_value=0).values * 100,
        width, label='R20 B$_{pos}$', color="#FDCF92", alpha=0.7)
ax1.bar([x + width * 1.5 for x in x_positions], r20_neg_means.reindex(all_days, fill_value=0).values * 100,
        width, label='R20 B$_{neg}$', color='#FFE7C1', alpha=0.7)

ax1.set_yscale('symlog')  # Set y-axis to logarithmic scale
ax1.set_ylim(1, 100)   # Set limits appropriately (avoid zero for log-scale)
ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())  # Maintain percentage format

ax1.set_xlabel('Day', fontsize=font_size)
ax1.set_ylabel('Germination Rate (%)', fontsize=font_size)
ax1.tick_params(axis='both', which='major', labelsize=font_size)
ax1.set_ylim(0, 100)
ax1.set_xticks(x_positions)
ax1.set_xticklabels(all_days)

ax1.grid(False)
plt.tight_layout()
#plt.show()
fig.savefig("Figure_2_e.png", bbox_inches='tight', dpi=300)


# Create separate legend figure
fig_leg = plt.figure(figsize=(1, 1), dpi=dpi_val)
handles, labels = ax1.get_legend_handles_labels()
fig_leg.legend(handles, labels, loc='center', frameon=False, fontsize=font_size, ncol=1)
fig_leg.canvas.draw()
plt.axis('off')
#plt.show()
fig_leg.savefig("Figure_2_e_legend.png", bbox_inches='tight', dpi=300)



#############
# Panel 2 f #
#############


df['True Positives (N)'] = df['True Positives (N)'].fillna(0)
df['True Positives (Pos)'] = df['True Positives (Pos)'].fillna(0)

# 3. Extract numeric day, define strain
df['Day_num'] = df['Day'].str.extract(r'(\d+)').astype(int)
df['Strain'] = df['Sample'].apply(
    lambda x: 'R02' if 'R02' in x else ('R20' if 'R20' in x else 'Other')
)

# 4. Group by (Day_num, Strain)
grouped = df.groupby(['Day_num', 'Strain'], as_index=False).agg(
    sumN=('True Positives (N)', 'sum'),     # Number of positives
    sumPos=('True Positives (Pos)', 'sum')  # Total tested
)

# 5. Calculate fractions (avoid divide-by-zero)
grouped['frac_pos'] = np.where(grouped['sumPos'] == 0,
                               0,
                               grouped['sumN'] / grouped['sumPos'])
grouped['frac_neg'] = 1 - grouped['frac_pos']

# 6. Pivot to separate R02 and R20 columns
pivoted = grouped.pivot(index='Day_num', columns='Strain', values=['sumPos','sumN','frac_pos','frac_neg'])
pivoted = pivoted.sort_index()

days = pivoted.index.values

# For R02
r02_frac_pos = pivoted[('frac_pos','R02')].fillna(0).values
r02_frac_neg = pivoted[('frac_neg','R02')].fillna(0).values

# For R20
r20_frac_pos = pivoted[('frac_pos','R20')].fillna(0).values
r20_frac_neg = pivoted[('frac_neg','R20')].fillna(0).values

# ========== PLOTTING ==========
fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=dpi_val)

x = np.arange(len(days))
width = 0.4
offset = width / 2

# === R02 bars (Positive on bottom in black; Negative on top in light gray) ===
ax.bar(
    x - offset,
    r02_frac_pos * 100,
    width,
    label='R2 Bacteria retained',
    color='#3D296D',
    alpha=0.8
)
ax.bar(
    x - offset,
    r02_frac_neg * 100,
    width,
    bottom=r02_frac_pos * 100,
    label='R2 Bacteria cleared',
    color='lightgrey',
    alpha=0.6
)

# === R20 bars (Positive on bottom in green; Negative on top in light green) ===
ax.bar(
    x + offset,
    r20_frac_pos * 100,
    width,
    label='R20 Bacteria retained',
    color='#FDCF92',
    alpha=0.8
)
ax.bar(
    x + offset,
    r20_frac_neg * 100,
    width,
    bottom=r20_frac_pos * 100,
    label='R20 Bacteria cleared',
    color='grey',
    alpha=0.6
)

# Axis labeling
ax.set_xlabel('Day', fontsize=font_size)
ax.set_ylabel('Positive Rate (%)', fontsize=font_size)
ax.set_xticks(x)
ax.set_xticklabels(days, fontsize=font_size)
ax.set_ylim(0, 110)  # a bit above 100% for clarity
ax.tick_params(axis='both', which='major', labelsize=font_size)

# ========== Legend as a separate figure ==========
legend_fig = plt.figure(figsize=(2, 1), dpi=dpi_val)
legend_ax = legend_fig.add_subplot(111)
legend_ax.axis('off')
handles, labels = ax.get_legend_handles_labels()
legend_ax.legend(handles, labels, loc='center', fontsize=font_size, frameon=False)
#legend_fig.savefig('legend_only.png', bbox_inches='tight', pad_inches=0.1)

# Remove legend from main plot
ax.legend().remove()

plt.tight_layout()
#plt.show()
fig.savefig("Figure_2_f.png", bbox_inches='tight', dpi=300)


#R02B-= #9E94B6; R02B+ = #3D296D; R20+=#FDCF92; R20- = #FFE7C1