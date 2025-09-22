# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 16:51:20 2025

@author: tgassler
"""

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FixedLocator, PercentFormatter
from matplotlib import rcParams

# --- font & figure params ----------------------------------------------------
font_type = 'Arial'
font_size = 8
rcParams['font.family'] = font_type
rcParams['font.size'] = font_size

def mm_to_inches(mm): return mm / 25.4

# --- plotting helper (unchanged logic) --------------------------------------
def plot_data(ax, df, label, color, show_ylabel):
    # line: average germination rate per round
    avg_germ_rate = df.groupby("Round")['germination rate day 2'].mean()
    avg_germ_rate.plot(ax=ax, kind='line', marker='o', color=color, linestyle='-',
                       linewidth=1, markersize=5, markerfacecolor=color,
                       markeredgewidth=1, alpha=1.0, label=f'{label}')

    # scatter: individual injections/replicates
    injection_numbers = df['Injection'].unique()
    for injection in injection_numbers:
        df_injection = df[df['Injection'] == injection]
        for round_num in df_injection['Round'].unique():
            y_values = df_injection[df_injection['Round'] == round_num]['germination rate day 2']
            ax.scatter([round_num] * len(y_values), y_values, color=color, s=8, alpha=0.6)

    # cosmetics
    ax.set_facecolor('white')
    ax.spines['top'   ].set_visible(True)
    ax.spines['right' ].set_visible(True)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'  ].set_color('black')
    ax.tick_params(axis='both', labelsize=font_size - 2)
    ax.set_xlabel('Round', fontsize=font_size)
    if show_ylabel:
        ax.set_ylabel('', fontsize=font_size)
    ax.set_xlim(0, 21)
    ax.set_ylim(0, 1)
    ax.set_xticks([1, 2, 3, 4, 5, 10, 11, 15, 20])
    ax.grid(True, which='both', axis='y', linestyle='--', linewidth=0.6, color='grey', alpha=0.8)
    ax.legend(fontsize=font_size - 2, loc='upper right')

    # percent y-axis
    yvalues = ax.get_yticks()
    ax.yaxis.set_major_locator(FixedLocator(yvalues))
    ax.yaxis.set_major_formatter(PercentFormatter(1))


df_all = pd.read_excel('sourcedata_Figure_2_b_c_d.xlsx', engine='openpyxl')

# filter exactly as before
df_all = df_all[(df_all['Line'] != '0') & ~df_all['Round'].isin([16, 19])]

# figure out how many injection groups we have (one row per injection)
injections = sorted(df_all['Injection'].dropna().unique().tolist())

# --- figure size (mm → inches) ----------------------------------------------
figure_width_mm  = 180
figure_height_mm = 120
fig_width, fig_height = mm_to_inches(figure_width_mm), mm_to_inches(figure_height_mm)

# --- subplots: rows = injections, cols = 4 categories -----------------------
fig, axes = plt.subplots(len(injections), 4, figsize=(fig_width, fig_height),
                         dpi=1200, sharey=True, sharex=True, squeeze=False)

for i, inj in enumerate(injections):
    df_inj = df_all[df_all['Injection'] == inj]

    # categories (as before)
    df_high_neg = df_inj[df_inj['Line'].str.contains('High Line Negative')]
    df_low_neg  = df_inj[df_inj['Line'].str.contains('Low Line Negative')]
    df_high     = df_inj[(df_inj['Line'] == 'High Line') & (df_inj['germination rate day 2'] >= 0)]
    df_low      = df_inj[(df_inj['Line'] == 'Low Line')  & (df_inj['germination rate day 2'] >= 0)]

    # plot into the 4 panels of this row
    plot_data(axes[i, 0], df_high_neg, 'Negatives from B++', 'red',   i == 0)  # y-label only on top row
    plot_data(axes[i, 1], df_high,     'B++',                'blue',  False)
    plot_data(axes[i, 2], df_low_neg,  'Negatives from B+',  'green', False)
    plot_data(axes[i, 3], df_low,      'B+',                 'orange',False)

# x-labels on bottom row
for ax in axes[-1, :]:
    ax.set_xlabel('Round', fontsize=font_size)

plt.subplots_adjust(wspace=0.3, hspace=0.4)
plt.tight_layout()
#plt.show()
plt.savefig("Supplement_Figure_5_c.png", bbox_inches='tight', dpi=300)
