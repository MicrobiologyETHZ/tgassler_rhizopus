# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 14:40:07 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl

# ========= USER SETTINGS FOR STYLING =========
width_mm = 40    # main plots: 45 mm wide
height_mm = 60   # main plots: 70 mm tall
legend_width_mm = 80  # legend: 80 mm wide
legend_height_mm = 10 # legend: 30 mm tall

font_size = 8
font_family = "Arial"
dpi = 2400

# Convert mm to inches
width_in = width_mm / 25.4
height_in = height_mm / 25.4
legend_width_in = legend_width_mm / 25.4
legend_height_in = legend_height_mm / 25.4

# Set global font & size
mpl.rcParams['font.family'] = font_family
plt.rcParams.update({'font.size': font_size})

# ========== EXAMPLE DATA & FILTERING (Use Your Actual Data) ==========
df = pd.read_excel('sourcedata_Figure_2_b_c_d.xlsx', engine='openpyxl')



df['Abundance'] *= 100
df['Fitness Index'] *= 100
df['Normalized Germ Rate'] *= 100

marker_dict = {1: 'o', 2: '^', 3: '*'}
palette = sns.color_palette("magma", len(df['Round'].unique()))

sns.set(style="ticks")
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 0.8

def add_full_border(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.8)

def plot_custom_bars(ax, df, y_column, ylabel, palette, lines=['B++', 'B+']):
    rounds = sorted(df['Round'].unique())
    bar_width = 0.05
    group_gap = 0.05
    x_positions = np.arange(len(lines)) * (0.5 + group_gap)

    round_colors = {rnd: palette[i % len(palette)] for i, rnd in enumerate(rounds)}
    offsets = np.linspace(-bar_width * len(rounds) / 2, bar_width * len(rounds) / 2, len(rounds))

    for i, rnd in enumerate(rounds):
        round_data = df[df['Round'] == rnd]
        for j, line in enumerate(['High Line', 'Low Line']):
            line_data = round_data[round_data['Line'] == line]
            x_pos = x_positions[j] + offsets[i]

            ax.bar(
                x_pos, line_data[y_column].mean(),
                width=bar_width, color=round_colors[rnd],
                edgecolor=None, linewidth=0, alpha=0.8, rasterized=True
            )
            # Overlay data points
            for injection, marker in marker_dict.items():
                injection_data = line_data[line_data['Injection'] == injection]
                ax.scatter(
                    [x_pos]*len(injection_data),
                    injection_data[y_column].values,
                    color='black', marker=marker,
                    edgecolor='none', alpha=0.7, s=8
                )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(lines, fontsize=font_size, rotation=45, ha='right')
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)
    add_full_border(ax)

# ========== MAIN PLOT WITH LEGEND ==========
fig1, ax1 = plt.subplots(figsize=(width_in + legend_width_in, height_in), dpi=dpi)
plot_custom_bars(ax1, df, 'Normalized Germ Rate', 'Normalized germination success (%)', palette)

# Create legend handles
rounds = sorted(df['Round'].unique())
handles = [
    plt.Line2D([0], [0], color=palette[i], lw=3, label=f'Round {rnd}')
    for i, rnd in enumerate(rounds)
]
marker_handles = [
    plt.Line2D([0], [0], color='black', marker=marker, linestyle='None', markersize=8, label=f'Injection {inj}')
    for inj, marker in marker_dict.items()
]
legend_handles = handles + marker_handles

# Add legend to the right of the plot
ax1.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fontsize=font_size,
    title_fontsize=font_size,
    frameon=False,
    labelspacing=0.4,
    handletextpad=0.3,
    borderpad=0.1,
    borderaxespad=0.4
)

plt.tight_layout()
#fig1.savefig("normalized_germination_rate_with_legend.png", bbox_inches='tight')
#plt.show()
plt.savefig("Supplementary_Figure_6.png", bbox_inches='tight', dpi=300)
