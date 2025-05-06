# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:33:43 2025

@author: tgassler
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from utilities import load_figure_data, save_figure_panel
# ========== USER SETTINGS FOR STYLING ==========
width_mm = 40   # figure width in mm
height_mm = 60  # figure height in mm
font_size = 8   # desired font size
font_family = "Arial"

# Convert mm -> inches
width_in = width_mm / 25.4
height_in = height_mm / 25.4

mpl.rcParams['font.family'] = font_family
plt.rcParams.update({'font.size': font_size})

# ========== READ THE CONSOLIDATED DATA ========== 

df_results = load_figure_data("Figure_3", "panel_c")

# ========== DEFINE CONDITION ORDER & COLORS ==========
conditions_list = ["NH strain", "R2 B+", "R11 B+", "R20 B+"]
colors = plt.cm.tab20.colors
color_map = {cond: colors[i] for i, cond in enumerate(conditions_list)}

# Create a color array for each row in df_results, default to 'gray' if not in color_map
bar_colors = [color_map.get(cond, 'gray') for cond in df_results['condition']]

# ========== CREATE THE BAR PLOT ==========
plt.figure(figsize=(width_in, height_in), dpi=1200)
df_results['condition'] = df_results['condition'].str.replace(" B+", "")
plt.bar(df_results['condition'], df_results['germination_rate'], color=bar_colors)
plt.ylabel('Germination Rate (%)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
save_figure_panel("Figure_3", "panel_c", format='png')
# plt.show()
