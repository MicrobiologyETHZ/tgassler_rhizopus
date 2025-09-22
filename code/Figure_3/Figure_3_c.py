# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:33:43 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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
df_results = pd.read_excel("Figure_3_c.xlsx", sheet_name="Figure_3_c", engine="openpyxl")

# ========== OPTIONAL: PRINT THE LOADED DATA ==========
print("Loaded Data for Plotting:")
print(df_results)

# ========== DEFINE CONDITION ORDER & COLORS ==========
conditions_list = ["NH strain", "R2 B+", "R11 B+", "R20 B+"]
colors = plt.cm.tab20.colors
color_map = {cond: colors[i] for i, cond in enumerate(conditions_list)}

# Create a color array for each row in df_results, default to 'gray' if not in color_map
bar_colors = [color_map.get(cond, 'gray') for cond in df_results['condition']]

# ========== CREATE THE BAR PLOT ==========
plt.figure(figsize=(width_in, height_in), dpi=1200)
plt.bar(df_results['condition'], df_results['germination_rate'], color=bar_colors)

plt.ylabel('Germination Rate (%)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

#plt.show()
plt.savefig("Figure_3_c.png",  bbox_inches='tight', dpi=300)
