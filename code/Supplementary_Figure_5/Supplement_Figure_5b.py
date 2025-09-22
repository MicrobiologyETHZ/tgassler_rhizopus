# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 16:22:17 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

# Replace 'your_file.xlsx' with the path to your actual Excel file
file_path = 'source_data_SF_5_b.xlsx'

# Read the Excel file
df = pd.read_excel(file_path)

# Filter the data for IND strains
df_ind = df[df['Strain'] == 'IND']

# Multiply germination rates by 100 to represent as percentages
df_ind['Germination rate day 2'] *= 100

# Set font size globally
plt.rcParams.update({'font.size': 8})

# Set the size of the plot
plt.figure(figsize=(2,2), dpi=600)

# Define the bar width and opacity
bar_width = 0.4
opacity = 0.5

# Define the positions for Day 1, Day 2, and Pooled
positions = [1, 2, 3]

# Plotting bars and single data points with a smoothed line for the distribution
for i, day in enumerate(sorted(df_ind['Day'].unique()) + ['Pooled']):
    # For 'Pooled', combine all data; else, filter the day
    if day == 'Pooled':
        day_data = df_ind['Germination rate day 2']
    else:
        day_data = df_ind[df_ind['Day'] == day]['Germination rate day 2']
    
    if len(day_data) > 0:  # Avoid KDE for empty datasets
        # Calculate the KDE for y-values
        kde = gaussian_kde(day_data)
        kde_positions = np.linspace(day_data.min(), day_data.max(), 100)
        kde_values = kde(kde_positions)
        
        # Scale KDE values to fit within the bar width
        kde_values_scaled = kde_values / max(kde_values) * bar_width
        
        # Align the KDE distribution to the right of the bar
        position = positions[i]
        kde_right_aligned = position + bar_width / 2 - kde_values_scaled
        
        # Plot the KDE as a filled area
        plt.fill_betweenx(kde_positions, kde_right_aligned, position + bar_width / 2, color='gray', alpha=0.5)

    # Plot the scatter points for each day or pooled
    jittered_x = position + np.random.uniform(-bar_width / 2, bar_width / 2, size=len(day_data))
    plt.scatter(jittered_x, day_data, alpha=0.6, color='black')

# Plot transparent bars for the means
means = df_ind.groupby('Day')['Germination rate day 2'].mean().tolist()
means.append(df_ind['Germination rate day 2'].mean())  # Append mean of pooled data
plt.bar(positions, means, color=['blue', 'orange', 'green'], width=bar_width, alpha=opacity, zorder=0)

# Add labels and title
plt.ylabel('Germination success (%)', fontsize=8)
#<plt.title('Germination rate day 2 for IND Strains on Day 1, Day 2, and Pooled', fontsize=8)
plt.xticks(positions, ['Week 1', 'Week 2', 'Pooled'], fontsize=8)
plt.ylim(0,80)
plt.yticks(fontsize=8)

# Show the plot
plt.tight_layout()
#plt.show()
plt.savefig("Supplement_Figure_5_b.png", bbox_inches='tight', dpi=300)
