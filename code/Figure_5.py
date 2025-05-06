# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 14:55:07 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.colors as mcolors
import sys
from utilities import load_figure_data, save_figure_panel, get_deg_gene_seqs

# Define inputs for figure size, font type, and font size
figure_width_mm = 50  # Width in mm
figure_height_mm = 50  # Height in mm
font_type = "Arial"
font_size = 8

# Convert figure size from mm to inches (1 inch = 25.4 mm)
figure_width_in = figure_width_mm / 25.4
figure_height_in = figure_height_mm / 25.4

# Set font properties globally for all plots
plt.rcParams.update({
    "font.family": font_type,
    "font.size": font_size
})

# Load the dataset
# Please replace '/path/to/your/IND_B3vsB1.xlsx' with the actual path to your file
#df = pd.read_csv("2025-01-24_24TG06_Rp_R20_pos_vs_Rp_R10_pos_l0a0.01_results.csv")

df = load_figure_data("Figure_5", "panel_b")

# p-value distributions

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=300)
plt.hist(df['pvalue'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('p-value')
plt.ylabel('Number of Genes')
plt.grid(True)

# plt.show()


plt.figure(figsize=(figure_width_in, figure_height_in), dpi=300)
plt.hist(df['padj'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('padj')
plt.ylabel('Number of Genes')
plt.grid(True)
# plt.show()

##############
# Panel 5 b #
#############

# Calculate absolute log2 fold change and -log10 of p-value

df['abs_log2FoldChange'] = abs(df['log2FoldChange'])
df['-log10(pvalue)'] = -np.log10(df['pvalue'])

# Create the volcano plot for different significance levels for up or downregulated genes
plt.figure(figsize=(figure_width_in, figure_height_in), dpi=2400)

# All points
plt.scatter(df['log2FoldChange'], -np.log10(df['pvalue']), edgecolor='grey', facecolors='none', alpha=0.5, s=7)

# Significance criteria for upregulated genes
significant_up = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)
plt.scatter(df[significant_up]['log2FoldChange'], -np.log10(df[significant_up]['pvalue']), edgecolor='#1F4FA0', facecolors='none', alpha=0.5, s=7)

# Significance criteria for downregulated genes
significant_down = (df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)
plt.scatter(df[significant_down]['log2FoldChange'], -np.log10(df[significant_down]['pvalue']), edgecolor='#CE3B6C', facecolors='none', alpha=0.5, s=7)

# Labels and Axes limits
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
# plt.show()

save_figure_panel("Figure_5", "panel_b",  format='png')


#  Filter for DEGs with fold change > 1 and padj < 0.01
deg_criteria = (df['padj'] < 0.05) & (df['abs_log2FoldChange'] > 0.6)
deg_df = df[deg_criteria]

# Upregulated genes
deg_up_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)]

# Downregulated genes
deg_down_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)]

######################
# Plotting the nDEGs #
######################

# Calculate the number of DEGs for each category
num_deg_up = len(deg_up_df)
num_deg_down = len(deg_down_df)

# Data for plotting
categories = ['DEGs']
counts_up = [num_deg_up]
counts_down = [num_deg_down]

# Create the bar plot as a stacked bar plot
plt.figure(figsize=(figure_width_in, figure_height_in), dpi=300)
bars_up = plt.bar(categories, counts_up, color='#0072B2', alpha=0.5, label='DEGs Up')
bars_down = plt.bar(categories, counts_down, color='#D55E00', alpha=0.5, bottom=counts_up, label='DEGs Down')

# Adding the counts above the bars for the total
total_counts = [x + y for x, y in zip(counts_up, counts_down)]
for idx, val in enumerate(total_counts):
    plt.text(idx, val + 0.05, val, ha='center', va='bottom')

# Adding the counts above the individual bars
for bar in bars_up:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, yval / 2, yval, ha='center', va='bottom')

for bar in bars_down:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, yval + counts_up[0] - yval / 2, yval, ha='center', va='bottom')

# Adding labels and title
plt.ylabel('nDEG')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout(rect=[0, 0, 0.85, 1])
# plt.show()

# Convert target to csv to get the target list for GO term enrichment

# deg_up_df.to_csv('Rp_R20_Pos_VS_R10_Pos_up.csv', index=False)

# deg_down_df.to_csv('Rp_R20_Pos_VS_R10_Pos_down.csv', index=False)

# Link to GO term analysis via String

get_deg_gene_seqs(deg_up_df, output_fasta_path="Rp_R20_Pos_VS_R10_Pos_up.fasta")
get_deg_gene_seqs(deg_down_df, output_fasta_path="Rp_R20_Pos_VS_R10_Pos_down.fasta")

# Submit output files against our proteome STRG0A18OCV

