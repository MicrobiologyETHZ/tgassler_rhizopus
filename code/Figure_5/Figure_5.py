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
import gzip

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

# 1) define your list of IDs (the 35 you combined)
#cell_wall_ids_y = [
 #   "g12504", "g7220", "g7221", "g7219", "g14612", "g18754",
  #  "g17424", "g19292", "g16845", "g22382", "g23080",
   # "g23576", "g24085", "g16271", "g11661", 
    #"g2788", "g14304", "g16391", "g17194", "g7251",
    #"g4445", "g21398", "g20162", "g9186", "g18773", "g16900",
    #"g1304", "g18816", "g1548", "g16312", "g21370"
#]

cell_wall_ids_y = [
    "g12504", "g7220", "g7221", "g7219", "g14612", "g18754",
    "g19292", "g16845", "g22382", "g23080",
    "g23576", "g24085", "g16271", "g11661", 
    "g2788", "g14304", "g17194", "g7251",
    "g21398", "g20162", "g18773", "g16900",
    "g1304", "g18816", "g1548", "g16312", "g21370"
]

cell_wall_ids_n = [
    "g7207", "g5475", "g5425","g19047","g9180"
]

# now your second list
ROS_ids_y = ["g9145"]
ROS_ids_n = ["g17875","g9311","g17858"]

# Load the dataset
# Please replace '/path/to/your/IND_B3vsB1.xlsx' with the actual path to your file
df = pd.read_csv("2025-05-19_24TG06_Rp_R20_pos_vs_Rp_R10_pos_l0a0.01_results.csv.gz")

# %% Choosing of FDR

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=300)
plt.hist(df['pvalue'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('p-value')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig("Figure_5_x.png")

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=300)
plt.hist(df['padj'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('padj')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()


# %% Calculate absolute log2 fold change and -log10 of p-value

df['abs_log2FoldChange'] = abs(df['log2FoldChange'])
df['-log10(pvalue)'] = -np.log10(df['pvalue'])

# %% Create the volcano plot for different significance levels for up or downregulated genes
plt.figure(figsize=(figure_width_in, figure_height_in), dpi=2400)

# All points
plt.scatter(df['log2FoldChange'], -np.log10(df['pvalue']), edgecolor='grey', facecolors='none', alpha=0.5, s=7)

# Significance criteria for upregulated genes
significant_up = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)
plt.scatter(df[significant_up]['log2FoldChange'], -np.log10(df[significant_up]['pvalue']), edgecolor='#1F4FA0', facecolors='none', alpha=0.5, s=7)

# Significance criteria for downregulated genes
significant_down = (df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)
plt.scatter(df[significant_down]['log2FoldChange'], -np.log10(df[significant_down]['pvalue']), edgecolor='#CE3B6C', facecolors='none', alpha=0.5, s=7)

# now highlight the cell wall genes that follow the trend
cell_wall_ids_y = df['Geneid'].isin(cell_wall_ids_y)
hl = df[cell_wall_ids_y]
plt.scatter(hl['log2FoldChange'], -np.log10(hl['pvalue']),
            c='green', edgecolor='none', s=8, alpha=0.8) # label='your genes')

# now highlight the cell wall genes that do not follow the trend
cell_wall_ids_n = df['Geneid'].isin(cell_wall_ids_n)
hl2 = df[cell_wall_ids_n]
plt.scatter(hl2['log2FoldChange'], -np.log10(hl2['pvalue']),
            c='none', edgecolor='green', s=8, alpha=0.8, linewidths=0.4) # label='your genes')


ROS_ids_y = df['Geneid'].isin(ROS_ids_y)
hl3 = df[ROS_ids_y]
plt.scatter(hl3['log2FoldChange'], -np.log10(hl3['pvalue']),
            c='yellow', edgecolor='none', s=8, alpha=0.8)

ROS_ids_n = df['Geneid'].isin(ROS_ids_n)
hl4 = df[ROS_ids_n]
plt.scatter(hl4['log2FoldChange'], -np.log10(hl4['pvalue']),
            c='none', edgecolor='yellow', s=8, alpha=0.8, linewidths=0.4)



# Labels and Axes limits
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
plt.tight_layout()
#plt.show()
plt.savefig("Figure_5_b.png", bbox_inches='tight', dpi=300)

# %% Filter for DEGs with fold change > 1 and padj < 0.01
deg_criteria = (df['padj'] < 0.05) & (df['abs_log2FoldChange'] > 0.6)
deg_df = df[deg_criteria]

# Upregulated genes
deg_up_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)]

# Downregulated genes
deg_down_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)]

# %% Plotting the nDEGs

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
#plt.show()


#%% Convert target to csv to get the target list for GO term enrichment

deg_up_df.to_csv('Rp_R20_Pos_VS_R10_Pos_up.csv', index=False)

deg_down_df.to_csv('Rp_R20_Pos_VS_R10_Pos_down.csv', index=False)

# Trial to link to GO term analysis via string

# Define the paths to your files
csv_path= 'Rp_R20_Pos_VS_R10_Pos_up.csv'
fasta_path = "../Figure_4/braker.faa.gz"
output_fasta_path = "Rp_R20_Pos_VS_R10_Pos_01_up.fasta"

# Read the Excel file
df = pd.read_csv(csv_path)

# Parse the FASTA file
fasta_sequences = {}
with gzip.open(fasta_path, 'rt') as fasta_handle:
    for record in SeqIO.parse(fasta_handle, "fasta"):
        header_parts = record.id.split(';')[0].split('=')
        if len(header_parts) > 1:
            gene_id = header_parts[1]
        else:
            # Handle the case where gene_id is not present
            gene_id = None  # or some default value
        isoform = record.id.split(';')[0]
        fasta_sequences[isoform] = str(record.seq)


# Prepare a list for the SeqRecord objects
target_sequences = []

# Iterate over each row in the DataFrame
for _, row in df.iterrows():
    gene_id = row['Geneid']
    
    # Get the .t1 and .t2 sequences if they exist
    seq_t1 = fasta_sequences.get(f"{gene_id}.t1")
    seq_t2 = fasta_sequences.get(f"{gene_id}.t2")

    # Create SeqRecord objects and add them to the list
    if seq_t1:
        target_sequences.append(SeqRecord(Seq(seq_t1), id=f"{gene_id}.t1", description=""))
    if seq_t2:
        target_sequences.append(SeqRecord(Seq(seq_t2), id=f"{gene_id}.t2", description=""))

# Write the target sequences to a new FASTA file
with open(output_fasta_path, 'w') as output_handle:
    SeqIO.write(target_sequences, output_handle, 'fasta')


# Submit output files against our proteome STRG0A18OCV

