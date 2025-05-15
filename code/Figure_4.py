# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 08:37:54 2025

@author: tgassler
"""
#########################
# Figure 4 a: PCA plots #
# #######################
  
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.transforms as transforms
from matplotlib_venn import venn3
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from matplotlib.lines import Line2D

from utilities import load_figure_data, save_figure_panel


def confidence_ellipse(x, y, ax, n_std=1.5, facecolor='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
                      width=ell_radius_x * 2,
                      height=ell_radius_y * 2,
                      facecolor=facecolor,
                      **kwargs)

    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


# Load the metadata
metadata = load_figure_data("Figure_4", "metadata", index_col=0)

# Load the normalized counts data
counts = load_figure_data("Figure_4", "panel_a", index_col=0)

# Transpose counts DataFrame to match the sample IDs in the metadata
counts = counts.set_index(counts.columns[0]).T.reset_index()
counts.columns = ['bmk_id'] + list(counts.columns[1:])

# Merge the metadata with the counts data on 'bmk_id'
merged_data = pd.merge(metadata, counts, on='bmk_id')
# Define the conditions to include
conditions_to_include = ['Rp_R10_pos', 'Rp_R20_pos', 'pIND']
# Include as needed 'pIND': 'NH strain','Rp_R10_pos': 'R10 B+','Rp_R10_neg': 'R10 B-',
# 'Rp_R20_pos': 'R20 B+','Rp_R20_neg': 'R20 B-' 

# Filter the merged data to include only the desired conditions
filtered_data = merged_data[merged_data['Condition'].isin(conditions_to_include)].copy()

# Map old condition names to new ones
condition_mapping = {
    'pIND': 'NH strain',
    'Rp_R10_pos': 'R10 B$_{pos}$',
    #'Rp_R10_neg': 'R10 B-',
    'Rp_R20_pos': 'R20 B$_{pos}$',
    #'Rp_R20_neg': 'R20 B-'
}

# Apply the mapping to the 'Condition' column
filtered_data['Condition'] = filtered_data['Condition'].map(condition_mapping)

# Select the count data for PCA
count_data = filtered_data.iloc[:, 8:]  # assuming first 8 columns are metadata

# Perform PCA
pca = PCA(n_components=6)
pca_result = pca.fit_transform(count_data)
filtered_data['PCA1'] = pca_result[:, 0]
filtered_data['PCA2'] = pca_result[:, 1]

# Calculate explained variance ratio for effect size
explained_variance = pca.explained_variance_ratio_

# Calculate the standard deviation for each principal component
std_devs = pca.explained_variance_ ** 0.5

# Define a color palette for the conditions
palette = sns.color_palette("colorblind", len(conditions_to_include))
color_map = dict(zip(condition_mapping.values(), palette))

##############
# Panel 4 a #
##############

# Define inputs for figure size, font type, and font size
figure_width_mm = 80  # Width in mm
figure_height_mm = 50  # Height in mm
font_type = "Arial"
font_size = 8

figure_width_in = figure_width_mm / 25.4
figure_height_in = figure_height_mm / 25.4

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=1200)
plt.rcParams.update({
    "font.family": font_type,
    "font.size": font_size
})

# Scatterplot with consistent colors
ax = sns.scatterplot(
    x='PCA1', y='PCA2',
    hue='Condition',
    palette=color_map,  # Use color map for consistent coloring
    data=filtered_data,
    legend="full"
)

# Calculate and plot ellipses for each condition
for condition in condition_mapping.values():
    subset = filtered_data[filtered_data['Condition'] == condition]
    confidence_ellipse(subset['PCA1'].values, subset['PCA2'].values, ax,
                       n_std=1.5, facecolor=color_map[condition],
                       edgecolor=color_map[condition], alpha=0.2)

# Labeling
plt.xlabel(f'PCA1 ({explained_variance[0]:.1%})')
plt.ylabel(f'PCA2 ({explained_variance[1]:.1%})')

# Place the legend right next to the plot
plt.legend(
    loc="center left", 
    bbox_to_anchor=(1, 0.5),  # Place legend to the right of the plot
    title="Condition",
    frameon=False  # Optional: remove legend frame for a cleaner look
)

# Adjust layout to ensure everything fits well
plt.tight_layout()
save_figure_panel("Figure_4", "panel_a", format='png')
# plt.show()

#############################
# Figure 4 b: Venn Diagramm #
#############################



# Load your dataframes
deg_RP1 = load_figure_data("Figure_4", "panel_b_part1")
deg_RP2 = load_figure_data("Figure_4", "panel_b_part2")

# Extract the "Geneid" columns and convert them to sets for all nDEGs
set_RP1 = set(deg_RP1['Geneid'])
set_RP2 = set(deg_RP2['Geneid'])

# Filter upregulated genes (Log2FoldChange > 0)
upregulated_RP1 = set(deg_RP1[deg_RP1['log2FoldChange'] > 0]['Geneid'])
upregulated_RP2 = set(deg_RP2[deg_RP2['log2FoldChange'] > 0]['Geneid'])

# Filter downregulated genes (Log2FoldChange < 0)
downregulated_RP1 = set(deg_RP1[deg_RP1['log2FoldChange'] < 0]['Geneid'])
downregulated_RP2 = set(deg_RP2[deg_RP2['log2FoldChange'] < 0]['Geneid'])

# Define colors to match PCA plot
# color_map = {
#     'pIND': "#2ca02c",       # Green
#     'Rp_R10_pos': "#1f77b4", # Blue
#     'Rp_R20_pos': "#ff7f0e"  # Orange
# }


color_map = {
    'pIND': "#0072B2",       # Green
    'Rp_R10_pos': "#E69F00", # Blue
    'Rp_R20_pos': "#009E73"  # Orange
}

set_colors = [color_map['Rp_R10_pos'], color_map['Rp_R20_pos'], color_map['pIND']]

# Create a figure with GridSpec layout
fig = plt.figure(figsize=(3.5, 2.1), dpi=600)
plt.rcParams.update({'font.size': 6})
gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1])

# Venn diagram for all nDEGs
ax0 = fig.add_subplot(gs[:, 0])
venn3([set_RP1, set_RP2, set()], ('RP_R10', 'RP_R20', ''), ax=ax0, set_colors=set_colors)
ax0.set_title("All DEGs")

# Venn diagram for upregulated genes
ax1 = fig.add_subplot(gs[0, 1])
venn3([upregulated_RP1, upregulated_RP2, set()], ('RP_R10', 'RP_R20', ''), ax=ax1, set_colors=set_colors)
ax1.set_title("Upregulated DEGs")

# Venn diagram for downregulated genes
ax2 = fig.add_subplot(gs[1, 1])
venn3([downregulated_RP1, downregulated_RP2, set()], ('RP_R10', 'RP_R20', ''), ax=ax2, set_colors=set_colors)
ax2.set_title("Downregulated DEGs")

# Adjust layout
plt.tight_layout()
save_figure_panel("Figure_4", "panel_b", format='png')
# plt.show()

###########################
# Figure 4 c: upregulated #
###########################

# Identify overlapping upregulated genes
overlapping_genes = upregulated_RP1.intersection(upregulated_RP2)

# Filter rows for the overlapping genes in both dataframes
overlap_RP1 = deg_RP1[deg_RP1['Geneid'].isin(overlapping_genes)]
overlap_RP2 = deg_RP2[deg_RP2['Geneid'].isin(overlapping_genes)]

# Merge dataframes to include details from both deg_RP1 and deg_RP2
overlap_df = pd.merge(
    overlap_RP1,
    overlap_RP2,
    on='Geneid',
    suffixes=('_RP1', '_RP2')
)

# Save the result to a CSV file for further analysis (optional)
# overlap_df.to_csv("overlapping_upregulated_genes.csv", index=False)


# Define parameters for figure size (in mm), font type, and font size
heatmap_width_mm = 15  # Heatmap width in mm
heatmap_height_mm = 75  # Heatmap height in mm
legend_width_mm = 8  # Legend width in mm
legend_height_mm = 2  # Legend height in mm
font_type = "Arial"  # Font type
font_size = 6  # Font size

# Convert figure size from mm to inches (1 inch = 25.4 mm)
heatmap_width_in = heatmap_width_mm / 25.4
heatmap_height_in = heatmap_height_mm / 25.4
legend_width_in = legend_width_mm / 25.4
legend_height_in = legend_height_mm / 25.4

# Update font type and size globally
plt.rcParams.update({"font.family": font_type, "font.size": font_size})

# Cap extreme values for better scaling in the heatmap
overlap_df["log2FoldChange_RP1"] = overlap_df["log2FoldChange_RP1"].clip(0, 5)
overlap_df["log2FoldChange_RP2"] = overlap_df["log2FoldChange_RP2"].clip(0, 5)

# Sort by log2FoldChange in RP1
overlap_df = overlap_df.sort_values("log2FoldChange_RP1", ascending=False)

# Create a heatmap dataframe with Geneid as the index
heatmap_data = overlap_df[["Geneid", "log2FoldChange_RP1", "log2FoldChange_RP2"]]
heatmap_data = heatmap_data.set_index("Geneid")

# Create the heatmap without annotations and without the color bar
plt.figure(figsize=(heatmap_width_in, heatmap_height_in), dpi=1200)
sns.heatmap(
    heatmap_data,
    annot=False,  # Remove annotations
    cmap="YlOrRd",  # Yellow to red color scale
    cbar=False,  # Remove the color bar
    linewidths=0.1,
    yticklabels=False,
    xticklabels=False
)
# plt.title("Heatmap of Overlapping Upregulated Genes (Ranked by Fold Change in RP1)")
plt.xlabel("Condition")
plt.ylabel("Genes (Ranked by Fold Change of R10 B+ vs. NH strain)")
plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_up", format="png")
# plt.show()

# Create the standalone colorbar
# Define the range for log2FoldChange (0 to 5)
vmin, vmax = 0, 5

# Create a standalone colorbar
fig, ax = plt.subplots(figsize=(legend_width_in, legend_height_in), dpi=1200)  # Use input sizes
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = plt.cm.YlOrRd  # Yellow to red color scale

# Add the colorbar to the figure
cb = mpl.colorbar.ColorbarBase(
    ax=ax,
    cmap=cmap,
    norm=norm,
    orientation="horizontal"
)

# Save or display the standalone colorbar
plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_up_colorbar", format="png")
# plt.show()

#######################################
# Figure 4 c: upregualted top10 genes # 
#######################################

# Define parameters for figure size (in mm), font type, and font size
heatmap_width_mm = 15  # Heatmap width in mm
heatmap_height_mm = 37.5  # Heatmap height in mm
legend_width_mm = 8  # Legend width in mm
legend_height_mm = 2  # Legend height in mm
font_type = "Arial"  # Font type
font_size = 6  # Font size

# Convert figure size from mm to inches (1 inch = 25.4 mm)
heatmap_width_in = heatmap_width_mm / 25.4
heatmap_height_in = heatmap_height_mm / 25.4
legend_width_in = legend_width_mm / 25.4
legend_height_in = legend_height_mm / 25.4

# Update font type and size globally
plt.rcParams.update({"font.family": font_type, "font.size": font_size})

# Select the top 25% based on log2FoldChange_RP1
# top_25_percent = overlap_df.head(len(overlap_df) // 10)

top_10_genes = overlap_df.nlargest(10, "log2FoldChange_RP1")
# Create a heatmap dataframe with Geneid as the index
heatmap_data = top_10_genes[['Geneid', 'log2FoldChange_RP1', 'log2FoldChange_RP2']]
heatmap_data = heatmap_data.set_index('Geneid')

# Create the heatmap without annotations
plt.figure(figsize=(heatmap_width_in, heatmap_height_in), dpi=1200)
sns.heatmap(
    heatmap_data,
    annot=False,  # Remove annotations
    cmap='YlOrRd',  # Yellow to red color scale
    cbar=False,  # Remove the color bar
    linewidths=0.25,
    yticklabels=False,
    xticklabels=False
)
# plt.title('Heatmap of Top 25% Overlapping Upregulated Genes (Ranked by Fold Change in RP1)')
plt.xlabel('Condition')
plt.ylabel('')
plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_up_top10", format="png")
# plt.show()

# Plot the legend as a separate figure
vmin, vmax = 0, 5  # Define the range for log2FoldChange

fig, ax = plt.subplots(figsize=(legend_width_in, legend_height_in), dpi=1200)  # Use input sizes
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = plt.cm.YlOrRd  # Yellow to red color scale

# Add the colorbar to the figure
cb = mpl.colorbar.ColorbarBase(
    ax=ax,
    cmap=cmap,
    norm=norm,
    orientation="horizontal",
    label=""
)

plt.tight_layout()
# save_figure_panel("Figure_4", "panel_c_up_legend", format="png")
# plt.show()

#############################
# Fogure 4 c: downregulated #
#############################

# Identify overlapping downregulated genes
overlapping_genes = downregulated_RP1.intersection(downregulated_RP2)

# Filter rows for the overlapping genes in both dataframes
overlap_RP1 = deg_RP1[deg_RP1['Geneid'].isin(overlapping_genes)]
overlap_RP2 = deg_RP2[deg_RP2['Geneid'].isin(overlapping_genes)]

# Merge dataframes to include details from both deg_RP1 and deg_RP2
overlap_df = pd.merge(
    overlap_RP1,
    overlap_RP2,
    on='Geneid',
    suffixes=('_RP1', '_RP2')
)

# Cap extreme values for better scaling in the heatmap
overlap_df["log2FoldChange_RP1"] = overlap_df["log2FoldChange_RP1"].clip(-5, 0)
overlap_df["log2FoldChange_RP2"] = overlap_df["log2FoldChange_RP2"].clip(-5, 0)

# Sort by log2FoldChange in RP1
overlap_df = overlap_df.sort_values("log2FoldChange_RP1")


# Save the result to a CSV file for further analysis (optional)
# overlap_df.to_csv("overlapping_downregulated_genes.csv", index=False)

# Define parameters for figure size (in mm), font type, and font size
heatmap_width_mm = 15  # Heatmap width in mm
heatmap_height_mm = 75  # Heatmap height in mm
legend_width_mm = 8  # Legend width in mm
legend_height_mm = 2  # Legend height in mm
font_type = "Arial"  # Font type
font_size = 6  # Font size

# Convert figure size from mm to inches (1 inch = 25.4 mm)
heatmap_width_in = heatmap_width_mm / 25.4
heatmap_height_in = heatmap_height_mm / 25.4
legend_width_in = legend_width_mm / 25.4
legend_height_in = legend_height_mm / 25.4

# Update font type and size globally
plt.rcParams.update({"font.family": font_type, "font.size": font_size})

# Create a heatmap dataframe with Geneid as the index
heatmap_data = overlap_df[["Geneid", "log2FoldChange_RP1", "log2FoldChange_RP2"]]
heatmap_data = heatmap_data.set_index("Geneid")

# Create the heatmap
plt.figure(figsize=(heatmap_width_in, heatmap_height_in), dpi=1200)
sns.heatmap(
    heatmap_data,
    annot=False,  # Remove annotations
    cmap="Blues_r",  # Blue color scale for downregulation
    cbar=False,  # Remove the color bar
    linewidths=0.1,
    yticklabels=False,
    xticklabels=False
)
plt.xlabel("Condition")
plt.ylabel("Genes (Ranked by Fold Change of R10 B+ vs. NH strain)")
plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_down", format="png")
# plt.show()

# Plot the legend as a separate figure
vmin, vmax = -5, 0  # Define the range for log2FoldChange

fig, ax = plt.subplots(figsize=(legend_width_in, legend_height_in), dpi=1200)  # Use input sizes
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = plt.cm.Blues_r  # Blue color scale for downregulation

# Add the colorbar to the figure
cb = mpl.colorbar.ColorbarBase(
    ax=ax,
    cmap=cmap,
    norm=norm,
    orientation="horizontal"
   
)

plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_down_colorbar", format="png")
# plt.show()


###################################
# Figure 4 c: dowregulated top 10 #
###################################

# Define parameters for figure size (in mm), font type, and font size
heatmap_width_mm = 15   # Heatmap width in mm
heatmap_height_mm = 37.5  # Heatmap height in mm
legend_width_mm = 10    # Legend width in mm
legend_height_mm = 5    # Legend height in mm
font_type = "Arial"     # Font type
font_size = 6           # Font size

# Convert figure size from mm to inches (1 inch = 25.4 mm)
heatmap_width_in = heatmap_width_mm / 25.4
heatmap_height_in = heatmap_height_mm / 25.4
legend_width_in = legend_width_mm / 25.4
legend_height_in = legend_height_mm / 25.4

# Update font type and size globally
plt.rcParams.update({"font.family": font_type, "font.size": font_size})

# Sort by log2FoldChange in RP1 (ascending so the most negative is at the top)
overlap_df = overlap_df.sort_values('log2FoldChange_RP1', ascending=True)

# Select the top 10 most negative genes
top_10_genes = overlap_df.nsmallest(10, "log2FoldChange_RP1")

# Create a heatmap dataframe with Geneid as the index
heatmap_data = top_10_genes[['Geneid', 'log2FoldChange_RP1', 'log2FoldChange_RP2']]
heatmap_data = heatmap_data.set_index('Geneid')

# Create the heatmap without annotations
plt.figure(figsize=(heatmap_width_in, heatmap_height_in), dpi=1200)
sns.heatmap(
    heatmap_data,
    annot=False,      # Remove annotations
    cmap='Blues_r',     # Blue color scale for downregulation
    cbar=False,       # Remove the color bar
    linewidths=0.25,
    yticklabels=False,
    xticklabels=False
)
plt.xlabel('Condition')
plt.ylabel('')
plt.tight_layout()
save_figure_panel("Figure_4", "panel_c_down_top10", format="png")
# plt.show()

# Plot the legend as a separate figure
vmin, vmax = -5, 0  # Define the range for log2FoldChange
fig, ax = plt.subplots(figsize=(legend_width_in, legend_height_in), dpi=1200)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = plt.cm.Blues

# Add the colorbar to the figure
cb = mpl.colorbar.ColorbarBase(
    ax=ax,
    cmap=cmap,
    norm=norm,
    orientation="vertical",
    label="Log2 Fold Change"
)

plt.tight_layout()
# save_figure_panel("Figure_4", "panel_c_down_legend", format="png")
# plt.show()

##############
# Figure 4 d #
##############

df1 = load_figure_data("Figure_4", "panel_d_part1")
df2 = load_figure_data("Figure_4", "panel_d_part2")

df1['abs_log2FoldChange'] = abs(df1['log2FoldChange'])
df2['abs_log2FoldChange'] = abs(df2['log2FoldChange'])

merged_df = pd.merge(df1, df2, on='Geneid', suffixes=('_RP_R10', '_RP_R20'))

significant_RP_R10 = (
    (merged_df['padj_RP_R10'] < 0.05) & 
    (merged_df['abs_log2FoldChange_RP_R10'] > 0.6)
)
significant_RP_R20 = (
    (merged_df['padj_RP_R20'] < 0.05) & 
    (merged_df['abs_log2FoldChange_RP_R20'] > 0.6)
)
significant_both = significant_RP_R10 & significant_RP_R20
non_significant = ~significant_RP_R10 & ~significant_RP_R20


def get_bubble_size(padj):
    if padj <= 0.00005:
        return 50
    elif padj <= 0.0005:
        return 40
    elif padj <= 0.005:
        return 25
    elif padj <= 0.05:
        return 18
    else:
        return 9


merged_df['bubble_size_RP_R10'] = merged_df['padj_RP_R10'].apply(get_bubble_size)
merged_df['bubble_size_RP_R20'] = merged_df['padj_RP_R20'].apply(get_bubble_size)
merged_df['bubble_size'] = merged_df[['bubble_size_RP_R10', 'bubble_size_RP_R20']].max(axis=1)

# ---------------------------------------------------------
# 2) DEFINE LEGEND ELEMENTS (unchanged)
# ---------------------------------------------------------
legend_elements_color = [
    Line2D([0], [0], marker='o', color='grey',  label='Non-significant', markersize=6),
    Line2D([0], [0], marker='o', color='#E69F00', label='Significant R10', markersize=6),
    Line2D([0], [0], marker='o', color='#009E73', label='Significant R20', markersize=6),
    Line2D([0], [0], marker='o', color='purple', label='Significant Both', markersize=6)
]

# For the bubble sizes, note that we take the sqrt of actual “area”
legend_elements_size = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(50), label='< 0.00005'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(40), label='< 0.0005'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(25), label='< 0.005'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(18), label='< 0.05'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(9), label='> 0.05')
]

# ---------------------------------------------------------
# 3) PLOT MAIN FIGURE WITHOUT ANY LEGENDS
# ---------------------------------------------------------
plt.figure(figsize=(3.2, 3.2), dpi=1200)
plt.rcParams.update({'font.size': 8})

# Plot non-significant
plt.scatter(
    merged_df[non_significant]['log2FoldChange_RP_R10'],
    merged_df[non_significant]['log2FoldChange_RP_R20'],
    s=merged_df[non_significant]['bubble_size'],
    alpha=0.4, color='grey',
    edgecolors='none',   # Remove circle outline
    linewidths=0      # No border
)

# Plot significant R10 only
plt.scatter(
    merged_df[significant_RP_R10 & ~significant_RP_R20]['log2FoldChange_RP_R10'],
    merged_df[significant_RP_R10 & ~significant_RP_R20]['log2FoldChange_RP_R20'],
    s=merged_df[significant_RP_R10 & ~significant_RP_R20]['bubble_size'],
    alpha=0.4, color='#E69F00',
    edgecolors='none',   # Remove circle outline
    linewidths=0      # No border
)

# Plot significant R20 only
plt.scatter(
    merged_df[significant_RP_R20 & ~significant_RP_R10]['log2FoldChange_RP_R10'],
    merged_df[significant_RP_R20 & ~significant_RP_R10]['log2FoldChange_RP_R20'],
    s=merged_df[significant_RP_R20 & ~significant_RP_R10]['bubble_size'],
    alpha=0.4, color='#009E73',
    edgecolors='none',   # Remove circle outline
    linewidths=0      # No border
)

# Plot significant both
plt.scatter(
    merged_df[significant_both]['log2FoldChange_RP_R10'],
    merged_df[significant_both]['log2FoldChange_RP_R20'],
    s=merged_df[significant_both]['bubble_size'],
    alpha=0.4, color='purple',
    edgecolors='none',   # Remove circle outline
    linewidths=0      # No border
)

plt.xlabel('Log2 FC R10 B+ vs NH strain')
plt.ylabel('Log2 FC R20 B+ vs NH strain')
save_figure_panel("Figure_4", "panel_d", format="png")
# plt.show()

# ---------------------------------------------------------
# 4) CREATE SEPARATE FIGURE FOR THE COLOR LEGEND
# ---------------------------------------------------------
fig_color_legend, ax_color_legend = plt.subplots(figsize=(1.8, 1.8), dpi=1200)
ax_color_legend.axis('off')  # no axes
color_legend = ax_color_legend.legend(
    handles=legend_elements_color,
    loc='center'
)
save_figure_panel("Figure_4", "panel_d_color_legend", format="png")
# plt.show()

# ---------------------------------------------------------
# 5) CREATE SEPARATE FIGURE FOR THE BUBBLE-SIZE LEGEND
# ---------------------------------------------------------
fig_size_legend, ax_size_legend = plt.subplots(figsize=(1.8, 1.8), dpi=1200)
ax_size_legend.axis('off')
size_legend = ax_size_legend.legend(
    handles=legend_elements_size,
    title='adjusted p-value',
    loc='center'
)
save_figure_panel("Figure_4", "panel_d_bubble_legend", format="png")
# plt.show()
