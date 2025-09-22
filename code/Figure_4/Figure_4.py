# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 08:37:54 2025

@author: tgassler
"""
# Figure 4 aPCA plots:

   
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.transforms as transforms
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.colors as mcolors
from matplotlib_venn import venn3
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib as mpl
import gzip




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


# --- read files exactly as before -------------------------------------------
metadata = pd.read_csv("2025-05-19_24TG06_Rp_metadata.csv", index_col=0)
metadata.index.name = "bmk_id"
metadata.reset_index(inplace=True)
metadata["bmk_id"] = metadata["bmk_id"].astype(str).str.strip()

counts   = pd.read_csv("2025-05-19_24TG06_vsd_Rp.csv.gz", index_col=0)

# --- reshape counts: just transpose, no extra set_index ---------------------
counts = (
    counts.T                     # samples -> rows, genes -> columns
          .rename_axis("bmk_id") # give the index a name
          .reset_index()         # …and turn it into a normal column
)

# --- (optional but wise) strip stray whitespace -----------------------------
for df in (metadata, counts):
    df["bmk_id"] = df["bmk_id"].astype(str).str.strip()

# --- merge ------------------------------------------------------------------
merged_data = pd.merge(metadata, counts, on="bmk_id", how="inner")


# Transpose counts DataFrame to match the sample IDs in the metadata
counts = counts.set_index(counts.columns[0]).T.reset_index()
counts.columns = ['bmk_id'] + list(counts.columns[1:])


# Define the conditions to include
conditions_to_include = ['Rp_R10_pos', 'Rp_R20_pos', 'pIND'] #include as needed 'pIND': 'NH strain','Rp_R10_pos': 'R10 B+','Rp_R10_neg': 'R10 B-',#'Rp_R20_pos': 'R20 B+','Rp_R20_neg': 'R20 B-' 

# Filter the merged data to include only the desired conditions
filtered_data = merged_data[merged_data['Condition'].isin(conditions_to_include)]

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
#plt.show()
plt.savefig("Figure_4_a.png", bbox_inches='tight', dpi=300)


###########################################################

# Define inputs for figure size, font type, and font size
figure_width_mm = 45  # Width in mm
figure_height_mm = 45  # Height in mm
nDEG_width_mm = 10
nDEG_heigth_mm = 45
font_type = "Arial"
font_size = 8
figure_resolution=1200

# Convert figure size from mm to inches (1 inch = 25.4 mm)
figure_width_in = figure_width_mm / 25.4
figure_height_in = figure_height_mm / 25.4
nDEG_width_in = nDEG_width_mm /25.4
nDEG_heigth_in = nDEG_heigth_mm /25.4


# Set font properties globally for all plots
plt.rcParams.update({
    "font.family": font_type,
    "font.size": font_size
})

# Load the dataset
df = pd.read_csv("2025-05-19_24TG06_Rp_R20_pos_vs_pIND_l0a0.01_results.csv.gz")


plt.figure(figsize=(figure_width_in, figure_height_in), dpi=figure_resolution)
plt.hist(df['pvalue'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('p-value')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()


plt.figure(figsize=(figure_width_in, figure_height_in), dpi=figure_resolution)
plt.hist(df['padj'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('padj')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()


df['abs_log2FoldChange'] = abs(df['log2FoldChange'])
df['-log10(pvalue)'] = -np.log10(df['pvalue'])

# %% Create the volcano plot for different significance level for up pr down regulated ones

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=figure_resolution)

# All points
plt.scatter(df['log2FoldChange'], -np.log10(df['pvalue']), edgecolor='grey', facecolors='none', alpha=0.5, s=7)

# Significance criteria for upregulated genes
significant_up = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)
plt.scatter(df[significant_up]['log2FoldChange'], -np.log10(df[significant_up]['pvalue']), edgecolor='blue', facecolors='none', alpha=0.5, s=7)

# Significance criteria for downregulated genes
significant_down = (df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)
plt.scatter(df[significant_down]['log2FoldChange'], -np.log10(df[significant_down]['pvalue']), edgecolor='orange', facecolors='none', alpha=0.5, s=7)

# Labels and Axes limits
plt.xlabel('Log2 FC')
plt.ylabel('-Log10(p-value)')
#plt.show()
plt.savefig("Sup_Figure_9_f_p1.png", bbox_inches='tight', dpi=300)

# Filter for DEGs with fold change > 1 and padj < 0.01
deg_criteria = (df['padj'] < 0.05) & (df['abs_log2FoldChange'] > 0.6)
deg_df = df[deg_criteria]
#upregulated ones
deg_criteria = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)
deg_up_df = df[deg_criteria]

#downregulated ones
deg_criteria = (df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)
deg_down_df = df[deg_criteria]


# Filter DataFrame to include only rows where 'Description' contains 'retro'
retro_deg_df = deg_df[deg_df['Description'].str.contains('retro', case=False, na=False)]




# %% Plotting the nDEGs

# Calculate the number of DEGs for each category
num_deg_up = len(deg_up_df)
num_deg_down = len(deg_down_df)

# Data for plotting
categories = ['DEGs']
counts_up = [num_deg_up]
counts_down = [num_deg_down]



# Create the bar plot as a stacked bar plot
plt.figure(figsize=(nDEG_width_in , nDEG_heigth_in), dpi=figure_resolution)
bars_up = plt.bar(categories, counts_up, color='blue', alpha=0.5, label='DEGs Up')
bars_down = plt.bar(categories, counts_down, color='orange', alpha=0.5, bottom=counts_up, label='DEGs Down')

# Adding the counts above the bars for the total
total_counts = [x + y for x, y in zip(counts_up, counts_down)]
for idx, val in enumerate(total_counts):
    plt.text(idx, val + 75, val, ha='center', va='bottom')

# Adding the counts above the individual bars
for bar in bars_up:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval/2, yval, ha='center', va='bottom')

for bar in bars_down:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + counts_up[bars_down.index(bar)] - yval/2, yval, ha='center', va='bottom')

# Adding labels and title
plt.ylabel('nDEG')
#plt.title('Total Number of DEGs with Upregulated and Downregulated Counts')
# Position the legend to the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show the plot
#plt.show()
plt.savefig("Sup_Figure_9_f_p2.png", bbox_inches='tight', dpi=300)

#%%
# Extract a list of Gene IDs for DEGs
deg_list = deg_df[['Geneid', 'Name', 'KEGG_Pathway', "Description", 'log2FoldChange']].values.tolist()

# Extract a list of upregulated Gene IDs for DEGs
deg_up_list = deg_up_df[['Geneid', 'Name', 'KEGG_Pathway', "Description", 'log2FoldChange']].values.tolist()

# Extract a list of downregulated Gene IDs for DEGs
deg_down_list = deg_down_df[['Geneid', 'Name', 'KEGG_Pathway', "Description", 'log2FoldChange']].values.tolist()

# Print the count and proportion of DEGs
number_of_degs = len(deg_list)
number_of_degs_up= len(deg_up_list)
number_of_degs_down= len(deg_down_list)
total_genes = df.shape[0]
proportion_of_degs = number_of_degs / total_genes
proportion_up = number_of_degs_up / total_genes
proportion_down = number_of_degs_down / total_genes

print(f"Number of DEGs with fold change > 1.5: {number_of_degs}")
print(f"Total number of genes analyzed: {total_genes}")
print(f"Proportion of DEGs with fold change > 1.5: {proportion_of_degs:.2%}")
print(f"Proportion of DEGs upregulated: {proportion_up:.2%}")
print(f"Proportion of DEGs downredulated: {proportion_down:.2%}")


#%% Making the CSV files

#deg_df.to_csv("nDEG_rp_R10.csv", index=False)


#############################
# Figure 4 b: Venn Diagramm

# Load your dataframes
deg_RP1 = pd.read_csv("nDEG_rp_R10.csv")
deg_RP2 = pd.read_csv("nDEG_rp_R20.csv")

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
color_map = {
    'pIND': "#0072B2",       # Green
    'Rp_R10_pos': "#E69F00", # Blue
    'Rp_R20_pos': "#009E73"  # Orange
}


set_colors = [color_map['Rp_R10_pos'], color_map['Rp_R20_pos'], color_map['pIND']]

# Create a figure with GridSpec layout
fig = plt.figure(figsize=(3.5, 2.1), dpi=1200)
plt.rcParams.update({'font.size': 7})
gs = GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1])

# Venn diagram for all nDEGs
ax0 = fig.add_subplot(gs[:, 0])
venn3([set_RP1, set_RP2, set()], ('', '', ''), ax=ax0, set_colors=set_colors)
#ax0.set_title("All nDEGs")

# Venn diagram for upregulated genes
ax1 = fig.add_subplot(gs[0, 1])
venn3([upregulated_RP1, upregulated_RP2, set()], ('', '',  ''), ax=ax1, set_colors=set_colors)
#venn3([upregulated_RP1, upregulated_RP2, set()], ('R10 B+', 'R20 B+',  ''), ax=ax1, set_colors=set_colors)
#ax1.set_title("Upregulated")

# Venn diagram for downregulated genes
ax2 = fig.add_subplot(gs[1, 1])
venn3([downregulated_RP1, downregulated_RP2, set()], ('', '', ''), ax=ax2, set_colors=set_colors)
#venn3([downregulated_RP1, downregulated_RP2, set()], ('R10 B+', 'R20 B+', ''), ax=ax2, set_colors=set_colors)
#ax2.set_title("Downregulated")

# Adjust layout
plt.tight_layout()
#plt.show()
plt.savefig("Figure_4_b.png", bbox_inches='tight', dpi=300)


# #%% #Figure 4 c


# Identify overlapping upregulated genes
overlapping_genes = upregulated_RP1.intersection(upregulated_RP2)

# # Filter rows for the overlapping genes in both dataframes
overlap_RP1 = deg_RP1[deg_RP1['Geneid'].isin(overlapping_genes)]
overlap_RP2 = deg_RP2[deg_RP2['Geneid'].isin(overlapping_genes)]

# Merge dataframes to include details from both deg_RP1 and deg_RP2
overlap_df = pd.merge(
    overlap_RP1,
    overlap_RP2,
    on='Geneid',
    suffixes=('_RP1', '_RP2')
)

# # Save the result to a CSV file for further analysis (optional)
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
#plt.title("Heatmap of Overlapping Upregulated Genes (Ranked by Fold Change in RP1)")
plt.xlabel("Condition")
plt.ylabel("Genes (Ranked by Fold Change of R10 B$_{pos}$ vs. NH strain)")
plt.tight_layout()
#plt.show()
plt.savefig("Figure_4_c_up.png", bbox_inches='tight', dpi=300)


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
plt.savefig("Figure_4_c_up_colorbar.png", bbox_inches='tight', dpi=300)
#plt.show()


# #%% Figure 4c Top10 genes upregulated


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

# Select the top 25% based on log2FoldChange_RP1
#top_25_percent = overlap_df.head(len(overlap_df) // 10)

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
#plt.title('Heatmap of Top 25% Overlapping Upregulated Genes (Ranked by Fold Change in RP1)')
plt.xlabel('Condition')
plt.ylabel('')
plt.tight_layout()
#plt.show()
plt.savefig("Figure_4_c_up_top10.png", bbox_inches='tight', dpi=300)
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
# plt.show()

# #%% Fogure 4c : Downregulated Genes


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

# Identify overlapping downregulated genes
overlapping_genes = downregulated_RP1.intersection(downregulated_RP2)

# Filter rows for the overlapping genes in both dataframes
overlap_RP1 = deg_RP1[deg_RP1["Geneid"].isin(overlapping_genes)]
overlap_RP2 = deg_RP2[deg_RP2["Geneid"].isin(overlapping_genes)]

# Merge dataframes
overlap_df = pd.merge(
    overlap_RP1,
    overlap_RP2,
    on="Geneid",
    suffixes=("_RP1", "_RP2")
)

# Cap extreme values for better scaling in the heatmap
overlap_df["log2FoldChange_RP1"] = overlap_df["log2FoldChange_RP1"].clip(-5, 0)
overlap_df["log2FoldChange_RP2"] = overlap_df["log2FoldChange_RP2"].clip(-5, 0)

# Sort by log2FoldChange in RP1
overlap_df = overlap_df.sort_values("log2FoldChange_RP1")

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
plt.ylabel("Genes (Ranked by Fold Change of R10 B$_{pos}$ vs. NH strain)")
plt.tight_layout()
plt.savefig("Figure_4_c_down.png", bbox_inches='tight', dpi=300)
#plt.show()

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
plt.savefig("Figure_4_c_down_colorbar.png", bbox_inches='tight', dpi=300)
# plt.show()

# #%% Figure 4 c: Top10 dowregulated genes

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
#plt.show()
plt.savefig("Figure_4_c_down_top10.png", bbox_inches='tight', dpi=300)


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
# plt.show()

# #%% This sod can be used to retrieve the AA sequences from the annotation file


# Depending on which df you are looking at, defines which ones will be extracted

df = overlap_df.nlargest(10, "log2FoldChange_RP1")
#df = overlap_df.nsmallest(10, "log2FoldChange_RP1")
#df = overlap_df #consider which one you are working with


fasta_path = "braker.faa.gz"
output_fasta_path = "Top10_Down.fasta" # "Overlap_Up.fasta"; "Top10_Down.fasta";  "Overlap_Up.fasta"; "Overlap_Down.fasta"

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

#%% Figure 4d

# ---------------------------------------------------------
# ONE-LINE KNOB – how much larger the bubbles appear inside
# the *size* legend only (1.0 = exact, >1 = larger)
# ---------------------------------------------------------
LEGEND_SCALE = 0.9
# ---------------------------------------------------------

# 0) HELPERS ------------------------------------------------
def get_bubble_area(padj):
    if padj <= 5e-5:
        return 50
    elif padj <= 5e-4:
        return 40
    elif padj <= 5e-3:
        return 25
    elif padj <= 5e-2:
        return 18
    else:
        return 9

def area_to_diameter(area):
    return 2 * np.sqrt(area / np.pi)          # points → points

# 1) READ + PREP DATA --------------------------------------
df1 = pd.read_csv("2025-05-19_24TG06_Rp_R10_pos_vs_pIND_l0a0.01_results.csv.gz")
df2 = pd.read_csv("2025-05-19_24TG06_Rp_R20_pos_vs_pIND_l0a0.01_results.csv.gz")

for df in (df1, df2):
    df["abs_log2FoldChange"] = df["log2FoldChange"].abs()

merged = pd.merge(df1, df2, on="Geneid", suffixes=("_R10", "_R20"))

sig_R10  = (merged["padj_R10"] < 0.05) & (merged["abs_log2FoldChange_R10"] > 0.6)
sig_R20  = (merged["padj_R20"] < 0.05) & (merged["abs_log2FoldChange_R20"] > 0.6)
sig_both = sig_R10 & sig_R20
nonsig   = ~sig_R10 & ~sig_R20

for side in ("R10", "R20"):
    merged[f"bubble_{side}"] = merged[f"padj_{side}"].apply(get_bubble_area)
merged["bubble_area"] = merged[["bubble_R10", "bubble_R20"]].max(axis=1)

# 2) LEGEND HANDLES ----------------------------------------
legend_color = [
    Line2D([0], [0], marker="o", color="grey", markeredgecolor="none", alpha=0.8, linestyle="None", label="Non-significant", markersize=8),
    Line2D([0], [0], marker="o", color="#E69F00", alpha=0.8, markeredgecolor="none", linestyle="None", label="Significant R10", markersize=8),
    Line2D([0], [0], marker="o", color="#009E73", alpha=0.8, markeredgecolor="none", linestyle="None", label="Significant R20", markersize=8),
    Line2D([0], [0], marker="o", color="purple", alpha=0.8, markeredgecolor="none", linestyle="None",  label="Significant both", markersize=8),
]

sizes  = [50, 40, 25, 18, 9]
labels = ["≤ 0.00005", "≤ 0.0005", "≤ 0.005", "≤ 0.05", "> 0.05"]
legend_size = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",alpha=0.4,markeredgecolor="none",
           markersize=area_to_diameter(a), label=lab)
    for a, lab in zip(sizes, labels)
]

# 3) MAIN SCATTER (no legends) -----------------------------
plt.rcParams.update({"font.size": 8,"font.family": "Arial" })

plt.figure(figsize=(3.2, 3.2), dpi=1200)

plt.scatter(merged.loc[nonsig,   "log2FoldChange_R10"],
            merged.loc[nonsig,   "log2FoldChange_R20"],
            s=merged.loc[nonsig, "bubble_area"],   color="grey",   alpha=0.4, edgecolors="none")

plt.scatter(merged.loc[sig_R10 & ~sig_R20, "log2FoldChange_R10"],
            merged.loc[sig_R10 & ~sig_R20, "log2FoldChange_R20"],
            s=merged.loc[sig_R10 & ~sig_R20, "bubble_area"], color="#E69F00", alpha=0.4, edgecolors="none")

plt.scatter(merged.loc[sig_R20 & ~sig_R10, "log2FoldChange_R10"],
            merged.loc[sig_R20 & ~sig_R10, "log2FoldChange_R20"],
            s=merged.loc[sig_R20 & ~sig_R10, "bubble_area"], color="#009E73", alpha=0.4, edgecolors="none")

plt.scatter(merged.loc[sig_both, "log2FoldChange_R10"],
            merged.loc[sig_both, "log2FoldChange_R20"],
            s=merged.loc[sig_both, "bubble_area"], color="purple", alpha=0.4, edgecolors="none")

plt.xlabel(r"Log$_2$ FC R10 $B_{pos}$ vs. NH strain")
plt.ylabel(r"Log$_2$ FC R20 $B_{pos}$ vs. NH strain")
plt.tight_layout()
plt.savefig("Figure_4_d.png", bbox_inches='tight', dpi=300)
#plt.show()

# 4) COLOR LEGEND -----------------------------------------
fig_c, ax_c = plt.subplots(figsize=(1.8, 1.8), dpi=1200)
ax_c.axis("off")
ax_c.legend(handles=legend_color, loc="center", frameon=False)
#plt.show()
plt.savefig("Figure_4_d_color_legend.png", bbox_inches='tight', dpi=300)

# 5) SIZE LEGEND ------------------------------------------
fig_s, ax_s = plt.subplots(figsize=(1.8, 1.8), dpi=1200)
ax_s.axis("off")
ax_s.legend(handles=legend_size,
            title="adjusted p-value",
            loc="center",
            frameon=False,
            markerscale=LEGEND_SCALE,   # ← tweak this one number
            labelspacing=1.2,
            handletextpad=0.8,
            borderpad=0.4)
#plt.show()
plt.savefig("Figure_4_d_size_legend.png", bbox_inches='tight', dpi=300)


#%% Exztracting the significant ones in both sets


# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord


# Please replace '/path/to/your/IND_B3vsB1.xlsx' with the actual path to your file
df1 = pd.read_csv("2025-05-19_24TG06_Rp_R10_pos_vs_pIND_l0a0.01_results.csv.gz")
df2 = pd.read_csv("2025-05-19_24TG06_Rp_R20_pos_vs_pIND_l0a0.01_results.csv.gz")

# Calculate absolute log2 fold change for significance criteria
df1['abs_log2FoldChange'] = abs(df1['log2FoldChange'])
df2['abs_log2FoldChange'] = abs(df2['log2FoldChange'])

# Merge the dataframes on a common column, assuming 'Geneid' is the common column
merged_df = pd.merge(df1, df2, on='Geneid', suffixes=('_RP_R10', '_RP_R20'))

# Determine significance based on given criteria
significant_RP_R10 = (merged_df['padj_RP_R10'] < 0.05) & (merged_df['abs_log2FoldChange_RP_R10'] > 0.6)
significant_RP_R20 = (merged_df['padj_RP_R20'] < 0.05) & (merged_df['abs_log2FoldChange_RP_R20'] > 0.6)
significant_both = significant_RP_R10 & significant_RP_R20
significant_RP_R10_only = significant_RP_R10 & ~significant_RP_R20
significant_RP_R20_only = significant_RP_R20 & ~significant_RP_R10

# Filter significant genes from both df1 and df2
significant_genes_both = merged_df[significant_both]['Geneid']
significant_genes_RP_R10_only = merged_df[significant_RP_R10_only]['Geneid']
significant_genes_RP_R20_only = merged_df[significant_RP_R20_only]['Geneid']

# Filter dataframes for each category
df1_filtered_both = df1[df1['Geneid'].isin(significant_genes_both)]
df2_filtered_both = df2[df2['Geneid'].isin(significant_genes_both)]
df1_filtered_RP_R10_only = df1[df1['Geneid'].isin(significant_genes_RP_R10_only)]
df2_filtered_RP_R20_only = df2[df2['Geneid'].isin(significant_genes_RP_R20_only)]

# Merge filtered dataframes on 'Geneid' and include relevant columns for each category
merged_filtered_df_both = pd.merge(
    df1_filtered_both[['Geneid', 'Name', 'KEGG_Pathway', 'Description', 'log2FoldChange']],
    df2_filtered_both[['Geneid', 'Name', 'KEGG_Pathway', 'Description', 'log2FoldChange']],
    on='Geneid', suffixes=('_RP_R10', '_RP_R20')
)

merged_filtered_df_RP_R10_only = df1_filtered_RP_R10_only[['Geneid', 'Name', 'KEGG_Pathway', 'Description', 'log2FoldChange']]
merged_filtered_df_RP_R20_only = df2_filtered_RP_R20_only[['Geneid', 'Name', 'KEGG_Pathway', 'Description', 'log2FoldChange']]

# Separate into upregulated and downregulated lists for both significant category
deg_up_df_both = merged_filtered_df_both[merged_filtered_df_both['log2FoldChange_RP_R10'] > 0]  # Upregulated genes
deg_down_df_both = merged_filtered_df_both[merged_filtered_df_both['log2FoldChange_RP_R10'] < 0]  # Downregulated genes

# Separate into upregulated and downregulated lists for MR_R10 only
deg_up_df_RP_R10_only = merged_filtered_df_RP_R10_only[merged_filtered_df_RP_R10_only['log2FoldChange'] > 0]  # Upregulated genes
deg_down_df_RP_R10_only = merged_filtered_df_RP_R10_only[merged_filtered_df_RP_R10_only['log2FoldChange'] < 0]  # Downregulated genes

# Separate into upregulated and downregulated lists for RP_R10 only
deg_up_df_RP_R20_only = merged_filtered_df_RP_R20_only[merged_filtered_df_RP_R20_only['log2FoldChange'] > 0]  # Upregulated genes
deg_down_df_RP_R20_only = merged_filtered_df_RP_R20_only[merged_filtered_df_RP_R20_only['log2FoldChange'] < 0]  # Downregulated genes

# Save upregulated and downregulated DEGs to CSV files for both significant category
#deg_up_df_both.to_csv('RP_R10_Pos_VS_RP_R20_Pos_IND_up.csv', index=False)
#deg_down_df_both.to_csv('RP_R10_Pos_VS_RP_R20_Pos_IND_down.csv', index=False)

# Save upregulated and downregulated DEGs to CSV files for MR_R10 only category
#deg_up_df_RP_R10_only.to_csv('RP_R10_only_up.csv', index=False)
#deg_down_df_RP_R10_only.to_csv('RP_R10_only_down.csv', index=False)

# Save upregulated and downregulated DEGs to CSV files for RP_R10 only category
#deg_up_df_RP_R20_only.to_csv('RP_R20_only_up.csv', index=False)
#deg_down_df_RP_R20_only.to_csv('RP_R20_only_down.csv', index=False)

#%% Get the target list for GO term enrichment


# Define the paths to your files -> take files generated above
csv_path= 'RP_R20_only_down.csv'
fasta_path = "braker.faa.gz"
output_fasta_path = "RP_R20_only_down.fasta"

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

#Submitt it to StringDb against our proteome STRG0A18OCV
 

#%% including the entires dfs
# Define a function to count occurrences of pooled keywords and track GeneIDs
def count_pooled_keywords_with_geneids(description_series, geneid_series, keywords):
    keyword_count = 0
    geneid_set = set()
    total_descriptions = 0
    for description, geneid in zip(description_series, geneid_series):
        if pd.notnull(description):
            total_descriptions += 1
            if any(keyword.lower() in description.lower() for keyword in keywords):
                keyword_count += 1
                geneid_set.add(geneid)
    return keyword_count, total_descriptions, geneid_set

# List of keywords to search for (pooled together)
keywords = ['retro', 'transposase', 'reverse']

# Apply the function to each category
count_both, total_both, geneids_both = count_pooled_keywords_with_geneids(merged_filtered_df_both['Description_RP_R10'], merged_filtered_df_both['Geneid'], keywords)
count_RP_R10_only, total_RP_R10_only, geneids_RP_R10_only = count_pooled_keywords_with_geneids(merged_filtered_df_RP_R10_only['Description'], merged_filtered_df_RP_R10_only['Geneid'], keywords)
count_RP_R20_only, total_RP_R20_only, geneids_RP_R20_only = count_pooled_keywords_with_geneids(merged_filtered_df_RP_R20_only['Description'], merged_filtered_df_RP_R20_only['Geneid'], keywords)

# Calculate keyword occurrences in the entire datasets
count_df1, total_df1, geneids_df1 = count_pooled_keywords_with_geneids(df1['Description'], df1['Geneid'], keywords)
count_df2, total_df2, geneids_df2 = count_pooled_keywords_with_geneids(df2['Description'], df2['Geneid'], keywords)

# Print the results with GeneIDs
def print_pooled_keyword_results(count, total, geneids, category_name):
    print(f"\nFrequency of pooled keywords in {category_name}:")
    geneid_list = list(geneids)
    print(f"Pooled keywords: {count} out of {total} ({count / total:.2%})")
    print(f"GeneIDs: {', '.join(geneid_list)}")

print_pooled_keyword_results(count_both, total_both, geneids_both, "both significant category")
print_pooled_keyword_results(count_RP_R10_only, total_RP_R10_only, geneids_RP_R10_only, "RP_R10 only category")
print_pooled_keyword_results(count_RP_R20_only, total_RP_R20_only, geneids_RP_R20_only, "RP_R20 only category")

print_pooled_keyword_results(count_df1, total_df1, geneids_df1, "entire df1")
print_pooled_keyword_results(count_df2, total_df2, geneids_df2, "entire df2")

# Determine which GeneIDs occur in RP_R10 only, RP_R20 only, or in both
def compare_geneids_between_categories(geneids_RP_R10, geneids_RP_R20, category_name):
    common_geneids = geneids_RP_R10 & geneids_RP_R20
    unique_RP_R10_geneids = geneids_RP_R10 - geneids_RP_R20
    unique_RP_R20_geneids = geneids_RP_R20 - geneids_RP_R10

    print(f"\nComparison of GeneIDs in {category_name}:")
    print(f"Common GeneIDs: {', '.join(common_geneids) if common_geneids else 'None'}")
    print(f"Unique GeneIDs for RP_R10: {', '.join(unique_RP_R10_geneids) if unique_RP_R10_geneids else 'None'}")
    print(f"Unique GeneIDs for RP_R20: {', '.join(unique_RP_R20_geneids) if unique_RP_R20_geneids else 'None'}")

compare_geneids_between_categories(geneids_RP_R10_only, geneids_RP_R20_only, "RP_R10 only vs RP_R20 only category")

#%%
import pandas as pd
from statsmodels.stats.proportion import proportions_ztest

# Observed counts and sample sizes
count_both = 4
total_both = 277
count_RP_R10_only = 25
total_RP_R10_only = 1500
count_RP_R20_only = 20
total_RP_R20_only = 591
count_df1 = 401
total_df1 = 21522
count_df2 = 401
total_df2 = 21522

# Pooled total counts and sample sizes for the entire dataset
count_entire = count_df1 + count_df2
total_entire = total_df1 + total_df2

# Proportion test for Both Significant Category vs. Entire Dataset
stat_both, pval_both = proportions_ztest([count_both, count_entire], [total_both, total_entire])

# Proportion test for RP_R10 Only Category vs. Entire Dataset
stat_RP_R10_only, pval_RP_R10_only = proportions_ztest([count_RP_R10_only, count_entire], [total_RP_R10_only, total_entire])

# Proportion test for RP_R20 Only Category vs. Entire Dataset
stat_RP_R20_only, pval_RP_R20_only = proportions_ztest([count_RP_R20_only, count_entire], [total_RP_R20_only, total_entire])

# Print the results
print(f"Proportion test for Both Significant Category vs. Entire Dataset: p-value = {pval_both:.4f}")
print(f"Proportion test for RP_R10 Only Category vs. Entire Dataset: p-value = {pval_RP_R10_only:.4f}")
print(f"Proportion test for RP_R20 Only Category vs. Entire Dataset: p-value = {pval_RP_R20_only:.4f}")

