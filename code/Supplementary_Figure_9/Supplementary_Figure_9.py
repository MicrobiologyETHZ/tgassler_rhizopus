# -*- coding: utf-8 -*-
"""
Created on Wed May 21 16:30:27 2025

@author: tgassler
"""

# PCA plots:
  
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.transforms as transforms

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
conditions_to_include = ['Rp_R20_pos', 'Rp_R20_neg', 'pIND'] #include as needed 'pIND': 'NH strain','Rp_R10_pos': 'R10 B+','Rp_R10_neg': 'R10 B-',#'Rp_R20_pos': 'R20 B+','Rp_R20_neg': 'R20 B-' 

# Filter the merged data to include only the desired conditions
filtered_data = merged_data[merged_data['Condition'].isin(conditions_to_include)]

# Map old condition names to new ones
condition_mapping = {
    'pIND': 'NH strain',
    #'Rp_R10_pos': 'R10 B$_{pos}$',
    #'Rp_R10_neg': 'R10 B$_{neg}$',
    'Rp_R20_pos': 'R20 B$_{pos}$',
    'Rp_R20_neg': 'R20 B$_{neg}$'
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

# Plot PCA with centroids and ellipses

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
plt.savefig("Supplementary_Figure_9_p1.png", bbox_inches='tight', dpi=300)
#%% Making the single comparisons and the lists
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.colors as mcolors

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
# Please replace '/path/to/your/IND_B3vsB1.xlsx' with the actual path to your file
df = pd.read_csv("2025-05-19_24TG06_Rp_R20_pos_vs_Rp_R10_pos_l0a0.01_results.csv.gz")

# %% Choosing of FDR

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=figure_resolution)
plt.hist(df['pvalue'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('p-value')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()
plt.savefig("Supplementary_Figure_9_p2.png", bbox_inches='tight', dpi=300)

plt.figure(figsize=(figure_width_in, figure_height_in), dpi=figure_resolution)
plt.hist(df['padj'], bins=100, color='skyblue', edgecolor='black')
plt.axvline(x=0.05, color='red', linestyle='--')  # Red dashed line at p-value = 0.05
plt.ylim(0, 6000)
plt.xlabel('padj')
plt.ylabel('Number of Genes')
plt.grid(True)
#plt.show()
plt.savefig("Supplementary_Figure_9_p3.png", bbox_inches='tight', dpi=300)
#%% # Calculate absolute log2 fold change and -log10 of p-value

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
plt.savefig("Supplementary_Figure_9_p4.png", bbox_inches='tight', dpi=300)

# %%

# Filter for DEGs with fold change > 1 and padj < 0.01
deg_criteria = (df['padj'] < 0.05) & (df['abs_log2FoldChange'] > 0.6)
deg_df = df[deg_criteria]
#upregulated ones
deg_criteria = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0.6)
deg_up_df = df[deg_criteria]

#downregulated ones
deg_criteria = (df['padj'] < 0.05) & (df['log2FoldChange'] < -0.6)
deg_down_df = df[deg_criteria]

#making the csv files

#deg_df.to_csv("R10pos_vs_IND.csv") (change as needed)



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
    plt.text(idx, val + 20, val, ha='center', va='bottom')

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
plt.savefig("Supplementary_Figure_9_p5.png", bbox_inches='tight', dpi=300)

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



