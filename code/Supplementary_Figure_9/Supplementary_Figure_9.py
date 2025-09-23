# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 08:41:59 2025

@author: tgassler
"""

# --- imports (kept similar) --------------------------------------------------
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.transforms as transforms
import os

# --- helpers -----------------------------------------------------------------
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

# ============================================================================ #
#                               CONFIG                                         #
# ============================================================================ #

METADATA_CSV   = "2025-05-19_24TG06_Rp_metadata.csv"
VSD_COUNTS_CSV = "2025-05-19_24TG06_vsd_Rp.csv.gz"

# your files look like:
# 2025-05-19_24TG06_Rp_<A>_vs_<B>_l0a0.01_results.csv
RESULTS_PREFIX = "2025-05-19_24TG06_Rp_"
RESULTS_SUFFIX = "_l0a0.01_results.csv.gz"

# figure style (identical to your version)
pca_figure_width_mm  = 80
pca_figure_height_mm = 50
font_type = "Arial"
font_size = 8
pca_dpi = 1200

figure_width_mm  = 45
figure_height_mm = 45
nDEG_width_mm    = 10
nDEG_heigth_mm   = 45
figure_resolution = 1200

# conversions
pca_fig_w_in  = pca_figure_width_mm / 25.4
pca_fig_h_in  = pca_figure_height_mm / 25.4
volcano_w_in  = figure_width_mm / 25.4
volcano_h_in  = figure_height_mm / 25.4
nDEG_w_in     = nDEG_width_mm / 25.4
nDEG_h_in     = nDEG_heigth_mm / 25.4

# fonts (kept)
plt.rcParams.update({"font.family": font_type, "font.size": font_size})

# thresholds (kept)
PADJ_THRESH = 0.05
LOG2FC_THR  = 0.6

# ============================================================================ #
#                               DATA LOAD                                      #
# ============================================================================ #
metadata = pd.read_csv(METADATA_CSV, index_col=0)
metadata.index.name = "bmk_id"
metadata.reset_index(inplace=True)
metadata["bmk_id"] = metadata["bmk_id"].astype(str).str.strip()

counts = pd.read_csv(VSD_COUNTS_CSV, index_col=0)
counts = counts.T.rename_axis("bmk_id").reset_index()

for df in (metadata, counts):
    df["bmk_id"] = df["bmk_id"].astype(str).str.strip()

merged_data = pd.merge(metadata, counts, on="bmk_id", how="inner")

# ============================================================================ #
#                         FILENAME BUILDER (FIXED)                             #
# ============================================================================ #
# observed convention from your folder:
#   results = RESULTS_PREFIX + A_token + "_vs_" + B_token + RESULTS_SUFFIX
#   where A_token = Condition with leading "Rp_" removed (if present)
#         B_token = Condition as-is (e.g., "Rp_R10_neg" or "pIND")
def strip_rp(tag):
    return tag[3:] if tag.startswith("Rp_") else tag

# explicit tokens matching your screenshot exactly
TOKEN_A = {
    'pIND':       'pIND',
    'Rp_R10_pos': 'R10_pos',
    'Rp_R10_neg': 'R10_neg',
    'Rp_R20_pos': 'R20_pos',
    'Rp_R20_neg': 'R20_neg',
}
TOKEN_B = {
    'pIND':       'pIND',
    'Rp_R10_pos': 'Rp_R10_pos',
    'Rp_R10_neg': 'Rp_R10_neg',
    'Rp_R20_pos': 'Rp_R20_pos',
    'Rp_R20_neg': 'Rp_R20_neg',
}

def build_results_path(A, B):
    candidates = [
        f"{RESULTS_PREFIX}{TOKEN_A.get(A, strip_rp(A))}_vs_{TOKEN_B.get(B, B)}{RESULTS_SUFFIX}",
        # fallbacks in case future files follow a different convention
        f"{RESULTS_PREFIX}{A}_vs_{B}{RESULTS_SUFFIX}",
        f"{RESULTS_PREFIX}{strip_rp(A)}_vs_{strip_rp(B)}{RESULTS_SUFFIX}",
        f"{RESULTS_PREFIX}{strip_rp(A)}_vs_{B}{RESULTS_SUFFIX}",
        f"{RESULTS_PREFIX}{A}_vs_{strip_rp(B)}{RESULTS_SUFFIX}",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    # default to the first naming we expect, but warn upstream
    return candidates[0]

# ============================================================================ #
#                               PCA PLOTTING                                   #
# ============================================================================ #
def make_pca(merged_df, conditions_to_include, condition_mapping):
    df = merged_df[merged_df['Condition'].isin(conditions_to_include)].copy()
    df['Condition'] = df['Condition'].map(condition_mapping)

    count_data = df.iloc[:, 8:]  # same assumption as your script

    pca = PCA(n_components=6)
    pca_result = pca.fit_transform(count_data)
    df['PCA1'] = pca_result[:, 0]
    df['PCA2'] = pca_result[:, 1]
    explained_variance = pca.explained_variance_ratio_

    label_order = [condition_mapping[c] for c in conditions_to_include]
    palette = sns.color_palette("colorblind", len(label_order))
    color_map = dict(zip(label_order, palette))

    fig = plt.figure(figsize=(pca_fig_w_in, pca_fig_h_in), dpi=pca_dpi)
    plt.rcParams.update({"font.family": font_type, "font.size": font_size})

    ax = sns.scatterplot(
        x='PCA1', y='PCA2',
        hue='Condition',
        hue_order=label_order,
        palette=color_map,
        data=df,
        legend="full"
    )

    for label in label_order:
        subset = df[df['Condition'] == label]
        if subset.shape[0] >= 2:
            confidence_ellipse(
                subset['PCA1'].values, subset['PCA2'].values, ax,
                n_std=1.5, facecolor=color_map[label],
                edgecolor=color_map[label], alpha=0.2
            )

    plt.xlabel(f'PCA1 ({explained_variance[0]:.1%})')
    plt.ylabel(f'PCA2 ({explained_variance[1]:.1%})')
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Condition", frameon=False)
    plt.tight_layout()
    #plt.show()
    return fig 

# ============================================================================ #
#                    VOLCANO + nDEG FOR A SINGLE COMPARISON                    #
# ============================================================================ #
def volcano_and_nDEG_for_comparison(results_csv_path):
    if not os.path.exists(results_csv_path):
        print(f"[WARN] Missing results file: {results_csv_path}")
        return

    df = pd.read_csv(results_csv_path)

    df['abs_log2FoldChange'] = df['log2FoldChange'].abs()
    safe_p = df['pvalue'].clip(lower=np.finfo(float).tiny)
    df['-log10(pvalue)'] = -np.log10(safe_p)

    # Volcano (identical aesthetics)
    volcano_fig = plt.figure(figsize=(volcano_w_in, volcano_h_in), dpi=figure_resolution)
    plt.scatter(df['log2FoldChange'], -np.log10(safe_p),
                edgecolor='grey', facecolors='none', alpha=0.5, s=7)
    significant_up = (df['padj'] < PADJ_THRESH) & (df['log2FoldChange'] >  LOG2FC_THR)
    plt.scatter(df.loc[significant_up, 'log2FoldChange'],
                -np.log10(safe_p[significant_up]),
                edgecolor='blue', facecolors='none', alpha=0.5, s=7)
    significant_down = (df['padj'] < PADJ_THRESH) & (df['log2FoldChange'] < -LOG2FC_THR)
    plt.scatter(df.loc[significant_down, 'log2FoldChange'],
                -np.log10(safe_p[significant_down]),
                edgecolor='orange', facecolors='none', alpha=0.5, s=7)
    plt.xlabel('Log2 FC')
    plt.ylabel('-Log10(p-value)')
    plt.tight_layout()
    #plt.show()

    # nDEG (identical look)
    deg_up   = (df['padj'] < PADJ_THRESH) & (df['log2FoldChange'] >  LOG2FC_THR)
    deg_down = (df['padj'] < PADJ_THRESH) & (df['log2FoldChange'] < -LOG2FC_THR)
    num_up, num_down = int(deg_up.sum()), int(deg_down.sum())

    ndeg_fig = plt.figure(figsize=(nDEG_w_in, nDEG_h_in), dpi=figure_resolution)
    bars_up = plt.bar(['DEGs'], [num_up],  color='blue',   alpha=0.5, label='DEGs Up')
    bars_dn = plt.bar(['DEGs'], [num_down], color='orange', alpha=0.5, bottom=[num_up], label='DEGs Down')

    total = num_up + num_down
    plt.text(0, total + 20, total, ha='center', va='bottom')
    # inside-bar labels
    for bar in bars_up:
        y = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, y/2, int(y), ha='center', va='bottom')
    for bar in bars_dn:
        y = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, num_up + y/2, int(y), ha='center', va='bottom')

    plt.ylabel('nDEG')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    #plt.show()
    return volcano_fig, ndeg_fig

# ============================================================================ #
#                              RUN THE SERIES                                  #
# ============================================================================ #

# --- R10 set -----------------------------------------------------------------
conditions_R10 = ['pIND', 'Rp_R10_pos', 'Rp_R10_neg']
mapping_R10 = {
    'pIND':       'NH strain',
    'Rp_R10_pos': 'R10 B$_{pos}$',
    'Rp_R10_neg': 'R10 B$_{neg}$'
}
pca_fig = make_pca(merged_data, conditions_R10, mapping_R10)
pca_fig.savefig("Supplementary_Figure_9_R10_pca.png", bbox_inches='tight', dpi=300)

comparisons_R10 = [
    ('Rp_R10_pos', 'pIND'),         # R10 Bpos vs. NH
    ('Rp_R10_pos', 'Rp_R10_neg'),   # R10 Bpos vs. R10 Bneg
    ('Rp_R10_neg', 'pIND')          # R10 Bneg vs. NH
]
for A, B in comparisons_R10:
    results_csv = build_results_path(A, B)
    vfig, ndegfig = volcano_and_nDEG_for_comparison(results_csv)
    vfig.savefig(f"Supplementary_Figure_9_R10_volcano_{A}_{B}.png", bbox_inches='tight', dpi=300)
    ndegfig.savefig(f"Supplementary_Figure_9_R10_ndeg_{A}_{B}.png", bbox_inches='tight', dpi=300)

# --- R20 set -----------------------------------------------------------------
conditions_R20 = ['pIND', 'Rp_R20_pos', 'Rp_R20_neg']
mapping_R20 = {
    'pIND':       'NH strain',
    'Rp_R20_pos': 'R20 B$_{pos}$',
    'Rp_R20_neg': 'R20 B$_{neg}$'
}
pca_fig = make_pca(merged_data, conditions_R20, mapping_R20)
pca_fig.savefig("Supplementary_Figure_9_R20_pca.png", bbox_inches='tight', dpi=300)

comparisons_R20 = [
    ('Rp_R20_pos', 'pIND'),         # R20 Bpos vs. NH
    ('Rp_R20_pos', 'Rp_R20_neg'),   # R20 Bpos vs. R20 Bneg
    ('Rp_R20_neg', 'pIND')          # R20 Bneg vs. NH
]
for A, B in comparisons_R20:
    results_csv = build_results_path(A, B)
    vfig, ndegfig = volcano_and_nDEG_for_comparison(results_csv)
    vfig.savefig(f"Supplementary_Figure_9_R20_volcano_{A}_{B}.png", bbox_inches='tight', dpi=300)
    ndegfig.savefig(f"Supplementary_Figure_9_R20_ndeg_{A}_{B}.png", bbox_inches='tight', dpi=300)
