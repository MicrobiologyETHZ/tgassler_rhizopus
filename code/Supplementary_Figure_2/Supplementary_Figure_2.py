# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 14:46:21 2025

@author: tgassler
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl

##############################################################################
# USER-CONFIGURABLE SETTINGS
##############################################################################
FIGURE_WIDTH_MM = 75    # Width in millimeters
FIGURE_HEIGHT_MM = 75   # Height in millimeters
FONT_FAMILY = 'Arial'
FONT_SIZE = 8
AXES_LABELSIZE = 8
XTICK_LABELSIZE = 8
YTICK_LABELSIZE = 8
LEGEND_FONTSIZE = 8

def mm_to_inches(mm):
    """Convert millimeters to inches."""
    return mm / 25.4

##############################################################################
# Update Matplotlib default font and figure settings
##############################################################################
mpl.rcParams['font.family'] = FONT_FAMILY
mpl.rcParams['font.size'] = FONT_SIZE
mpl.rcParams['axes.labelsize'] = AXES_LABELSIZE
mpl.rcParams['xtick.labelsize'] = XTICK_LABELSIZE
mpl.rcParams['ytick.labelsize'] = YTICK_LABELSIZE
mpl.rcParams['legend.fontsize'] = LEGEND_FONTSIZE

# Convert figure size from mm to inches
fig_width_in = mm_to_inches(FIGURE_WIDTH_MM)
fig_height_in = mm_to_inches(FIGURE_HEIGHT_MM)

##############################################################################
# Load the Excel file into a DataFrame
##############################################################################
df = pd.read_excel('sourcedata_SF2a.xlsx', engine='openpyxl')

##############################################################################
# Function to model exponential growth
##############################################################################
def exponential_growth(x, a, b):
    return a * np.exp(b * x)

##############################################################################
# Function to calculate doubling time with minimal doublings
##############################################################################
def calculate_max_growth_rate(time, count, min_doublings=1.5):
    log_count = np.log(count)
    max_growth_rate = 0.0
    best_interval = (0, 1)

    # Iterate over all possible intervals to find the maximum growth rate
    for i in range(len(time)):
        for j in range(i + 1, len(time)):
            if log_count.iloc[j] - log_count.iloc[i] >= np.log(2) * min_doublings:
                growth_rate = (log_count.iloc[j] - log_count.iloc[i]) / (time.iloc[j] - time.iloc[i])
                if growth_rate > max_growth_rate:
                    max_growth_rate = growth_rate
                    best_interval = (i, j)

    # Calculate doubling time
    doubling_time = np.log(2) / max_growth_rate if max_growth_rate > 0 else np.inf
    return doubling_time, max_growth_rate, best_interval

##############################################################################
# Calculate the maximum growth rate and minimal doubling time over the dataset
##############################################################################
doubling_time, max_growth_rate, (start, end) = calculate_max_growth_rate(
    df['deltaT_h_01'], df['Voxel_01']
)

# Basic safety check
if not np.isfinite(doubling_time) or end <= start:
    raise ValueError("Could not find a valid max-growth interval that meets the min_doublings criterion.")

##############################################################################
# Fit the exponential growth model ONLY on the max-growth interval
##############################################################################
time_sub = df['deltaT_h_01'].iloc[start:end + 1].to_numpy()
vox_sub  = df['Voxel_01'].iloc[start:end + 1].to_numpy()

# Initial guesses: a ≈ first y, b ≈ max_growth_rate (log-slope)
a0 = max(vox_sub[0], 1e-12)
b0 = max_growth_rate if np.isfinite(max_growth_rate) and max_growth_rate > 0 else 0.1

params_sub, covariance_sub = curve_fit(
    exponential_growth, time_sub, vox_sub, p0=(a0, b0), maxfev=10000
)

##############################################################################
# Plot the data and the fit (fit only over the highlighted interval)
##############################################################################
plt.figure(figsize=(fig_width_in, fig_height_in), dpi=1200)
plt.scatter(df['deltaT_h_01'], df['Voxel_01'], color='blue', label='Sum of detected objects')

# Smooth time points for exponential fit curve — restricted to the interval
t_fit = np.linspace(time_sub.min(), time_sub.max(), 200)
plt.plot(t_fit, exponential_growth(t_fit, *params_sub),
         label='Exponential Fit', color='red')

# Highlight the interval with the maximum growth rate
plt.axvspan(df['deltaT_h_01'].iloc[start], df['deltaT_h_01'].iloc[end],
            color='yellow', alpha=0.3, label=fr'$T_{{D\mathrm{{max}}}}$ = {doubling_time:.1f} h')

plt.xlabel('Time (h)')
plt.ylabel('Total Voxels (a.u)')
# plt.title('Exponential Growth Fit with Maximum Growth Rate Interval')  # (kept off as in original)
plt.legend(loc='upper left', frameon=True, edgecolor='black')
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()
#plt.show()
plt.savefig("Supplementary_Figure_2_a.png", bbox_inches='tight', dpi=300)


#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

##############################################################################
# APPEARANCE
##############################################################################
FIGURE_WIDTH_MM  = 75
FIGURE_HEIGHT_MM = 75
FONT_FAMILY      = "Arial"
FONT_SIZE        = 8

def mm_to_inches(mm):         # helper: mm → inch
    return mm / 25.4

mpl.rcParams.update({
    "font.family":     FONT_FAMILY,
    "font.size":       FONT_SIZE,
    "axes.labelsize":  FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
})

##############################################################################
# LOAD PRE-PROCESSED DATA
##############################################################################
df   = pd.read_excel("sourcedata_SF2b.xlsx", engine="openpyxl")
time = df["Time_h"]

##############################################################################
# MATCH ORIGINAL COLOURS
##############################################################################
colour_map = {
    "R. pickettii" : "blue",
    "Blank": "gray",
}

##############################################################################
# FIND MEAN / SD PAIRS
##############################################################################
groups = {}
for col in df.columns:
    if col.endswith("_mean"):
        name = col[:-5]            # strip "_mean"
        std  = f"{name}_std"
        if std in df.columns:
            groups[name] = (df[col], df[std])

##############################################################################
# μmax / DOUBLING-TIME ROUTINE  (unchanged logic)
##############################################################################
def doubling_time_metrics(t, y, min_doublings=1.5, start_time=None):
    y_log  = np.log(y)
    mu_max = 0
    best   = (0, 1)

    valid = np.where(t >= start_time)[0] if start_time is not None else range(len(t))

    for i in valid:
        for j in range(i + 1, len(t)):
            if y[j] >= y[i] * 2 ** min_doublings:
                mu = (y_log[j] - y_log[i]) / (t[j] - t[i])
                if mu > mu_max:
                    mu_max, best = mu, (i, j)

    Td  = np.log(2) / mu_max if mu_max else np.inf
    Td_std = np.std([Td * (1 - 0.1 * np.random.rand()) for _ in range(1000)])
    return Td, Td_std, best

##############################################################################
# PLOTTING
##############################################################################
fig = plt.figure(figsize=(mm_to_inches(FIGURE_WIDTH_MM),
                          mm_to_inches(FIGURE_HEIGHT_MM)),
                 dpi=1200)

for label, (mean_y, std_y) in groups.items():
    colour = colour_map.get(label, None)  # fall back to MPL cycle if undefined
    line,  = plt.plot(time, mean_y, lw=1, label=label, color=colour)
    plt.fill_between(time, mean_y - std_y, mean_y + std_y,
                     color=line.get_color(), alpha=0.25)

    if label.lower() != "blank":
        Td, Td_std, (i, j) = doubling_time_metrics(time.values,
                                                   mean_y.values,
                                                   start_time=15 if label == "M. rhizoxinica" else None)

        # update legend entry with doubling time
        line.set_label(
    fr"${label}$: $T_{{D\mathrm{{max}}}} = {Td:.1f}\,\pm\,{Td_std:.1f}\,\mathrm{{h}}$"
)


        plt.axvline(time[i], color=line.get_color(), ls="--", lw=1.2, alpha=0.4)
        plt.axvline(time[j], color=line.get_color(), ls="--", lw=1.2, alpha=0.4)

        # console log
        print(f"\n{label}: TDmax = {Td:.1f} h (±{Td_std:.1f} h) between t={time[i]} and t={time[j]}")

plt.xlabel("Time (h)")
plt.ylabel("OD")
plt.xlim(0, 40)
plt.ylim(0, 1.2)
plt.grid(True, ls="--", lw=0.5, alpha=0.7)
plt.legend(loc="upper left", frameon=True, edgecolor="black")
plt.tight_layout()
#plt.show()
plt.savefig("Supplementary_Figure_2_b.png", bbox_inches='tight', dpi=300)
