# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:57:27 2025

@author: tgassler
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from statannotations.Annotator import Annotator  # Statistical annotations
from utilities import load_figure_data, save_figure_panel

#####################
# Panel i
#####################

data = load_figure_data("Figure_1", "panel_i")
# Define pairs for statistical comparison
pairs = [('B+', 'B++')]  # Example pair names; replace with actual conditions

# Convert plot size to mm (1 inch = 25.4 mm)
width_mm = 30  # Adjusted width for better scaling
height_mm = 35  # Adjusted height for better scaling
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Calculate Y-axis max value
max_y = data['Background Normalized RawIntDen'].max()  # Get the max value
y_padding = 50000  # Add some padding above violins
y_max = int((max_y + y_padding) // 25000) * 25000  # Round up to the nearest 25000

# Initialize the plot
plt.figure(figsize=(width_in, height_in), dpi=1200)

# Create violin plot with quartiles
sns.violinplot(
    data=data,
    x='Condition',
    y='Background Normalized RawIntDen',
    inner='quartile',  # Include quartiles inside the violin
    scale='width',  # Scale violins by dataset size
    palette=['orange', 'purple'],  # Custom colors
    linewidth=1  # Edge width of violins
)

# Overlay individual data points
sns.stripplot(
    data=data,
    x='Condition',
    y='Background Normalized RawIntDen',
    color='black',  # Points in black
    jitter=True,  # Add jitter to avoid overlap
    size=1.5,  # Point size
    alpha=0.8  # Transparency
)

# Adding statistical annotation
annotator = Annotator(
    plt.gca(), pairs, data=data, x='Condition', y='Background Normalized RawIntDen'
)
annotator.configure(test='t-test_ind', text_format='star', loc='outside', fontsize=6)
annotator.apply_and_annotate()

# Adjust y-axis ticks and range
plt.ylim(-5000, y_max + 60000)  # Ensure space above violins for annotations

# Set custom y-axis ticks
y_ticks = range(0, y_max + 50000, 100000)  # Adjust tick intervals as needed
plt.yticks(ticks=y_ticks, fontsize=6)  # Explicitly set y-ticks and their font size

# Configure scientific notation for y-axis
ax = plt.gca()  # Get the current axis
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))  # Always use scientific notation
ax.yaxis.set_major_formatter(formatter)

# Adjust the font size of the scientific notation label
ax.yaxis.get_offset_text().set_fontsize(6)  # Set the font size of the ×10⁵ label

# Adjust x-axis ticks
plt.xticks(rotation=30, fontsize=6, ha='right')  # Rotate x-axis labels properly

# Remove axis labels
plt.xlabel("")  # Remove x-axis label
plt.ylabel("")  # Remove y-axis label

# Tighten layout to prevent overlap
plt.tight_layout()

save_figure_panel("Figure_1", "panel_i", format='png', dpi=300)
# Show the plot
#plt.show()

################################
# Panel j
################################

# Define the file path
#file_path = "R4_Load_Germ.xlsx"

# Load the Excel file into a pandas DataFrame
#data = pd.read_excel(file_path)
data = load_figure_data("Figure_1", "panel_j")
# Convert germination rates to percentages
data['germination rate day 2'] *= 100

# Define custom colors for each condition
custom_colors = {
    "Originating B++": "grey",
        "Sorted B+": "orange",
    "Sorted B++": "purple"
}

# Convert plot size to mm (1 inch = 25.4 mm)
width_mm =35  # Adjusted width for better scaling
height_mm = 40  # Adjusted height for better scaling
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Initialize the plot
plt.figure(figsize=(width_in, height_in), dpi=1200)

# Create a bar plot for mean germination rate per condition
sns.barplot(
    data=data,
    x='Condition',
    y='germination rate day 2',
    ci=None,  # No error bars
    palette=[custom_colors[condition] for condition in data['Condition'].unique()]
)

# Overlay individual data points
sns.stripplot(
    data=data,
    x='Condition',
    y='germination rate day 2',
    color='black',  # Points in black
    jitter=True,  # Add jitter to avoid overlap
    size=3,  # Point size
    alpha=0.8  # Transparency
)

# Adjust y-axis ticks and range
max_y = data['germination rate day 2'].max()  # Get the max value
y_padding = 5  # Add some padding above the bars
plt.ylim(0, max_y + y_padding)  # Ensure space above bars
plt.yticks(range(0, int(max_y) + int(y_padding) + 10, 20), fontsize=6)  # Y-axis ticks in steps of 10%
plt.ylabel("")  # Remove y-axis label

# Customize font sizes and labels
plt.xticks(rotation=30, fontsize=6, ha='right')  # Rotate x-axis labels
plt.xlabel("")  # Remove x-axis label

# Tighten layout to prevent overlap
plt.tight_layout()
save_figure_panel("Figure_1", "panel_j", format='png', dpi=300)
# Show the plot
#plt.show()


