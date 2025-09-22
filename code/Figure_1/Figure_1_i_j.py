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
from scipy.stats import ttest_ind


#############
# Panel 1 i #
#############

# Define the file path
file_path = "sourcedata_Figure_1_i.xlsx"

# Load the Excel file into a pandas DataFrame
data = pd.read_excel(file_path)

# Define pairs for statistical comparison
pairs = [('B+', 'B++')]  # Example pair names; replace with actual conditions

# Convert plot size to mm (1 inch = 25.4 mm)
width_mm = 30  # Adjusted width for better scaling
height_mm = 40  # Adjusted height for better scaling
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Calculate Y-axis max value
max_y = data['Area_red'].max()  # Get the max value
y_padding = 5000  # Add some padding above violins
y_max = int((max_y + y_padding) // 2500) * 2500  # Round up to the nearest 25000

# Initialize the plot
plt.figure(figsize=(width_in, height_in), dpi=2400)

# Create violin plot with quartiles
sns.violinplot(
    data=data,
    x='Condition',
    y='Area_red',
    hue='Condition',
    inner='quartile',  # Include quartiles inside the violin
    density_norm='width',  # Scale violins by dataset size
    palette=['orange', 'purple'],  # Custom colors
    linewidth=1  # Edge width of violins
)

# Overlay individual data points
sns.stripplot(
    data=data,
    x='Condition',
    y='Area_red',
    color='black',  # Points in black
    jitter=True,  # Add jitter to avoid overlap
    size=1.5,  # Point size
    alpha=0.8  # Transparency
)

# Adding statistical annotation
annotator = Annotator(
    plt.gca(), pairs, data=data, x='Condition', y='Area_red'
)
annotator.configure(test='t-test_ind', text_format='star', loc='outside', fontsize=6)
#annotator.apply_and_annotate()

# Adjust y-axis ticks and range
plt.ylim(-1000, y_max + 500)  # Ensure space above violins for annotations

# Set custom y-axis ticks
y_ticks = range(0, y_max + 2500, 5000)  # Adjust tick intervals as needed
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

# Show the plot
#plt.show()
plt.savefig("Figure_1i.png",  bbox_inches='tight', dpi=300)
plt.close()




# Extract the two groups
group1 = data[data['Condition'] == 'B+']['Area_red']
group2 = data[data['Condition'] == 'B++']['Area_red']

# Perform independent t-test
t_stat, p_val = ttest_ind(group1, group2, equal_var=False)  # Welch's t-test

# Print result
print(f"Comparison: ('B+', 'B++'), p-value: {p_val:.4e}, t-statistic: {t_stat:.4f}")

#############
# Panel 1 j #
#############

# Define the file path
file_path = "sourcedata_Figure_1_j.xlsx"

# Load the Excel file into a pandas DataFrame
data = pd.read_excel(file_path)

# Convert germination rates to percentages
data['germination rate day 2'] *= 100

# Define custom colors for each condition
custom_colors = {
    "Sorted B+": "orange",
    "Sorted B++": "purple"
}

# Convert plot size to mm (1 inch = 25.4 mm)
width_mm =30  # Adjusted width for better scaling
height_mm = 42  # Adjusted height for better scaling
width_in = width_mm / 25.4
height_in = height_mm / 25.4

# Initialize the plot
plt.figure(figsize=(width_in, height_in), dpi=2400)

# Create a bar plot for mean germination rate per condition
sns.barplot(
    data=data,
    x='Condition',
    y='germination rate day 2',
    errorbar=None,  # No error bars
    hue='Condition',
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

# Show the plot
#plt.show()
plt.savefig("Figure_1j.png", bbox_inches='tight', dpi=300)
plt.close()
