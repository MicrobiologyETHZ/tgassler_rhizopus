import yaml
import pandas as pd
import json
import os


def load_figure_data(figure_name, panel_name):
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures_config.yaml")
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    if 'figures' not in config or figure_name not in config['figures']:
        raise ValueError(f"Figure '{figure_name}' not found in configuration")
    
    if panel_name not in config['figures'][figure_name]:
        raise ValueError(f"Panel '{panel_name}' not found in figure '{figure_name}'")
    
    panel_info = config['figures'][figure_name][panel_name]
    data_file = panel_info.get('data_file')
    if not data_file:
        raise ValueError(f"No data file specified for {figure_name}/{panel_name}")

    file_extension = os.path.splitext(data_file)[1].lower()
    try:
        if file_extension in ['.csv', '.tsv']:
            delimiter = ',' if file_extension == '.csv' else '\t'
            return pd.read_csv(data_file, delimiter=delimiter)
        
        elif file_extension == '.json':
            with open(data_file, 'r') as f:
                return json.load(f)
            
        elif file_extension in ['.xls', '.xlsx']:
            return pd.read_excel(data_file)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")
            
    except Exception as e:
        raise IOError(f"Error loading data from {data_file}: {str(e)}")


import yaml
import os
import matplotlib.pyplot as plt

def save_figure_panel(figure_name, panel_name, format='svg', dpi=300):
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures_config.yaml")
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    if 'figures' not in config or figure_name not in config['figures']:
        raise ValueError(f"Figure '{figure_name}' not found in configuration")
    if panel_name not in config['figures'][figure_name]:
        raise ValueError(f"Panel '{panel_name}' not found in figure '{figure_name}'")
    panel_info = config['figures'][figure_name][panel_name]
    figure_file = panel_info.get('figure_file')
    if not figure_file:
        raise ValueError(f"No data file specified for {figure_name}/{panel_name}")
    output_dir = os.path.dirname(figure_file)
    os.makedirs(output_dir, exist_ok=True)
    try:
        # Save the current matplotlib figure
        plt.savefig(figure_file, format=format, bbox_inches='tight', dpi=dpi)
        plt.close()  # Close the figure to free memory
        return figure_file
    except Exception as e:
        raise IOError(f"Error saving figure to {figure_file}: {str(e)}")


# Example usage
if __name__ == "__main__":
    import numpy as np
    
    # Create a sample plot
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, 'b-')
    plt.title('Sample Sine Wave')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    
    try:
        # Save the figure
        output_path = save_figure_panel("Figure_1", "panel_a")
        print(f"Figure saved to: {output_path}")
    except Exception as e:
        print(f"Error: {e}")