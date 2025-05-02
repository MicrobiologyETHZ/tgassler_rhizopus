import yaml
import pandas as pd
import json
import os


def load_figure_data(figure_name, panel_name, index_col=None):
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
            return pd.read_csv(data_file, delimiter=delimiter, index_col=index_col)
        
        elif file_extension == '.json':
            with open(data_file, 'r') as f:
                return json.load(f)
            
        elif file_extension in ['.xls', '.xlsx']:
            return pd.read_excel(data_file, index_col=index_col)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")
            
    except Exception as e:
        raise IOError(f"Error loading data from {data_file}: {str(e)}")


import yaml
import os
import matplotlib.pyplot as plt

def save_figure_panel(figure_name, panel_name,format='svg', dpi=300):
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures_config.yaml")
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    if 'figures' not in config or figure_name not in config['figures']:
        raise ValueError(f"Figure '{figure_name}' not found in configuration")
    if panel_name not in config['figures'][figure_name]:
        raise ValueError(f"Panel '{panel_name}' not found in figure '{figure_name}'")
    panel_info = config['figures'][figure_name][panel_name]
    figure_file = panel_info.get('figure_file') + "." + format
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


# This sod can be used to retrieve the AA sequences from the annotation file
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Depending on which df you are looking at, defines which ones will be extracted

#df = overlap_df.nlargest(10, "log2FoldChange_RP1")
#df = overlap_df.nsmallest(10, "log2FoldChange_RP1")
#df = overlap_df #consider which one you are working with

def get_deg_gene_seqs(df, fasta_path="braker.fasta", output_fasta_path="Overlap_Down.fasta"):

    # "Top10_Up.fasta"; "Top10_Down.fasta";  "Overlap_Up.fasta"; "Overlap_Down.fasta"
    # Parse the FASTA file
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
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