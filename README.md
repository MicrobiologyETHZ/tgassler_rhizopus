# Induced endosymbiosis between Rhizopus microsporus and Ralstonia pickettii indicates a shift from antagonism to commensalism 
### Accompanying analyses scripts

## Citation
This repository contains a collection of Python scripts for generating figure panels for the following publication: 

```
Gassler, T., ... (2025). Induced endosymbiosis between Rhizopus microsporus and Ralstonia pickettii indicates a shift from antagonism to commensalism. Journal Name, Volume(Issue), Pages.
DOI: https://doi.org/xx.xxxx/xxxxx

```

## Requirements

All required packages are specified in the `environment.yml` file. The scripts were developed and tested with Python 3.13.


## Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/MicrobiologyETHZ/tgassler_rhizopus.git
   cd tgassler_rhizopus
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate rhizopus-env
   ```
   
   Alternatively, you can use mamba for faster dependency resolution:
   ```bash
   mamba env create -f environment.yml
   mamba activate rhizopus-env
   ```

## Configuration

The `code/figures_config.yml` file contains the mapping between scripts, input data files, and output figure names. It specifies which data files are required for each script and the naming conventions for the output figures. Please make sure the paths specifed in the config file are pointing to the correct data files.


## Usage

You can run individual scripts to generate specific figure panels:

```bash
python code/Figure_1_i_j.py
```

Or run all scripts sequentially to generate all figures:

```bash
for script in code/*.py; do python "$script"; done
```

## License

This project is licensed under the GNU GENERAL PUBLIC License - see the [LICENSE](LICENSE) file for details.