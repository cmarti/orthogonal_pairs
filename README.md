# Evolutionary paths that link orthogonal pairs of binding proteins

This repository contains scripts to reproduce the analysis of the computational binding landscape comprising all possible combination of mutations that differentiate orthogonal colicin toxin-antitoxin colicin systems E2/Im2 and E9/Im9 at their binding interface.

### Requirements

- [Rosetta](https://rosettacommons.org/software/). Version?
- [gpmap_tools](https://github.com/cmarti/gpmap-tools). Run with commit id 54eadf1048acea529402e4a9ab1ac69fd1ec01e5

### Installation

Create and activate new python environment
```bash
conda create -n orthogonal_pairs python=3.8
conda activate orthogonal_pairs
```

Install dependencies
```bash
git clone https://github.com/cmarti/gpmap-tools.git
cd gpmap-tools
git checkout 54eadf1048acea529402e4a9ab1ac69fd1ec01e5
pip install -r requirements.txt
python setup.py install
cd ..
```

To install Rosetta ...



Clone repository
```bash
git clone https://github.com/cmarti/orthogonal_pairs.git
cd orthogonal_pairs
```

Add repository dir to PYTHONPATH
```bash
path=`pwd`
export PYTHONPATH=$PYTHONPATH:$path
```

In total, the process is expected to take about 2 minutes

### Input data

Input data comprises is provided in the `data` folder and comprises two files:

- `Input PDB files for the scripts below`
- `colicins.pq`: a parquet file containing the predicted binding energies for every possible pair computed in the context of both E2/Im2 E9/Im9 crystal structures
- `experimental_data.csv`: a csv file containing experimental measurements for the 110 pairs used to calibrated the computationally predicted energies


### Computing E/Im energies with Rosetta

To create energetic calculations for structures of binding proteins:

First, create refined model with `Refinement.xml` and `flag_ref` :

```bash 
~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol rosetta_xmls/Refinement.xml -s <path_to_PDB_file> @rosetta_xmls/flag_ref -parser:script_vars cst_full_path=<path_to_coordinate_constraint_file> -nstruct 20
```

Second, create the model with the wanted mutations using `Making_model.xml` and `flags_Mm`:
```bash
~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @rosetta_xmls/flags_Mm
```

### Visualizing the E/Im binding landscape

Scripts to reproduce each figure can be found in the `scripts` directory one by one or running the following bash scripts

To make calculations required to plot the visualizations.

```bash
bash 01_make_calculations.sh
```

To reproduce only the figures once the computations are done run the following script

```bash
export MPLBACKEND='agg'
bash 02_make_figures.sh
```

### Expected running times

Rosetta was run in ...
- Model refinement: 
- Energy computatin in mutant sequences: 


Visualization scripts were run in Ubuntu 20.04.6 OS with a Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz processor and required about 12GB of memory. Running times were:
- Calculations:   19 minutes 16.110 seconds
- Plots:           9 minutes  6.902 seconds


### Expected output

- Calculation of binding energies of the E/Im pairs on the two backbones is stored at `Output file from running Rosetta`
- Intermediate files required to reproduce the figures will be stored in the `visualization` folder.
- Individual panels will be saved at the `figures` folder, which were assembled manually into the different figures using [Inkscape](https://github.com/cmarti/gpmap-tools)


### Citation

Ziv Avizemer, Carlos Martí‐Gómez, Shlomo Yakir Hoch et al. Evolutionary paths that link orthogonal pairs of binding proteins, 20 April 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2836905/v1]
