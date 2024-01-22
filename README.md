# Evolutionary paths that link orthogonal pairs of binding proteins

This repository contains scripts to reproduce the analysis of the computational binding landscape comprising all possible combination of mutations that differentiate orthogonal colicin toxin-antitoxin colicin systems E2/Im2 and E9/Im9 at their binding interface.

### Requirements

- [gpmap_tools](https://github.com/cmarti/gpmap-tools). Run with commit id e57dee1121b2142d4f9b85d3e35f3edf60f40fcd

### Installation

Create and activate new python environment
```bash
conda create -n orthogonal_pairs python=3.8
conda activate orthogonal_pairs
```

Install dependencies
```bash
git clone https://github.com/cmarti/gpmap-tools.git
cd gpmap_tools
git checkout e57dee1121b2142d4f9b85d3e35f3edf60f40fcd
python setup.py install
cd ..
```

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

### Input data

Input data comprises is provided in the `data` folder and comprises two files:

- `colicins.pq`: a parquet file containing the predicted binding energies for every possible pair computed in the context of both E2/Im2 E9/Im9 crystal structures
- `experimental_data.csv`: a csv file containing experimental measurements for the 110 pairs used to calibrated the computationally predicted energies

### Scripts

Scripts to reproduce each figure can be found in the `scripts` directory one by one or just by running

```bash
bash make_figures.sh
```

### Output figures

Individual panels will be saved at the `figures` folder, which were assembled manually into the different figures using [Inkscape](https://github.com/cmarti/gpmap-tools)

### Citation

Ziv Avizemer, Carlos Martí‐Gómez, Shlomo Yakir Hoch et al. Evolutionary paths that link orthogonal pairs of binding proteins, 20 April 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2836905/v1]
