# Evolutionary paths that link orthogonal pairs of binding proteins

This repository contains scripts to reproduce the analysis of the computational binding landscape comprising all possible combination of mutations that differentiate orthogonal colicin toxin-antitoxin colicin systems E2/Im2 and E9/Im9 at their binding interface.

### Requirements

- [Rosetta](https://github.com/RosettaCommons/rosetta). Run with commit id 1cdfc600cfe87cb1f733259dbc36c6734721770d
- [gpmap_tools](https://github.com/cmarti/gpmap-tools). Run with commit id 54eadf1048acea529402e4a9ab1ac69fd1ec01e5

### Installation

To install Rosetta, follow these [instructions](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build)

```bash
git clone https://github.com/RosettaCommons/rosetta
cd rosetta
git checkout 1cdfc600cfe87cb1f733259dbc36c6734721770d
cd source

cd ../..
```


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

Input data comprises is provided in the `data` folder and comprises a series of files:

- Input files for rosetta modeling
	- PDB files: `E2Im2_3u43.pdb.gz` and `E9Im9_1emv.pdb.gz`
	- Constraint files: `cst_E2Im2_3u43` and `cst_E9Im9_1emv`
- `colicins.pq`: a parquet file containing the predicted binding energies for every possible pair computed in the context of both E2/Im2 E9/Im9 crystal structures
- `experimental_data.csv`: a csv file containing experimental measurements for the 110 pairs used to calibrated the computationally predicted energies


### Computing E/Im energies with Rosetta

To create energetic calculations for structures of binding proteins:

First, we create a refined model of each wild-type pair, E2/Im2 and E9/Im9 with `Refinement.xml` and `flag_ref` :

- input structure: PDB file
- coordinate constriant file: cst file for the PDB you running

```bash 
~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol rosetta_xmls/Refinement.xml -s data/E2Im2_3u43.pdb.gz @rosetta_xmls/flag_ref -parser:script_vars cst_full_path=data/cst_E2Im2_3u43 -nstruct 20

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol rosetta_xmls/Refinement.xml -s data/E9Im9_1emv.pdb.gz @rosetta_xmls/flag_ref -parser:script_vars cst_full_path=data/cst_E9Im9_1emv -nstruct 20

```

These scripts will generate a series of files ...



Second, create the model with the desired mutations using `Making_model.xml` and `flags_Mm`:

- input structure: refined model with lowest score
- coordinate constraint file: same as above
- residues to fix: critical residues for function (for colicins binding: 54A,55A,86B)

For each combination of mutations from the resfile the following command should be created:
```bash
~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @rosetta_xmls/flags_Mm -parser:script_vars target0=24A new_res0=ARG pac0=true target1=26A new_res1=GLU pac1=true target2=27A new_res2=GLY pac2=true target3=28A new_res3=ALA pac3=true target4=29A new_res4=THR pac4=true target5=32A new_res5=ASP pac5=true target6=33A new_res6=ASP pac6=true target7=34A new_res7=ASN pac7=true target8=38A new_res8=ARG pac8=true target9=58A new_res9=ASP pac9=true target10=72B new_res10=LYS pac10=true target11=73B new_res11=PRO pac11=false target12=77B new_res12=SER pac12=false target13=78B new_res13=ASN pac13=true target14=83B new_res14=LYS pac14=true target15=97B new_res15=LYS pac15=false target16=98B new_res16=ARG pac16=true
```

More examples for more combinations:

```bash
~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @<path_to_flag_file> -parser:script_vars target0=24A new_res0=ARG p
ac0=true target1=26A new_res1=GLU pac1=true target2=27A new_res2=GLY pac2=true target3=28A new_res3=ALA pac3=true target4=29A new_res
4=THR pac4=true target5=32A new_res5=ASP pac5=true target6=33A new_res6=ASP pac6=true target7=34A new_res7=ASN pac7=true target8=38A
new_res8=ARG pac8=true target9=58A new_res9=ASP pac9=true target10=72B new_res10=LYS pac10=true target11=73B new_res11=PRO pac11=fals
e target12=77B new_res12=SER pac12=false target13=78B new_res13=ASN pac13=true target14=83B new_res14=LYS pac14=true target15=97B new
_res15=GLU pac15=true target16=98B new_res16=MET pac16=false

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @<path_to_flag_file> -parser:script_vars target0=24A new_res0=ARG p
ac0=true target1=26A new_res1=GLU pac1=true target2=27A new_res2=GLY pac2=true target3=28A new_res3=ALA pac3=true target4=29A new_res
4=THR pac4=true target5=32A new_res5=ASP pac5=true target6=33A new_res6=ASP pac6=true target7=34A new_res7=ASN pac7=true target8=38A
new_res8=ARG pac8=true target9=58A new_res9=ASP pac9=true target10=72B new_res10=LYS pac10=true target11=73B new_res11=PRO pac11=fals
e target12=77B new_res12=SER pac12=false target13=78B new_res13=ASN pac13=true target14=83B new_res14=LYS pac14=true target15=97B new
_res15=GLU pac15=true target16=98B new_res16=VAL pac16=false

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @<path_to_flag_file> -parser:script_vars target0=24A new_res0=ARG p
ac0=true target1=26A new_res1=GLU pac1=true target2=27A new_res2=GLY pac2=true target3=28A new_res3=ALA pac3=true target4=29A new_res
4=THR pac4=true target5=32A new_res5=ASP pac5=true target6=33A new_res6=ASP pac6=true target7=34A new_res7=ASP pac7=false target8=38A
 new_res8=THR pac8=false target9=58A new_res9=GLU pac9=false target10=72B new_res10=ASN pac10=false target11=73B new_res11=GLY pac11=
true target12=77B new_res12=SER pac12=false target13=78B new_res13=SER pac13=false target14=83B new_res14=ASN pac14=false target15=97
B new_res15=LYS pac15=false target16=98B new_res16=VAL pac16=false

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @<path_to_flag_file> -parser:script_vars target0=24A new_res0=ARG p
ac0=true target1=26A new_res1=GLU pac1=true target2=27A new_res2=GLY pac2=true target3=28A new_res3=ALA pac3=true target4=29A new_res
4=THR pac4=true target5=32A new_res5=ASP pac5=true target6=33A new_res6=ASP pac6=true target7=34A new_res7=ASP pac7=false target8=38A
 new_res8=THR pac8=false target9=58A new_res9=GLU pac9=false target10=72B new_res10=ASN pac10=false target11=73B new_res11=GLY pac11=
true target12=77B new_res12=SER pac12=false target13=78B new_res13=SER pac13=false target14=83B new_res14=ASN pac14=false target15=97
B new_res15=LYS pac15=false target16=98B new_res16=MET pac16=false
```

Run each command to get the model and energetic terms for that mutational combination


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

Rosetta was run in (see below for the format, which processor and memory requirements). Running times were
- Model refinement: ~7 hrs for 20 refined models for each wild-type structure
- Energy computatin in mutant sequences for one model : ~2min for each E/IM pair in each backbone (sequence space includes >2M mutants)

Visualization scripts were run in Ubuntu 20.04.6 OS with a Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz processor and required about 12GB of memory. Running times were:
- Calculations:   19 minutes 16.110 seconds
- Plots:           9 minutes  6.902 seconds


### Expected output

- Calculation of binding energies of the E/Im pairs on the two backbones is stored at `Output file from running Rosetta`
- Intermediate files required to reproduce the figures will be stored in the `visualization` folder.
- Individual panels will be saved at the `figures` folder, which were assembled manually into the different figures using [Inkscape](https://github.com/cmarti/gpmap-tools)

### Citation

Ziv Avizemer, Carlos Martí‐Gómez, Shlomo Yakir Hoch et al. Evolutionary paths that link orthogonal pairs of binding proteins, 20 April 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2836905/v1]
