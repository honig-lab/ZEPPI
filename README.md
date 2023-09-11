# Welcome to ZEPPI 
## Overview

**ZEPPI** (***Z-score Evaluation of Protein-Protein Interfaces***) is a method designed for proteome-scale sequence-based evaluation of protein-protein interfaces as defined by 3D atomic-level models of protein-protein interaction (PPI) direct binary complexes. 

Structural models for the PPI complexes may be derived from a number of sources, including:

- Experimentally determined structure in [The Protein Data Bank (PDB)](https://www.rcsb.org)
- Template-based modeling as in [PrePPI](https://honiglab.c2b2.columbia.edu/PrePPI/)
- Protein docking studies
- Deep learning algorithms such as those based on AlphaFold Multimer. 

ZEPPI addresses the evaluation challenge by capitalizing on sequence co-evolution and conservation of the residues that come into interfacial contact. The hallmark of ZEPPI lies in its computation of the interface Z-score, **the ZEPPI score**, achieved through the comparison of interface-derived contact metrics against metrics derived from randomly selected surface residues. Thus, ZEPPI negates the need for factoring in indirect interactions, as the structural model delineates the interacting residues.

This repository contains code for ZEPPI and the needed input data for the tutorial examples. ZEPPI can be run either on a PC for a small number of protein-protein structural models or, for larger-scale explorations, as batched jobs on SGE- or SLURM-equipped clusters.

## Setup

ZEPPI is implemented in Python 3 and requires the following libraries: *biopython, numpy, numba, scipy, pandas,* and the below published packages:

- [**Surfv**](https://honig.c2b2.columbia.edu/surface-algorithms) to calculate solvent accessible surface area for protein structures
- [**HMMER**](http://hmmer.org/) to make multiple sequence alignments for sequence homologs of a query sequence
- [**HH-suite**](https://vogdb.org/research/hh-suite) to search for protein sequences similar to a query sequence in protein sequence databases

The install time varies but typically should not exceed 30 minutes. This program has been tested on MacOS 12.6.7 and 13.4.1, Springdale Linux 7.9 (Verona), and Rocky Linux 8.5 (Green Obsidian) with bash-4.4.20, python-3.9.7, numpy-1.23.4, biopython-1.79, numba-0.56.2, scipy-1.9.1, and pandas-1.5.0.


## Run ZEPPI with examples

### Step 1. Download ZEPPI

Download the zipped GitHub codebase and unzip it. To run the tutorial examples, make sure the above-mentioned Python modules are installed on your machine.


### Step 2. Run ZEPPI
Go to the unzipped folder, and edit the *`Run_ZEPPI.sh`* file with the *ZEPPI_base* variable being your directory path. Run it with *`bash`* to see the usage and available options for running ZEPPI.

```properties
bash Run_ZEPPI.sh
```

```properties
Usage: bash Run_ZEPPI.sh PPI_list.csv Output.csv -option
```
The parameters are:
- `PPI_list`  Points to the *csv* file containing the PPI name, chain names, and UniProt identifiers of your query PPIs.
- `Output`  Contains the computed ZEPPI metrics for the input PPI list.
- `-m`  To calculate ZEPPI based on mutual information and conservation (recommended for heterodimers).
- `-md` To calculate ZEPPI based on mutual information, conservation, and direct coupling analysis (recommended for homodimers; memory-consuming).
- Output files will be named after the `PPI_list` file with different endings.

In the *`Demo`* folder, we provided the input and output files of 10 bacterial dimers as an example. Try to run the below command in the *`Scratch`* folder, and compare your result with the provided output from the *`Demo`* folder. To run ZEPPI for your own structure models, the input files required and their formats are described in the `Run_ZEPPI.sh` file provided with the codebase. The expected run time for the demo PPIs is about 2 min for the `-m` option or 8 min for the `-md` option on a standard desktop computer. The expected output files are stored in the *`Demo`* folder with details in the  *`Demo/Metrics`* folder.

```properties
bash Run_ZEPPI.sh ./Demo/bacteria_PDBdimer_demo.csv ./Scratch/bacteria_PDBdimer_demo_ZEPPI_m.csv -m
```

The final output file *bacteria_PDBdimer_demo_ZEPPI.csv* contains the `ZEPPI` score (last column) and other columns as described in the below:

| Column    | Meaning |
| -------- | ------- |
| PPI  | the name of PPI  |
| N_MSA | the MSA depth   |
| Zmean_MI    | Z-score using the mean Mutual Information of all the interface contacts   |
| Zmean_Con   | Z-score using the mean Conservation score of all the interface contacts  |
| Ztop_MI    | Z-score using the top Mutual Information of all the interface contacts  |
| Ztop_Con   | Z-score using the top Conservation score of all the interface contacts |

For more details, please read our paper. 

## Citation
## Contact

Haiqing Zhao (<hz2592@cumc.columbia.edu>) or Barry Honig (<bh6@cumc.columbia.edu>)


