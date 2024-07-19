# Welcome to ZEPPI 
## Overview

**ZEPPI** (***Z-score Evaluation of Protein-Protein Interfaces***) is a method designed for proteome-scale sequence-based evaluation of protein-protein interfaces as defined by 3D atomic-level models of protein-protein interaction (PPI) direct binary complexes. This method was first introduced in [this paper](https://www.pnas.org/doi/abs/10.1073/pnas.2400260121). 

Structural models for the PPI complexes may be derived from a number of sources, including:

- Experimentally determined structure in [The Protein Data Bank (PDB)](https://www.rcsb.org)
- Template-based modeling as in [PrePPI](https://honiglab.c2b2.columbia.edu/PrePPI/)
- Protein docking studies
- Deep learning algorithms such as those based on AlphaFold Multimer. 

ZEPPI addresses the evaluation challenge by capitalizing on sequence co-evolution and conservation of the residues that come into interfacial contact. The hallmark of ZEPPI lies in its computation of the interface Z-score, **the ZEPPI score**, achieved through the comparison of interface-derived contact metrics against metrics derived from randomly selected surface residues. Thus, ZEPPI negates the need for factoring in indirect interactions, as the structural model delineates the interacting residues.

This repository contains code for ZEPPI and the needed input data for the tutorial examples. ZEPPI can be run either on a PC for a small number of protein-protein structural models or, for larger-scale explorations, as batched jobs on SGE- or SLURM-equipped clusters.

## Setup
ZEPPI is implemented in Python 3 and requires the following modules: ***biopython, numpy, numba, scipy, pandas***, and the below published packages:

- [**Surfv**](https://honig.c2b2.columbia.edu/surface-algorithms) ([github](https://github.com/honig-lab/SURFace-Algorithms)) to calculate solvent accessible surface area for protein structures
- [**HMMER**](http://hmmer.org/) to make multiple sequence alignments for sequence homologs of a query sequence
- [**HH-suite**](https://vogdb.org/research/hh-suite) to search for protein sequences similar to a query sequence in protein sequence databases

The installation time varies but typically should not exceed 30 minutes. This program has been tested on MacOS 12.6.7 and 13.4.1, Springdale Linux 7.9 (Verona), and Rocky Linux 8.5 (Green Obsidian) with bash-4.4.20, python-3.9.7, numpy-1.23.4, biopython-1.79, numba-0.56.2, scipy-1.9.1, and pandas-1.5.0.


## Run ZEPPI with examples

### Step 1. Download ZEPPI

Download the zipped GitHub codebase and unzip it. To run the tutorial examples, make sure the above-mentioned Python modules are installed on your machine.


### Step 2. Run ZEPPI
Go to the unzipped folder, and configure your directory path and python path in the *`Run_ZEPPI.sh`* file. Run the command with *`bash`* to see the needed arguments and available options.

```properties
bash Run_ZEPPI.sh
```

```properties
Usage: bash Run_ZEPPI.sh PPI_list.csv Output.csv -option
where:
    PPI_list  A csv file containing your query PPIs.
    -m  calculate ZEPPI on mutual information & conservation (recommended for heterodimers).
    -md calculate ZEPPI on mutual information & conservation and direct coupling analysis (recommended for homodimers).
Output files will be named after the input file with different endings.
```
The parameters are:
- `PPI_list`  is the input *csv* file containing the PPI name, chain names, and UniProt identifiers.
- `Output`  is output ZEPPI metrics calculated for the input PPI list.
- `-m`  is the option to calculate ZEPPI based on mutual information and conservation (recommended for heterodimers).
- `-md` is the option to calculate ZEPPI based on mutual information, conservation, and direct coupling analysis (recommended for homodimers; memory-consuming).
- Output files will be named after the `PPI_list` file with different endings.

As a demonstration, the *`Demo`* folder provides the input files for the example PPI_list that contains ten bacterial PDB dimer structures. Try to test ZEPPI in the *`Scratch`* folder with the below command.  The expected running time is ~2 min for the `-m` option or ~8 min for the `-md` option. The expected output files of both option `-m` or `-md` are provided in the *`Demo`* folder. Detailed score files are in the  *`Demo/Metrics`* folder. To run ZEPPI for your own structure models, the required input files are described in the `Run_ZEPPI.sh` file and the needed codes are stored in the *`Methods`* folder.

```properties
bash Run_ZEPPI.sh $YOUR_PATH/Demo/bacteria_PDBdimer_demo.csv $YOUR_PATH/Scratch/bacteria_PDBdimer_demo_ZEPPI_m.csv -m
```

### Step 3. Understand ZEPPI results

In the above example, the final output *bacteria_PDBdimer_demo_ZEPPI.csv* contains the `ZEPPI` score (last column). The detailed results are described in the below:

| Column    | Meaning |
| -------- | ------- |
| PPI  | the name of PPI  |
| N_MSA | the MSA depth   |
| N_IFR | the number of interface residues   |
| Zmean_MI    | Z-score using the mean Mutual Information of all the interface contacts   |
| Zmean_Con   | Z-score using the mean Conservation score of all the interface contacts  |
| Zmean_DCA   | Z-score using the mean DCA score of all the interface contacts  |
| Ztop_MI    | Z-score using the top Mutual Information of all the interface contacts  |
| Ztop_Con   | Z-score using the top Conservation score of all the interface contacts |
| Ztop_DCA   | Z-score using the mean DCA score of all the interface contacts  |
| ZEPPI   | the largest Z-score among all the above metrics  |

For more details, please read the ZEPPI paper. 

## Citation
Zhao, Haiqing, et al. "ZEPPI: Proteome-scale sequence-based evaluation of protein–protein interaction models." Proceedings of the National Academy of Sciences 121.21 (2024): e2400260121.

```latex
@article{zhao2024zeppi,
  title={ZEPPI: Proteome-scale sequence-based evaluation of protein--protein interaction models},
  author={Zhao, Haiqing and Petrey, Donald and Murray, Diana and Honig, Barry},
  journal={Proceedings of the National Academy of Sciences},
  volume={121},
  number={21},
  pages={e2400260121},
  year={2024},
  publisher={National Acad Sciences}
}
```
## Contact

Haiqing Zhao (<hz2592@cumc.columbia.edu>) or Barry Honig (<bh6@cumc.columbia.edu>)


