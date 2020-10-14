# SMILESMerge

SMILESMerge is a cheminformatic, Python program for that crosses two parent ligands (SMILES string) to create a child compound. SMILESMerge uses the crossover algorithm from the program AutoGrow4.

SMILESMerge accepts a list of SMILES (in .smi format) and will attempt to produce a user defined number of unique children SMILES. SMILESMerge also allows users to apply chemical property filters. Additionally, SMILESMerge options the automated conversion of the children from 1/2D SMILES to 3D representations (PDB and 3D-SDF) using Gypsum-DL.

![Image of crossover](https://github.com/Jacob-Spiegel/SMILESMerge/blob/main/Figures/crossover.png)

# Using SMILESMerge

SMILESMerge is primarily a python based program which can be run through command-line and GUI interface. More details can be found within `/PATH/SMILESMerge/tutorial/tutorial.md` or by running:

`python /PATH/SMILESMerge/RunSMILESMerge.py -h`

## Command-Line Interface

In a terminal with a modern Python environment run either:

`python /PATH/SMILESMerge/RunSMILESMerge.py -j /PATH/To/JSON_Paramerter_File.json`

or

`python /PATH/SMILESMerge/RunSMILESMerge.py --source_compound_file  /PATH/To/SMILES_File.smi`

Where `/PATH/To/JSON_Paramerter_File.json` is a JSON file containing all user parameters and `/PATH/To/SMILES_File.smi` is a file of all SMILES to be crossed


## GUI Interface

The GUI interface requires the additional dependency of GOOEY (https://github.com/chriskiehl/Gooey). To use the GUI, please run the following command from a terminal with a modern Python environment:

`python /PATH/SMILESMerge/RunSMILESMerge_GUI.py `


# Dependencies

SMILESMerge requires the following dependencies:
- RDKit
- func_timeout
- scipy
- numpy
- mpi4py (optional but required for MPI multiprocessing)
- Gooey (optional but required for GUI interface)

SMILESMerge was tested on a Python 3.7.7 environment on a Ubuntu OS. A docker image is provided for users with unsupported operating systems, such as Windows.

# Acknowledgment and Citing

Much of this code is take directly and/or adapted from AutoGrow4. This program also relies on Gypsum-DL and Dimorphite-DL for ligand handling and multiprocessing. Please remember to cite the following papers:

## SMILESMerge Citation:

- Spiegel, J.O. SMILESMerge: an open-source program for automated de novo ligand design using a crossover method. (2020)

## GlauconiteFilter Citation:

- Spiegel, J.O. GlauconiteFilter: an open-source program for automated ADME-PK filtering. (2020) https://doi.org/10.5281/zenodo.4087648
[![DOI](https://zenodo.org/badge/303535253.svg)](https://zenodo.org/badge/latestdoi/303535253)

## AutoGrow4 Citation:

- Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm for de novo drug design and lead optimization. J Cheminform 12, 25 (2020). https://doi.org/10.1186/s13321-020-00429-4

## Gypsum-DL Citation:

- Ropp, P.J., Spiegel, J.O., Walker, J.L. et al. Gypsum-DL: an open-source program for preparing small-molecule libraries for structure-based virtual screening. J Cheminform 11, 34 (2019). https://doi.org/10.1186/s13321-019-0358-3

## AutoGrow4 Dimorphite-DL:

- Ropp, P.J., Kaminsky, J.C., Yablonski, S. et al. Dimorphite-DL: an open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11, 14 (2019). https://doi.org/10.1186/s13321-019-0336-9

