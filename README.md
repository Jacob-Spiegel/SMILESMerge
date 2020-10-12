# SMILESMerge

SMILESMerge is a cheminformatic, Python program for that crosses two parent ligands (SMILES string) to create a child compound. SMILESMerge uses the crossover algorithm from the program AutoGrow4.

SMILESMerge accepts a list of SMILES (in .smi format) and will attempt to produce a user defined number of unique children SMILES. SMILESMerge also allows users to apply chemical property filters. Additionally, SMILESMerge options the automated conversion of the children from 1/2D SMILES to 3D representations (PDB and 3D-SDF) using Gypsum-DL.

![Image of crossover](https://github.com/Jacob-Spiegel/SMILESMerge/blob/main/Figures/crossover.png)

# Acknowledgment and Citing

Much of this code is take directly and/or adapted from AutoGrow4. This program also relies on Gypsum-DL and Dimorphite-DL for ligand handling and multiprocessing. Please remember to cite the following papers:

## SMILESMerge Citation:
    
- Citation Pending

## AutoGrow4 Citation:

- Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm for de novo drug design and lead optimization. J Cheminform 12, 25 (2020). [doi: 10.1186/s13321-020-00429-4]

## Gypsum-DL Citation:

- Ropp, P.J., Spiegel, J.O., Walker, J.L. et al. Gypsum-DL: an open-source program for preparing small-molecule libraries for structure-based virtual screening. J Cheminform 11, 34 (2019). https://doi.org/10.1186/s13321-019-0358-3

## AutoGrow4 Dimorphite-DL:

- Ropp, P.J., Kaminsky, J.C., Yablonski, S. et al. Dimorphite-DL: an open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11, 14 (2019). https://doi.org/10.1186/s13321-019-0336-9

