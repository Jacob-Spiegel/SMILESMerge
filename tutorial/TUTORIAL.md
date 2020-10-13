# Welcome to SMILESMerge

This document will break down how to run SMILESMerge. It will also cover what
dependencies are required.

Please note that there are several paths used in this tutorial that need to be
replaced with the proper paths on one's own system. These paths will be
indicated by a string of ALL_CAPS. For example (from Section: Installing
AutoGrow):

`cd /PATH_TO/DESIRED_DIR/`

You must replace `/PATH_TO/DESIRED_DIR/` with the path to the `SMILESMerge`
directory on your own system. On a Ubuntu OS, this may look like:

`cd /home/jacob/Documents/`

For brevity, we simplify the main `SMILESMerge` directory to `/SMILESMerge/`
throughout this tutorial. You may need to supply the `/PATH_TO/DESIRED_DIR/`
described above before `/SMILESMerge/`.

## Computer Requirements

SMILESMerge has been tested on Ubuntu 16.04 and higher, as well as MacOS 10.13
High Sierra. It has been verified to work on an HPC cluster using SMP
multithreading (RedHat Enterprise Server release 7.3 Maipo).

SMILESMerge has not been configured for Windows OS, but a script capable of
running SMILESMerge within a docker container on Windows can be found:

`/SMILESMerge/docker/SMILESMerge_in_docker.py`

This script should run on any docker-enabled machine, and should be capable of
multithreading. Details on running SMILESMerge within a docker container can be
found below, in the Section: Docker submission.

## Installing SMILESMerge

Download a copy of SMILESMerge from Jacob Spiegel's github at
[https://github.com/Jacob-Spiegel/SMILESMerge.git](https://github.com/Jacob-Spiegel/SMILESMerge.git).

You can also install SMILESMerge using the git clone command:

```bash
cd /PATH_TO/DESIRED_DIR/
git clone https://github.com/Jacob-Spiegel/SMILESMerge.git
```

## Dependencies

SMILESMerge has several dependencies that may need to be installed separately.

### Bash (Required)

A modern installation of bash is required to run SMILESMerge. SMILESMerge has been
tested using GNU bash, version 4.4.19. macOS and Linux come with Bash
preinstalled.

### Coreutils (Required For macOS)

Most Linux operating systems include the `timeout` tool (part of the
`coreutils` package) that SMILESMerge requires. Use on macOS requires the
separate installation of the `coreutils` package, available through
`homebrew`, which provides the equivalent `gtimeout` binary.

```bash
sudo brew install coreutils
```

### Python Installation (Required)

SMILESMerge is primarily written in python. A modern version of python can be
installed using `conda`:

- [https://docs.conda.io/projects/conda/en/latest/user-guide/install/](https://docs.conda.io/projects/conda/en/latest/user-guide/install/),
  or
- [http://www.python.org/getit/](http://www.python.org/getit/).

SMILESMerge has been tested with python 2.7, 3.6, and 3.7. Future support and
updates will focus on 3.7. We recommend using the most current version of
python available, 3.7 or newer.

### Python APIs (Required)

SMILESMerge also uses several python API libraries beyond those found in the
standard library. These must be installed separately. Most can be installed
via `conda` or `pip`.

#### Mandatory Installations

RDKit: Cheminformatic library. RDKit can be downloaded via `conda`/`pip`. To
install using `conda` use the command:

```bash
conda install -c rdkit rdkit
```

We use the following RDKit sub-libraries in SMILESMerge:

```python
import rdkit
from rdkit import RDLogger, Chem, DataStructs
from rdkit.Chem import MolSurf, Crippen, rdFMCS, Descriptors, AllChem, FilterCatalog,  Lipinski, rdDepictor
from rdkit.Chem.Draw import PrepareMolForDrawing, rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem.rdchem import BondStereo
```

NumPy (mathematical functions) can be downloaded via `conda`/`pip`. It can be
`conda` installed using the command `conda install -c anaconda numpy`.
SMILESMerge has been tested using `numpy` version 1.15.0.

SciPy (mathematical functions) can be downloaded via `conda`/`pip`. It can be
`conda` installed using the command `conda install -c anaconda scipy`.
SMILESMerge has been tested using `scipy` version 1.1.0.

func_timeout (pythonic timeout tool) can be downloaded via `pip`. It can be
`pip` installed using the command `pip install func-timeout`. SMILESMerge has
been tested using `func_timeout` version 4.3.5.

#### Optional Installations

Gooey API is required to run the GUI interface `/SMILESMerge/RunSMILESMerge_GUI.py`.
It can be obtained via `conda`/`pip`. Details for installation can be found at:
`https://github.com/chriskiehl/Gooey`


mpi4py (MPI multithreading python library) is required for MPI multithreading.
It can be downloaded via `conda`/`pip`. It can be `conda` installed using the
command `conda install -c anaconda mpi4py`. SMILESMerge has been tested using
`mpi4py` version 3.0.1. This may require a preinstallation of `mpich`: `sudo
apt install mpich`

SMILESMerge requires `mpi4py` version 2.1.0 and higher. To check the version:

1. open a python window.
2. enter into the window:

```python
>>> import mpi4py
>>> mpi4py.__version__
    3.0.1
```

MPI mode also requires an MPI-enabled computer environment. The authors use
OpenMPI. OpenMPI installation instructions can be found:
http://lsi.ugr.es/jmantas/pdp/ayuda/datos/instalaciones/Install_OpenMPI_en.pdf

Quick OpenMPI installation is possible in a bash terminal:

```bash
sudo apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.3 libopenmpi-dbg libopenmpi-dev
```

Establishing a fully MPI-enabled computer network is complicated and should
only be attempted by qualified technicians. The authors used an Intel’s
Omni-Path communication architecture that was established by experts at the
University of Pittsburgh’s Center for Research Computing. The authors DO NOT
RECOMMEND ATTEMPTING THIS ON YOUR OWN.

### Pre-Installed Python and Binary Dependencies

SMILESMerge comes with several dependencies preinstalled, requiring no
additional effort by the user. These packages have licenses that allow them to
be freely redistributed. If a dependency updates, please feel free to contact
us, and we will do our best to make our code future-compatible.

#### SMILES Conversion to 3D and Protonation Adjustments

SMILESMerge performs most of its ligand handling using 2D SMILES. SMILESMerge uses
the free and open-source program Gypsum-DL to convert from SMILES to 3D SDF
format. Gypsum-DL is prepackaged in SMILESMerge. Gypsum-DL itself also includes
the MolVS and Dimorphite-DL packages.

- Gypsum-DL:
  - Version: 1.1.2
  - Location: `/SMILESMerge/smilesmerge/operators/convert_files/gypsum_dl/`
  - Citation: Ropp PJ, Spiegel JO, Walker JL, Green H, Morales GA, Milliken
    KA, Ringe JJ, Durrant JD. Gypsum-DL: An Open-Source Program for Preparing
    Small-Molecule Libraries for Structure-Based Virtual Screening. J
    Cheminform. 11(1):34, 2019. [PMID: 31127411] [doi:
    10.1186/s13321-019-0358-3]
  - License: Apache version 2.0

- Dimorphite-DL:
  - Version: 1.2.2
  - Location:
    `/SMILESMerge/smilesmerge/operators/convert_files/gypsum_dl/gypsum_dl/Steps/SMILES/dimorphite_dl`
  - Citation: Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019)
    Dimorphite-DL: An open-source program for enumerating the ionization
    states of drug-like small molecules. J Cheminform 11:14.
    doi:10.1186/s13321-019-0336-9.
  - License: Apache version 2.0

- MolVS:
  - Version: v0.1.1 2019 release
  - Location:
    `/SMILESMerge/smilesmerge/operators/convert_files/gypsum_dl/gypsum_dl/molvs`
  - Citation: https://molvs.readthedocs.io; Take from
    https://github.com/mcs07/MolVS
  - License: MIT License

## Running SMILESMerge

There are two ways of running SMILESMerge: (1) through command-line interface; (2)
through GUI interface.

### Command-line Interface

To run SMILESMerge, use the python script `RunSMILESMerge.py`, located in the top
SMILESMerge directory, from the command line. SMILESMerge accepts user input via
two methods:

1. Command-line submission: executing directly from the command line.

```bash
cd /PATH_TO/SMILESMerge/

python RunSMILESMerge.py \
    --source_compound_file /SMILESMerge/source_compounds/naphthalene_smiles.smi \
    --root_output_folder /PATH_TO/output_directory/ \
    --number_of_crossovers 50 \
    --number_of_processors -1 \
    --LipinskiLenientFilter \
    --max_variants_per_compound 5 \
    >  /PATH_TO/OUTPUT/text_file.txt 2>  /PATH_TO/OUTPUT/text_errormessage_file.txt
```

2. json file submission: store SMILESMerge parameters in a .json file

```bash
cd /PATH_TO/SMILESMerge/
python RunSMILESMerge.py -j /PATH_TO/json_file_with_variable.json
```

Examples of the json files can be found in the folder
`/SMILESMerge/sample_sub_scripts/`.

### GUI Interface

The GUI interface requires the additional dependency of GOOEY (https://github.com/chriskiehl/Gooey). To use the GUI, please run the following command from a terminal with a modern Python environment:

`python /SMILESMerge/RunSMILESMerge_GUI.py `

If all dependencies are installed properly this will prompt you to select all required parameters.

## Understanding SMILESMerge Parameters

An explanation of every parameter can be retrieved by running:

```bash
python /SMILESMerge/RunSMILESMerge.py --help
```

Custom options such as custom filters,
are described in other parts of the tutorial.
Details for preparing source compound files are
provided directly below.

### Source Compound Files

Source compound files simply tab-delineated SMILES files (.SMI). Specify the
path using the parameter `--source_compound_file`.

Examples of source compound files can be found at
`/SMILESMerge/source_compounds/`

A detail log of how the examples files were prepared is located at
`/SMILESMerge/source_compounds/Example_source_compound_notes.txt`

An accessory script that converts a folder of PDB files to a tab-delineated
.SMI file is provided at
`/SMILESMerge/accessory_scripts/convert_directory_ligands_pdb_to_smi.py`

Details for using this accessory script are provided near the bottom of this
document, in the section
"/SMILESMerge/accessory_scripts/convert_directory_ligands_pdb_to_smi.py".

## Docker Submission

The `/SMILESMerge/docker/` directory contains the scripts to run SMILESMerge
within a docker container. These scripts are useful when using an OS that is
not compatible with SMILESMerge or its dependencies, such as Windows.

Prior to running these scripts, please install the docker software. Also, be
sure to ***also always run these scripts with sudo (linux/macOS) or
administrator privileges (Windows).***.

Running SMILESMerge via docker the first time will take a few minutes longer
because docker must install the dependencies. The same is true if docker
images have been purged.

Depending on the SMILESMerge settings, processor speed/count, etc., SMILESMerge
may complete within minutes or may take as long as multiple days. Please make
sure to use settings that are appropriate for your system. Using `nohup` may
be a useful wrapper for longer runs or when running jobs remotely (i.e., over
ssh).

More details are provided directly below and in the
`/SMILESMerge/docker/README.md` section.

### How to Setup SMILESMerge in Docker

Dockerized SMILESMerge requires the user to specify parameters via a JSON file
(not the command line).

To run the `SMILESMerge_in_docker.py` script:

Linux/MacOS:

1. Change into the `/SMILESMerge/docker/` directory in a bash terminal: `cd
   /SMILESMerge/docker/`
2. Run `SMILESMerge_in_docker.py` with `sudo` and supply a json file using the
   normal pathing of your system.
3. Execute `SMILESMerge_in_docker.py` with `sudo` privileges, providing it with a
   JSON file (MUST EXECUTE FROM `/SMILESMerge/docker/`): `sudo python
   SMILESMerge_in_docker.py -j ./examples/sample_submit_SMILESMerge_docker.json`
4. Results will appear in the output directory specified by the
   `--root_output_folder` parameter.

Windows OS:

1. Open a docker-enabled and bash-enabled terminal with administrator
   privileges.
2. Change into the `/SMILESMerge/docker/` directory in a bash enabled terminal:
   `cd /SMILESMerge/docker/`
3. Execute `SMILESMerge_in_docker.py` with `sudo` privileges, providing it with a
   JSON file (MUST EXECUTE FROM `/SMILESMerge/docker/`): `python
   SMILESMerge_in_docker.py -j ./examples/sample_submit_SMILESMerge_docker.json`
4. Results will appear in the output directory specified by the
   `--root_output_folder` parameter.

## Providing Custom Plugins

SMILESMerge was designed to be modular. This allows for the easy swapping of
code. SMILESMerge is intended to be a living codebase. If you have added good
custom code and would like to make it open source, please contact the authors
so that we can grow the user options.

Many of the SMILESMerge functions can be supplemented with custom options. These
functions include:

1. Custom Ligand Filters ***

*** Indicates that when using this feature, the code is automatically copied
into the appropriate SMILESMerge directory. This is only done once, so please
unittest the code prior to incorporating it into SMILESMerge. A print message
will indicate where the file has been copied. That file can be manually
deleted or overwritten by the user. Restart SMILESMerge after the custom files
have been automatically copied into the proper locations. After that the new
script should be integrated into SMILESMerge.

SMILESMerge ASSUMES ALL CUSTOM CODE HAS BEEN TESTED AND FUNCTIONS WITH SPECIFIED
I/O.

### 1. Custom Ligand Filters ***

This feature allows the user to incorporate custom python scripts for
filtering ligands. These filters are applied to ligands after they are created
by crossover but before Gypsum-DL conversion to 3D.

This custom code will be copied to the directory:
`/SMILESMerge/smilesmerge/operators/filter/filter_classes/filter_children_classes/`

#### Script Formatting

These filters use a class-based inheritance architecture with filter classes
that must

1. Inherit ParentFilterClass located at
   `/SMILESMerge/smilesmerge/operators/filter/filter_classes/parent_filter_class.py`
2. Have a unique name: `class unique_name(ParentFilter)` (`unique_name` cannot
   match one of the predefined filters)
3. Have at least one function called `run_filter` (`run_filter` takes a single
   variable which must be an rdkit molecule object).

#### Running Custom Filters

Because parameters can be specified via command-line or JSON file, we provide
an example of each when submitting custom filters.

1. Command-line submission:
   - Custom file is located at `/PATH_TO/custom_filter_1.py`
   - Unique class name is `custom_filter_1` (this will be how it is called in
     future versions)
   - To run multiple custom filters replace
       `[["custom_filter_1","/PATH_TO/custom_filter_1.py"]]` with:
       `[["custom_filter_1","/PATH_TO/custom_filter_1.py"],["custom_filter_2","/PATH_TO/custom_filter_2.py"]]`

```bash
python RunSMILESMerge.py \
    ... \
    --alternative_filter [["custom_filter_1","/PATH_TO/custom_filter_1.py"]]
```

2. JSON file submission:
    - Where the custom file is located at `/PATH_TO/custom_filter_1.py`
    - Unique class name is `custom_filter_1` (this will be how it is called in
      future submissions)
    - To run multiple files  Replace
      `[["custom_filter_1","/PATH_TO/custom_filter_1.py"]]` with:
      `[["custom_filter_1","/PATH_TO/custom_filter_1.py"],["custom_filter_2","/PATH_TO/custom_filter_2.py"]]`

```json
{
    ... ,
    "alternative_filter": [["custom_filter_1","/PATH_TO/custom_filter_1.py"]]
}
```

Submit in terminal: `python RunSMILESMerge.py -j
/PATH_TO/json_file_with_variable.json`

## Other Factors for Consideration Prior to Running SMILESMerge

### Processors and Multiprocessing Style

SMILESMerge can be run on a local computer such as a laptop or PC,
as well as larger high-performance clusters and servers.

#### If Running on a Laptop or PC

We recommend lowering some SMILESMerge parameters to reduce the computational
overhead for smaller machines.

- Lower the population size and number of generations. This will mean a less
  intense search of chemistry space but will make run times more reasonable.
- Lower the `max_variation` to 1. This means for every ligand created by
  SMILESMerge, we will only create 1 conformer and thus only dock once per
  ligand. This of course means a trade-off of getting more useful information
  for each ligand for computational efficiency.

We also recommend considering how long you can allow the computer to run. If
you need to continually use the computer while running SMILESMerge then you want
to fix the `number_of_processors` to leave several available to perform other
activities.

If you can leave the computer to run undisturbed for an extended period we
recommend setting `number_of_processors = -1`, which will use all available
processors.

#### If Running On A Larger Super Computer

We recommend fixing the `number_of_processors` to however many processors you
will be dedicating to SMILESMerge. If `number_of_processors = -1` than all
available processors will be use to run SMILESMerge.

#### If Running On A Cluster

We recommend setting the `number_of_processors = -1` and defining the number
of processors in an SBATCH-type submission script.

## Multiprocessing/MPI/Parallelization/Parallelizer

SMILESMerge uses the `Parallelizer.py` script from Gypsum-DL
(`/SMILESMerge/smilesmerge/operators/convert_files/gypsum_dl/gypsum_dl/Parallelizer.py`).

This script creates a Parallelizer class object which can divide jobs in three
manners:

1. Serial: run all jobs 1 at a time
2. Multiprocessing: dynamically allocated distribution of jobs across multiple
   cpus on the same device
3. MPI: static allocation of jobs across many cpus across multiple machines.

### Important Notes when Running on Clusters Using SLURM

1. Multiprocessing: When running SMILESMerge in **Multiprocessing mode** using
   SLURM, one should:
   1. 1st run the `cache_prerun` option on a single processor. `srun -n 1
      python RunSMILESMerge.py -c`
      - USE `srun` or `mpirun` for the `cache_prerun`. This limits the
        `prerun` to a single processor thus preventing errors caused by race
        conditions when creating pycache files.
   2. Then run SMILESMerge as intended. `python RunSMILESMerge.py -j
      custom_parameters.json`
      - Do not use `srun` or `mpirun` for the production run. cpu/job
        distribution is handled internally. Using `srun` or `mpirun` can cause
        errors with the `mpi4py` universe.
2. MPI: When running SMILESMerge in **MPI mode** using SLURM, one should:
    1. 1st run the `cache_prerun` option on a single processor. `srun -n 1
       python RunSMILESMerge.py -c`
       - USE `srun` or `mpirun` for the `cache_prerun`. This limits the prerun
         to a single processor thus preventing errors caused by race
         conditions when creating pycache files.
    2. Then run the simulation as intended.
        - `mpirun -n num_processors python -m mpi4py RunSMILESMerge.py -j
          custom_parameters.json`
        - Make sure to provide the `-m mpi4py` before `RunSMILESMerge.py`. This
          tells python how to handle Exceptions.

## Accessory Scripts

SMILESMerge provides several accessory scripts for preparing files, processing
data, and analyzing data.

These files can be found within the `/SMILESMerge/accessory_scripts/` folder.

### Preparation Scripts Pre-Run

#### /SMILESMerge/accessory_scripts/remove_duplicates_from_smi.sh

This script accepts a file path to a tab-delineated .smi file. It then filters
the file for redundancies in the 1st and 2nd columns of the file.

The output file is the input file + '_no_dup.smi'

This script uses Bash rather than Python because it is less memory intensive
when dealing with large .smi files in the millions-of-compounds range. This is
important when filtering through large databases such as ZINC15.

This script takes one input variable (`filename` str: Required). This is the
path to the tab-delineated .smi file to remove any redundancies.

Example submit:

```bash
bash /SMILESMerge/accessory_scripts/remove_duplicates_from_smi.sh \
    /PATH_TO/TO/SMILES.smi
```

#### /SMILESMerge/accessory_scripts/convert_directory_ligands_pdb_to_smi.py

This script converts a directory of pdb files (small molecules only, not
proteins) to SMILES and creates a single .smi file with all SMILES.

This script takes 3 input arguments:

1. `--source_folder` str (-s) Required. Path to folder containing .pdb files to
    convert. File must contain a single small molecules. Without proteins.
    Files must end with either .pdb or .PDB'
2. `--output_folder` str (-o) Required. Path to folder where we will output a
    .smi file of converted .pdb files.
3. `--number_of_processors` int (-p). Number of processors to use for parallel
    calculations. This script is not MPI enable but is able to multithread
    using SMP architecture. Set to -1 for all available CPUs.

Example run:

```bash
python /SMILESMerge/accessory_scripts/convert_directory_ligands_pdb_to_smi.py \
    --source_folder /PATH_TO/OF/PDBS/ \
    --output_folder /PATH_TO/TO/OUTPUT/ \
    --number_of_processors -1
```

#### /SMILESMerge/accessory_scripts/fragmenter_of_smi_mol.py

This script will fragment compounds from a .smi file. It is useful for lead
optimization. This script was used for the PARPi lead-optimization runs in the
AutoGrow4 paper.

This can fragment compounds in two manners:

1. BRICS decomposition: This fragments along synthesizable bonds
2. Fragment rotatable bonds: This breaks compounds along rotatable bonds.
   There is an option to skip carbon-carbon single bonds.

For each molecule, all permutation of fragments are calculated. For example,
fragment rotatable bonds `C-O-C1CCCC1` could produce the following fragments:

- `C-O-C1CCCC1`, not breaking any bonds
- `C` and `O-C1CCCC1`, breaking the 1st bond
- `C-O` and `C1CCCC1`, breaking the 2nd bond
- `C` and `O` and `C1CCCC1`, breaking the 1st bond and 2nd bond

A limit on maximum number of fragments per compound and a minimum number of
atoms per fragment can be set.

This script takes seven input arguments:

1. `--smi_file` str Required. Path to tab-delineated .smi file to fragment
2. `--output_smi_file` str (-o). Path to output tab-delineated .smi file of
   fragments. If not provided it will play a file in the same directory as
   smi_file titled smi_file + _Fragmented.smi
3. `--frags_per_seed_lig` int. Number of fragments to create per input SMILES.
   default is -1 which mean all possible fragments.
4. `--run_brics` bool. Whether to fragment ligands using BRICS fragmentation.
   This fragments along synthesizable bonds. Default is True.
5. `--run_frag` bool. Whether to fragment ligands over all rotatable bonds.
   Default is True.
6. `--c_c_bonds_off` bool. Whether to exclude fragmenting carbon-carbon single
   bonds. Default is True. If True it will ignore fragments on C-C bonds; if
   False it will fragment.
7. `--number_of_processors` int (-p). Number of processors to use for parallel
   calculations. This script is not MPI enable but is able to multithread
   using SMP architecture. Set to -1 for all available CPUs.

Example run:

```bash
python /SMILESMerge/accessory_scripts/fragmenter_of_smi_mol.py \
    -smi_file /PATH_TO/OF/SMILES.smi
```

### Graph Generation For Post-Run Analysis

#### /SMILESMerge/accessory_scripts/make_lineage_figures.py

This script creates figures that list all ligands that parented a given
ligand.

All compounds for the entire SMILESMerge run will be compiled into a dictionary
that is used to trace lineages. We pickle these dictionaries so these
dictionaries do not need to be recreated every time the script is run. For
this reason the 1st time running this script will take longer than future
runs. A pre-run option will compile these data sets without generating
figures.

1. `--output_dir` str (-o): Required. Path to folder to output files. will be
    created if does not exist
2. `--input_dir` str (-i): Required. Path to input folder containing the
    SMILESMerge run. This should be the top folder which contains the vars.json
    file.
3. `--mol_name` str: Required unless prerun. This is the name of the molecule
    whose lineage will be traced back. If not provided or None, the script
    will simply compile the necessary dictionaries/picklefiles and then
    terminate. These pickle files are stored in the input folder containing
    the vars.json file from the SMILESMerge run. Example mol_name:
    `Gen_1_Cross_539621` or `Gen_1_Cross_539621__1`. Can also be provided as
    full-name ie: `(naphthalene_45+naphthalene_90)Gen_1_Cross_539621`
4. `--source_compound_file` str: Required. This is the source .smi file used to
    seed generation zero of the SMILESMerge run. This is an essential file.
5. `--pre_run` bool. If True this will compile the necessary
    dictionaries/picklefiles and then terminate. These pickle files are stored
    in the input folder containing the vars.json file from the SMILESMerge run.

Example submit:

```bash
python /SMILESMerge/accessory_scripts/make_lineage_figures.py \
    -o /PATH_TO/output_folder/ \
    -i /PATH_TO/INPUT_RUN/Run_3/ \
    --source_compound_file /SMILESMerge/source_compounds/PARPI_BRICS_frags.smi \
    --mol_name Gen_17_Cross_727024
```
