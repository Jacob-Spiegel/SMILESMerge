"""
Populates an AutoGrow generation via mutation and crossover.
Also filters and converts SMILES to 3d SDFS.
"""
import __future__

import os
import random
import copy
import sys

import rdkit
import rdkit.Chem as Chem

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import SMILESMerge.operators.filter.execute_filters as Filter
import SMILESMerge.operators.mutation.execute_mutations as Mutation
import SMILESMerge.operators.crossover.execute_crossover as execute_crossover
import SMILESMerge.operators.convert_files.conversion_to_3d as conversion_to_3d
import SMILESMerge.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


def get_usable_format(infile):
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles



#############
# Main run Autogrow operators to make a generation
#############
def populate_generation(vars):
    """
    This will run all of the mutations, crossovers, and filters for a single
        generation. Populates a new generation of ligands.

    Inputs:
    :param dict vars: a dictionary of all user variables

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_generation_smiles_list: list with the new population
        of ligands
    :returns: bool None: returns None twice if any step failed. This will
        result in the program ending
    """
    number_of_processors = int(vars["number_of_processors"])

    # Determine which generation it is and how many mutations and crossovers
    # to make

    num_crossovers = vars["number_of_crossovers"]
    num_mutations = vars["number_of_mutants"]

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.

    # Total Population size of this generation
    total_num_desired_new_ligands = num_crossovers + num_mutations

    print("MAKE MUTATIONS")
    # Making Mutations

    # Package user vars specifying the Reaction library to use for mutation
    rxn_library_variables = [
        vars["rxn_library"],
        vars["rxn_library_file"],
        vars["function_group_library"],
        vars["complementary_mol_directory"],
    ]

    # List of SMILES from mutation
    new_mutation_smiles_list = []

    seed_list = get_complete_list_prev_gen_or_source_compounds(vars)
    # Save seed list
    save_ligand_list(
        vars["output_directory"],
        "",
        "Seed_List",
    )

    seed_list_mutations = copy.deepcopy(seed_list)

    # Make all the required ligands by mutations
    while len(new_mutation_smiles_list) < num_mutations:
        sys.stdout.flush()

        num_mutants_to_make = num_mutations - len(new_mutation_smiles_list)

        # Make all mutants
        new_mutants = Mutation.make_mutants(
            vars,
            1,
            number_of_processors,
            num_mutants_to_make,
            seed_list_mutations,
            new_mutation_smiles_list,
            rxn_library_variables,
        )
        if new_mutants is None:
            # try once more
            new_mutants = Mutation.make_mutants(
                vars,
                1,
                number_of_processors,
                num_mutants_to_make,
                seed_list_mutations,
                new_mutation_smiles_list,
                rxn_library_variables,
            )

        if new_mutants is None:
            break

        # Remove Nones:
        new_mutants = [x for x in new_mutants if x is not None]

        for i in new_mutants:
            new_mutation_smiles_list.append(i)
            if len(new_mutation_smiles_list) == num_mutations:
                break
    sys.stdout.flush()

    # save new_mutation_smiles_list
    save_ligand_list(
        vars["output_directory"],
        new_mutation_smiles_list,
        "Chosen_Mutants",
    )

    if (
            new_mutation_smiles_list is None
            or len(new_mutation_smiles_list) < num_mutations
    ):
        print("")
        print("")
        print("We needed to make {} ligands through Mutation".format(num_mutations))
        print(
            "We only made {} ligands through Mutation".format(
                len(new_mutation_smiles_list)
            )
        )
        print("")
        print("")
        raise Exception("Mutation failed to make enough new ligands.")

    print("FINISHED MAKING MUTATIONS")

    # Get starting compounds to seed Crossovers
    seed_list_crossovers = copy.deepcopy(seed_list)

    print("MAKE CROSSOVERS")
    sys.stdout.flush()

    # Making Crossovers
    # List of smiles from crossover
    new_crossover_smiles_list = []

    # Make all the required ligands by Crossover
    while len(new_crossover_smiles_list) < num_crossovers:
        sys.stdout.flush()
        num_crossovers_to_make = num_crossovers - len(new_crossover_smiles_list)

        # Make all crossovers
        new_crossovers = execute_crossover.make_crossovers(
            vars,
            1,
            number_of_processors,
            num_crossovers_to_make,
            seed_list_crossovers,
            new_crossover_smiles_list,
        )
        if new_crossovers is None:
            # try once more
            new_crossovers = execute_crossover.make_crossovers(
                vars,
                1,
                number_of_processors,
                num_crossovers_to_make,
                seed_list_crossovers,
                new_crossover_smiles_list,
            )
        if new_crossovers is None:
            break

        # Remove Nones:
        new_crossovers = [x for x in new_crossovers if x is not None]

        # append those which passed the filter
        for i in new_crossovers:
            new_crossover_smiles_list.append(i)
            if len(new_crossover_smiles_list) == num_crossovers:
                break

    # save new_crossover_smiles_list
    save_ligand_list(
        vars["output_directory"],
        new_crossover_smiles_list,
        "Chosen_Crossovers",
    )

    if (
            new_crossover_smiles_list is None
            or len(new_crossover_smiles_list) < num_crossovers
    ):
        print("")
        print("")
        print("We needed to make {} ligands through Crossover".format(num_crossovers))
        print(
            "We only made {} ligands through Crossover".format(
                len(new_crossover_smiles_list)
            )
        )
        print("")
        print("")
        raise Exception("Crossover failed to make enough new ligands.")

    print("FINISHED MAKING CROSSOVERS")
    sys.stdout.flush()


    # make a list of all the ligands from mutations, crossovers, and from the
    # last generation
    new_generation_smiles_list = []
    full_generation_smiles_list = []
    for i in new_mutation_smiles_list:
        new_generation_smiles_list.append(i)
        full_generation_smiles_list.append(i)

    for i in new_crossover_smiles_list:
        new_generation_smiles_list.append(i)
        full_generation_smiles_list.append(i)

    if len(full_generation_smiles_list) < total_num_desired_new_ligands:
        print("We needed ", total_num_desired_new_ligands)
        print("We made ", len(full_generation_smiles_list))
        print(
            "population failed to make enough mutants or crossovers... \
            Errors could include not enough diversity, too few seeds to \
            the generation, the seed mols are unable to cross-over due \
            to lack of similariy, or all of the seed lack functional groups \
            for performing reactions"
        )
        return None, None, None

    # Save the Full Generation
    smiles_to_convert_file, new_gen_folder_path = save_generation_smi(
        vars["output_directory"],
        new_generation_smiles_list,
        "New_SMILES",
    )

    sys.stdout.flush()

    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with rdkit
    # This will output sdf files into a folder. The .smi.0.sdf file is not a
    # valid mol, but all the others will be valid the 1st Smiles in the
    # original .smi file is saved as .smi.1.sdf and 2nd file is saved as
    # .smi.2.sdf
    if vars["convert_to_3D"] is True:
        conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, new_gen_folder_path)
    sys.stdout.flush()

    return smiles_to_convert_file, full_generation_smiles_list

#############
# Get seeds
#############
def test_source_smiles_convert(smile_info):
    """
    This attempts to convert a SMILES string to an rdkit.Chem.rdchem.Mol
    object
        -done in a try statement so that it is some bad SMILES string which is
        incapable of being converted
        - it also checks that the SMILES string is able to be sanitized

    Inputs:
    :param list smile_info: a list containing the SMILES of a ligand, its ID
        and potentially additional information about the ligand

    Returns:
    :returns: list smile_info: If it passed the test, it returns the list
        containing the SMILES of a ligand, its ID and potentially additional
        information about the ligand
    :returns: str printout: If it failed to convert it returns the error
        message. This passess out to prevent MPI print issues
    """
    if smile_info is None or len(smile_info) == 0:
        printout = (
            "REMOVING SMILES FROM SOURCE LIST: Blank "
            + "entry in source compound list.\n"
        )
        printout = printout + "\tRemoving: {}".format(smile_info)
        return printout

    if len(smile_info) == 1:
        printout = "REMOVING SMILES FROM SOURCE LIST: Unformatted or blank "
        printout = printout + "entry in source compound list.\n"
        printout = printout + "\tRemoving: {}".format(smile_info)
        return printout

    # separate out SMILES str and ID
    smile_str = smile_info[0]
    smile_id = str(smile_info[1])

    if type(smile_str) is not type(""):
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string is not a "
        printout = printout + "String. Check for formatting errors. \n"
        printout = printout + "\tIgnored SMILES is: {}".format(smile_str)
        return printout

    # Try importing it into RDKit with Sanitization off. Tests for errors in
    # having the wrong data type
    try:
        mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    except:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string failed "
        printout = printout + "to import into RDKit.\n\t "
        printout = printout + "Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # This will fail if there are valence errors. We won't try to correct
    # someones source compound list Although the MOH.check_sanitization will
    # do that. try sanitizing, which is necessary later
    try:
        Chem.SanitizeMol(mol)
    except:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES "
        printout = printout + "string failed to Sanitize in RDKit.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Make the mol again fresh and try running it through MOH.handleHs() This
    # will try protanating and Deprotanating the mol. If it can't handle that
    # We reject it as many functions will require this sort of manipulation.
    # More advanced sanitization issues will also be removed in this step
    mol = Chem.MolFromSmiles(str(smile_str), sanitize=False)
    mol = MOH.handleHs(mol, True)

    if mol is None:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string failed \
                    to be protanated or deprotanated.\n"
        printout = (
            printout
            + "\t This is often an issue with valence and sanitization "
            + "issues with the SMILES string."
        )
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Check there are no * which are atoms with atomic number=0
    mol = MOH.check_for_unassigned_atom(mol)
    if mol is None:
        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string contained "
        printout = printout + "an unassigned atom type labeled as *.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # Check for fragments.
    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) != 1:

        printout = "REMOVING SMILES FROM SOURCE LIST: SMILES string was fragmented.\n"
        printout = printout + "\t Removed SMILE string is: {} \n".format(smile_str)
        printout = printout + "\t Removed SMILE ID is: {}".format(smile_id)
        return printout

    # the ligand is good enough to use throughout the program!
    return smile_info


def get_complete_list_prev_gen_or_source_compounds(vars):
    """
    Get the source compounds list from the source
    compound list

    This also filters the list to ensure mols can be imported to RDKit and
    that they pass the drug-likliness filters.

    Uses the User specified starting compounds
    This takes an .smi file

    Inputs:
    :param dict vars: a dictionary of all user variables

    Returns:
    :returns: list usable_list_of_smiles: a list with SMILES strings, names,
        and information about the smiles from the source compound list
    """
    # This will be the full length list of starting molecules as the seed
    source_file = str(vars["source_compound_file"])
    usable_list_of_smiles = get_usable_format(source_file)

    if len(usable_list_of_smiles) == 0:
        print(
            "\nThere were no available ligands in source compound. Check formatting\n"
        )
        raise Exception(
            "There were no available ligands in source compound. Check formatting"
        )

    # Test that every SMILES in the usable_list_of_smiles is a valid SMILES
    # which will import and Sanitize in RDKit. SMILES will be excluded if they
    # are fragmented, contain atoms with no atomic number (*), or do not
    # sanitize
    job_input = tuple([tuple([i]) for i in usable_list_of_smiles])

    usable_list_of_smiles = vars["parallelizer"].run(
        job_input, test_source_smiles_convert
    )
    usable_list_of_smiles = [x for x in usable_list_of_smiles if x is not None]
    print_errors = [x for x in usable_list_of_smiles if type(x) is str]
    usable_list_of_smiles = [x for x in usable_list_of_smiles if type(x) is list]
    for x in print_errors:
        print(x)

    if len(usable_list_of_smiles) == 0:
        printout = "\nThere were no ligands in source compound which could sanitize.\n"
        print(printout)
        raise Exception(printout)

    if vars["filter_source_compounds"] is True:

        prefilter_list = copy.deepcopy(usable_list_of_smiles)
        print("")
        print("Running Filter on the Source Compounds")
        usable_list_of_smiles = Filter.run_filter(vars, usable_list_of_smiles)

        # Remove Nones:
        usable_list_of_smiles = [x for x in usable_list_of_smiles if x is not None]

        if len(usable_list_of_smiles) == 0:
            printout = "\nThere were no ligands in source compound which \
                        passed the User-selected Filters.\n"
            print(printout)
            raise Exception(printout)

        for lig in usable_list_of_smiles:
            failed_filter_list = []
            if lig not in prefilter_list:
                failed_filter_list.append(lig[1])

        if len(failed_filter_list) != 0:
            printout = "\n THE FOLLOWING LIGANDS WERE REMOVED FROM THE\
                        SOURCE LIST: Failed the User-selected Filters\n"
            printout = printout + "\t{}".format(failed_filter_list)
            print(printout)

    random.shuffle(usable_list_of_smiles)

    return usable_list_of_smiles


def make_seed_list(vars, source_compounds_list):
    """
    Get the starting compound list for running the Mutation and Crossovers


    Inputs:
    :param dict vars: a dictionary of all user variables
    :param list source_compounds_list: a list with SMILES strings, names, and
        information about the smiles from the source compound list
    :param int num_seed_diversity: the number of seed molecules which come
        from diversity selection

    Returns:
    :returns: list usable_list_of_smiles: a list with SMILES strings, names,
        and information about the smiles
    """
    usable_list_of_smiles = copy.deepcopy(source_compounds_list)

    # This will be the full length list of starting molecules as the seed
    random.shuffle(usable_list_of_smiles)

    return usable_list_of_smiles

#############
# Saving Output files for generations and seeds
#############
def save_generation_smi(output_directory,
                        formatted_smile_list, nomenclature_tag):
    """"
    This function saves a list of newly generated population of ligands as an
    .smi file. .smi file column 1 is the SMILES string and column 2 is its
    smile ID


    Inputs:
    :param dict output_directory: the directory of the run to save the
        generation
    :param list formatted_smile_list: list of the newly generated population
        of ligands
    :param str nomenclature_tag: The str describing the ligand list if None
        than don't add tag. It is the full list of all ligands for the
        generation.

    Returns:
    :returns: str output_file_name: name of the output file
    :returns: str output_directory: the path to the folder containing files
    """
    # folder for this new generation
    output_directory = output_directory + os.sep
    # make the name for the new file
    output_file_name = output_directory + "{}.smi".format(
        nomenclature_tag
    )

    # write as a tab delineated .smi file
    with open(output_file_name, "w") as f:
        for smile in formatted_smile_list:
            # smile_string = smile[0]
            # smile_id = smile[1]
            x = str(smile[0] + "\t" + str(smile[1]) + "\n")
            f.write(x)

    sys.stdout.flush()
    return output_file_name, output_directory


def save_ligand_list(output_directory, list_of_chosen_ligands, nomenclature_tag):
    """
    Save the list of ligands. nomenclature_tag is a string such as "Mutation"
    or "Crossover" describing what this data is used
    for.

    Inputs:
    :param dict output_directory: the directory of the run to save the
    :param list list_of_chosen_ligands: The formatted list of ligands to seed
    :param str nomenclature_tag: The str describing the ligand list
        -ie seeding_mutations is the list that seeded the mutations while
            mutations would be the list of mutations generated from the
            seeding_mutations list
        -ie. mutation, crossover
    """

    # make a folder for the Seed files
    seed_folder_path = output_directory + os.sep + "SeedFolder" + os.sep

    # check if folders exist, if not make them
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    if not os.path.isdir(seed_folder_path):
        os.makedirs(seed_folder_path)

    output_file_name = "{}{}.smi".format(
        seed_folder_path, nomenclature_tag
    )

    # save to a new output smiles file. ie. save to ranked_smiles_file
    with open(output_file_name, "w") as output:
        for line in list_of_chosen_ligands:
            output_line = "\t".join(line) + "\n"
            output.write(output_line)

    sys.stdout.flush()
