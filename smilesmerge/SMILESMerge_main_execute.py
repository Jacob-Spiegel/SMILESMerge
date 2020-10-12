"""
Top level for running SMILESMerge.
Runs all population generation (operations).
Runs plotting at end.
"""
import __future__

import os
import glob
import sys
import shutil

import smilesmerge.operators.operations as operations

def main_execute(vars):
    """
    This function takes the user variables and runs Autogrow

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """

    # Unpack necessary variables
    # output_directory is the root output folder for the run
    output_directory = vars["output_directory"]

    # This will run operations which will:
    # 1)  generate new ligands
    # 2) optionally filter ligands
    # 3) optionally convert from 1D smiles to 3D (mol2/PDB)
    
    sys.stdout.flush()


    smile_file_new_gen, new_gen_ligands_list = operations.populate_generation(vars)
    sys.stdout.flush()

    if new_gen_ligands_list is None:
        raise ValueError("Population failed to make enough crossovers... \
                            Errors could include not enough diversity, too few seeds to the generation, \
                            setting number_of_crossovers too high, or the seed mols \
                            are unable to cross-over due to lack of similarity.")

    sys.stdout.flush()
#


def find_last_generation(folder_path_string_no_gen):
    """
    This will take a folder path which is missing an interger at the end, and
    find if there are any folders which exist with that path with an interger.
    If there are it will return the highest interger that when added to the
    path exists as a directory.

    If no directories exist with that path+0 then we return None. This causes
    the starting generation of this attempt to run to be generation_0.
    Starting fresh.

    folder_path_string_no_gen = output_directory + "generation_"

    Inputs:
    :param str folder_path_string_no_gen: the folder to check.

    Returns:
    :returns: int last_gen_number: the int of the last generation number or
        None if no previous runs.
    """

    path_exists = True
    i = 1
    while path_exists is True:
        folder_path = "{}{}{}".format(folder_path_string_no_gen, i, os.sep)
        if os.path.exists(folder_path):
            i = i + 1

        else:
            path_exists = False

    if i == 1:
        # Check to see if there's a Run 0 based on the seed.
        i = 0
        folder_path = "{}{}{}".format(folder_path_string_no_gen, i, os.sep)
        if os.path.exists(folder_path) is False:
            return None

        # There are no previous runs in this directory
        last_gen_number = None

    else:
        last_gen_number = i - 1

    return last_gen_number
#

def determine_if_gen_completed(gen_dir_path, gen_number):
    """
    Check if this generation has completed or if it failed. Every generation
    which completes has a .smi file title generation_0_ranked.smi (with the
    number of the generation between the word generation and ranked).
    -If a Run failed due to either a hard crash or a soft crash, there should
        not be a ranked .smi file.

    Inputs:
    :param str gen_dir_path: is the path of the generation folder within a Run
        folder.
    :param int gen_number: The generation number of the folder.

    Returns:
    :returns: bool os.path.isfile(file_path): Returns True if the gen_dir_path
        has a ranked.smi file. Returns False if the gen_dir_path does not have a
        ranked.smi file
    """

    ranked_file_name = "generation_{}_ranked.smi".format(gen_number)
    file_path = "{}{}{}".format(gen_dir_path, os.sep, ranked_file_name)

    return os.path.isfile(file_path)
#

def delete_temporary_files_and_folders(file_or_folder):
    """
    This deletes all temporary files.

    Inputs:
    :param str file_or_folder: the file or folder to delete

    """
    if os.path.exists(file_or_folder) is True:
        if os.path.isdir(file_or_folder) is True:
            try:
                shutil.rmtree(file_or_folder)
            except:
                pass
        else:
            try:
                os.remove(file_or_folder)
            except:
                pass

        # If it failed to delete try via bash command
        if os.path.exists(file_or_folder) is True:
            command = "rm -rf {}".format(file_or_folder)
            try:
                os.system(command)
            except:
                pass
    else:
        pass
#
