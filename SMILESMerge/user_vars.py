"""user_vars
This should contain the functions for defining input variables.
Both the default variables and the user input variables.
This should also validate them.
"""

import __future__

import os
import copy
import datetime
import json
import sys
import platform
from shutil import copyfile


def program_info():
    """
    Get the program version number, etc.

    Returns:
    :returns: str program_output: a string for the print of the program information
    """
    program_output = "\nAutoGrow Version 4.0.2\n"
    program_output = program_output + " ================== \n"
    program_output = (
        program_output
        + "If you use AutoGrow 4.0.2 in your research, please cite the following reference:\n"
    )
    program_output = program_output + "Spiegel, J.O., Durrant, J.D. \n"
    program_output = program_output + "AutoGrow4: an open-source genetic algorithm "
    program_output = program_output + "for de novo drug design and lead optimization. \n"
    program_output = program_output + "J Cheminform 12, 25 (2020). \n"
    program_output = program_output + "[doi: 10.1186/s13321-020-00429-4]\n"
    program_output = program_output + " ================== \n\n"

    return program_output


#
def save_vars_as_json(vars):
    """
    This function saves the vars dictionary as a json file. This can be used
    later to track experiments and is necessary for several of the utility
    scripts.
    It saves all variables except the parallelizer class object.

    It saves the file to the output_directory + "vars.json"
        -If AutoGrow has been run multiple times for the same directory it
        will save the new vars file as append a number to the file name
        starting with 2. The util scripts will only look at the original "vars.json"
            ie) output_directory + "vars_2.json"

    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    """
    output_directory = vars["output_directory"]

    vars_file = output_directory + os.sep + "vars.json"
    if os.path.exists(vars_file):
        # vars.json already exist
        # lets make the next file
        path_exists = True
        i = 2
        while path_exists is True:
            vars_file = "{}{}vars_{}.json".format(output_directory, os.sep, i)
            if os.path.exists(vars_file):
                i = i + 1
            else:
                path_exists = False

    temp_vars = {}
    for k in vars.keys():
        if "parallelizer" in k or k == "filter_object_dict":
            continue

        temp_vars[k] = copy.deepcopy(vars[k])

    with open(vars_file, "w") as fp:
        json.dump(temp_vars, fp, indent=4)


def multiprocess_handling(vars):
    """
    This function handles the multiprocessing functions. It establishes a Paralellizer object
    and adds it to the vars dictionary.
    Inputs:
    :param dict vars: dict of user variables which will govern how the programs runs
    Returns:
    :returns: dict vars: dict of user variables which will govern how the programs runs
    """

    # Handle Serial overriding number_of_processors
    # serial fixes it to 1 processor
    if vars["multithread_mode"].lower() == "serial":
        vars["multithread_mode"] = "serial"
        if vars["number_of_processors"] != 1:
            print(
                "Because --multithread_mode was set to serial, "
                + "this will be run on a single processor."
            )
        vars["number_of_processors"] = 1

    # Handle mpi errors if mpi4py isn't installed
    if vars["multithread_mode"].lower() == "mpi":
        vars["multithread_mode"] = "mpi"
        try:
            import mpi4py
        except:
            printout = "mpi4py not installed but --multithread_mode is set to"
            printout = printout + " mpi. \n Either install mpi4py or switch "
            printout = printout + "multithread_mode to multithreading or serial"
            raise ImportError(printout)

        try:
            import func_timeout
            from func_timeout import func_timeout, FunctionTimedOut
        except:
            printout = "func_timeout not installed but --multithread_mode is "
            printout = printout + "set to mpi. \n Either install func_timeout "
            printout = printout + "or switch multithread_mode to"
            printout = printout + " multithreading or serial"
            raise ImportError(printout)

    # # # launch mpi workers
    if vars["multithread_mode"] == "mpi":
        # Avoid EOF error
        from SMILESMerge.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import (
            Parallelizer,
        )

        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"]
        )

        if vars["parallelizer"] is None:
            printout = "EOF ERRORS FAILED TO CREATE A PARALLIZER OBJECT"
            print(printout)
            raise Exception(printout)

    else:
        # Lower level mpi (ie making a new Parallelizer within an mpi)
        #   has problems with importing the MPI environment and mpi4py
        #   So we will flag it to skip the MPI mode and just go to multithread/serial
        # This is a saftey precaution
        from SMILESMerge.operators.convert_files.gypsum_dl.gypsum_dl.Parallelizer import Parallelizer

        vars["parallelizer"] = Parallelizer(
            vars["multithread_mode"], vars["number_of_processors"], True
        )

    return vars


############################################
###### Variables Handlining Settings #######
############################################
def check_for_required_inputs(input_params):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict input_params: The parameters. A dictionary of {parameter name: value}.
    """
    keys_from_input = list(input_params.keys())

    list_of_required_inputs = [
        "root_output_folder",
        "source_compound_file"
    ]

    missing_variables = []
    for variable in list_of_required_inputs:
        if variable in keys_from_input:
            continue
        missing_variables.append(variable)

    if len(missing_variables) != 0:
        printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunAutogrow -h"
        printout = printout + "\nThe following required variables are missing: "
        for variable in missing_variables:
            printout = printout + "\n\t" + variable
        print("")
        print(printout)
        print("")
        raise NotImplementedError("\n" + printout + "\n")


    #######################################
    # Check that all required files exist #
    #######################################

    # convert paths to abspath, in case necessary
    input_params["root_output_folder"] = os.path.abspath(
        input_params["root_output_folder"]
    )
    input_params["source_compound_file"] = os.path.abspath(
        input_params["source_compound_file"]
    )

    # Check root_output_folder exists
    if os.path.exists(input_params["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(input_params["root_output_folder"])
        except:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

        if os.path.exists(input_params["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(input_params["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(input_params["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file can not be found. \
            File must be a tab delineated .smi file."
        )
    if ".smi" not in input_params["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )



def determine_bash_timeout_vs_gtimeout():
    """
    This function tests whether we should use the BASH command "timeout" (for linux)
     or the coreutils function "gtimeout" for MacOS which can be obtained
     through homebrew

    Returns:
    :returns: str timeout_option: A string either "timeout" or "gtimeout" describing
     whether the bash terminal is able to use the bash function timeout or gtimeout
    """

    if sys.platform.lower() in ["linux", "linux2"]:
        # Should be true and default installed in all Linux machines
        return "timeout"

    command = 'timeout 1 echo " "'
    # Running the os.system command for command will return 0,1, or 32512
    # 0 means that the timeout function works (most likely this is a linux os)
    # 32512 means that the timeout function DOES NOT Work (most likely this is MacOS)

    try:  # timeout or gtimeout
        timeout_result = os.system("g" + command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported \
            by Autogrow or you may need to execute through Bash."
        )
    if timeout_result == 0:
        timeout_option = "gtimeout"
        return timeout_option
    print("gtimeout failed to run, we will check timeout")

    try:  # timeout or gtimeout
        timeout_result = os.system(command)
    except:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
            Autogrow or you may need to execute through Bash."
        )

    if timeout_result == 0:
        timeout_option = "timeout"
        return timeout_option

    printout = "Need to install GNU tools for Bash to work. \n"
    printout = (
        printout
        + "This is essential to use Bash Timeout function in SMILESMerge. \n"
    )
    printout = printout + "\t This will require 1st installing homebrew. \n"
    printout = printout + "\t\t Instructions found at: https://brew.sh/ \n"
    printout = printout + "\t Once brew is installed, please run:"
    printout = printout + " sudo brew install coreutils \n\n"
    print(printout)
    raise Exception(printout)


def check_dependencies():
    """
    This function will try to import all the installed dependencies that will be
    used in SMILESMerge. If it fails to import it will raise an ImportError
    """

    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option not in ["timeout", "gtimeout"]:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
        Autogrow or you may need to execute through Bash."
        )

    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem.Draw import PrepareMolForDrawing
        from rdkit.Chem import rdFMCS
        from rdkit.Chem import FilterCatalog
        from rdkit.Chem.FilterCatalog import FilterCatalogParams
        import rdkit.Chem.Lipinski as Lipinski
        import rdkit.Chem.Crippen as Crippen
        import rdkit.Chem.Descriptors as Descriptors
        import rdkit.Chem.MolSurf as MolSurf

    except:
        print("You need to install rdkit and its dependencies.")
        raise ImportError("You need to install rdkit and its dependencies.")

    # molvs is prepackaged within gypsum_dl
    # try:
    #     from molvs import standardize_smiles as ssmiles
    # except:
    #     print("You need to install molvs and its dependencies.")
    #     raise ImportError("You need to install molvs and its dependencies.")

    try:
        import numpy
    except:
        print("You need to install numpy and its dependencies.")
        raise ImportError("You need to install numpy and its dependencies.")

    try:
        from scipy.cluster.vq import kmeans2
    except:
        print("You need to install scipy and its dependencies.")
        raise ImportError("You need to install scipy and its dependencies.")

    try:
        import os
        import sys
        import glob
        import subprocess
        import multiprocessing
        import time
    except:
        print(
            "Missing a Python Dependency. Could be import: os,sys,glob,\
                subprocess,multiprocess, time."
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
                os,sys,glob,subprocess,multiprocess, time."
        )

    try:
        import copy
        import random
        import string
        import math
    except:
        print("Missing a Python Dependency. Could be import: copy,random, string,math")
        raise ImportError(
            "Missing a Python Dependency. Could be import: copy,random, string,math"
        )

    try:
        from collections import OrderedDict
        import webbrowser
        import argparse
        import itertools
        import unittest
    except:
        print(
            "Missing a Python Dependency. Could be import: collections,\
            webbrowser,argparse,itertools,unittest"
        )
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            collections,webbrowser,argparse,itertools,unittest"
        )

    try:
        import textwrap
        import pickle
        import json
    except:
        print("Missing a Python Dependency. Could be import: textwrap, pickle,json")
        raise ImportError(
            "Missing a Python Dependency. Could be import: \
            textwrap, pickle,json"
        )


def define_defaults():
    """
    Sets the command-line parameters to their default values.

    Returns:
    :returns: dict vars: a dictionary of all default variables
    """

    vars = {}

    # where we are currently (absolute filepath from route)
    # used for relative pathings
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Crossover function
    vars["max_time_mcs_prescreen"] = 1
    vars["max_time_mcs_thorough"] = 1
    vars["min_atom_match_mcs"] = 4
    vars["protanate_step"] = False

    # processors
    vars["number_of_processors"] = 1
    vars["multithread_mode"] = "multithreading"

    # Populations settings
    vars["filter_source_compounds"] = True
    vars["num_generations"] = 1
    vars["number_of_crossovers"] = 10

    # Filters
    vars["LipinskiStrictFilter"] = False
    vars["LipinskiLenientFilter"] = False
    vars["GhoseFilter"] = False
    vars["GhoseModifiedFilter"] = False
    vars["MozziconacciFilter"] = False
    vars["VandeWaterbeemdFilter"] = False
    vars["PAINSFilter"] = False
    vars["NIHFilter"] = False
    vars["BRENKFilter"] = False
    vars["No_Filters"] = False
    vars["alternative_filter"] = None

    # gypsum # max variance is the number of conformers made per ligand
    vars["convert_to_3D"] = True
    vars["max_variants_per_compound"] = 3
    vars["gypsum_thoroughness"] = 3
    vars["min_ph"] = 6.4
    vars["max_ph"] = 8.4
    vars["pka_precision"] = 1.0
    vars["gypsum_timeout_limit"] = 10

    # Other vars
    vars["debug_mode"] = False
    vars["reduce_files_sizes"] = False
    
    # Check Bash Timeout function (There's a difference between MacOS and linux)
    # Linux uses timeout while MacOS uses gtimeout
    timeout_option = determine_bash_timeout_vs_gtimeout()
    if timeout_option in  ["timeout", "gtimeout"]:
        vars["timeout_vs_gtimeout"] = timeout_option
    else:
        raise Exception(
            "Something is very wrong. This OS may not be supported by \
             Autogrow or you may need to execute through Bash."
        )

    return vars

############################################
######## Input Handlining Settings #########
############################################
def convert_json_params_from_unicode(params_unicode):
    """
    Set the parameters that will control this ConfGenerator object.

    :param dict params_unicode: The parameters. A dictionary of {parameter name:
                value}.
    Returns:
    :returns: dict params: Dictionary of User variables
    """
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    if sys.version_info < (3,):
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, unicode):
                val = str(val).encode("utf8")
            key = param.encode("utf8")
            params[key] = val
    else:
        for param in params_unicode:
            val = params_unicode[param]
            key = param
            params[key] = val
    return params


def check_value_types(vars, argv):
    """
    This checks that all the user variables loaded in use that same or comparable
    datatypes as the defaults in vars. This prevents type issues later in the
    simulation.

    Given the many uservars and the possibility for intentional differences,
    especially as the program is developed, this function tries to be
    NOT OPINIONATED, only correcting for several obvious and easy to correct issues
    of type discrepancies occur between argv[key] and vars[key]
        ie
            1) argv[key] = "true" and vars[key] = False
                this script will not change argv[key] to False... it will
                convert "true" to True
                ---> argv[key]=True
            2) argv[key] = "1.01" and vars[key] = 2.1
                this script will change argv[key] from "1.01" to float(1.01)

    Inputs:
    :param dict vars: Dictionary of program defaults, which will later be
        overwritten by argv values
    :param dict argv: Dictionary of User specified variables
    Returns:
    :returns: dict vars: Dictionary of program defaults, which will later
        be overwritten by argv values
    :returns: dict argv: Dictionary of User specified variables
    """
    for key in list(argv.keys()):
        if key not in list(vars.keys()):
            # Examples may be things like root_output_folder
            #   Just skip these
            continue

        if type(argv[key]) != type(vars[key]):
            # Several variable default is None which means checks are
            # processed elsewhere...
            if vars[key] is None:
                # check argv[key] is "none" or "None"
                if type(argv[key]) == str:
                    if argv[key].lower() == "none":
                        argv[key] = None
                else:
                    continue

            # Handle number types
            elif type(vars[key]) == int or type(vars[key]) == float:
                if type(argv[key]) == int or type(argv[key]) == float:
                    # this is fine
                    continue
                elif type(argv[key]) == str:
                    try:
                        temp_item = float(argv[key])
                        if type(temp_item) == float:
                            argv[key] = temp_item
                        else:
                            printout = "This parameter is the wrong type.\n \t Check : "
                            printout = printout + "{} type={}\n".format(
                                key, type(argv[key])
                            )
                            printout = printout + "\t Should be type={}\n\t".format(
                                type(vars[key])
                            )
                            printout = (
                                printout
                                + "Please check Autogrow documentation using -h"
                            )
                            raise IOError(printout)
                    except:
                        printout = "This parameter is the wrong type. \n \t Check :"
                        printout = printout + " {} type={}\n".format(
                            key, type(argv[key])
                        )
                        printout = printout + "\t Should be type={}\n\t".format(
                            type(vars[key])
                        )
                        printout = (
                            printout + "Please check Autogrow documentation using -h"
                        )
                        raise IOError(printout)
                else:
                    printout = "This parameter is the wrong type. \n \t Check :"
                    printout = printout + " {} type={}\n".format(key, type(argv[key]))
                    printout = printout + "\t Should be type={}\n\t".format(
                        type(vars[key])
                    )
                    printout = printout + "Please check Autogrow documentation using -h"
                    raise IOError(printout)
            elif type(vars[key]) == bool:
                if argv[key] is None:
                    # Do not try to handle this. May make sense.
                    continue
                if type(argv[key]) == str:
                    if argv[key].lower() in ["true", "1"]:
                        argv[key] = True
                    elif argv[key].lower() in ["false", "0"]:
                        argv[key] = False
                    elif argv[key].lower() in ["none"]:
                        argv[key] = None
                    else:
                        printout = "This parameter is the wrong type. \n \t Check :"
                        printout = printout + " {} type={}\n".format(
                            key, type(argv[key])
                        )
                        printout = printout + "\t Should be type={}\n\t".format(
                            type(vars[key])
                        )
                        printout = (
                            printout + "Please check Autogrow documentation using -h"
                        )
                        raise IOError(printout)
                else:
                    printout = "This parameter is the wrong type. \n \t Check :"
                    printout = printout + " {} type={}\n".format(key, type(argv[key]))
                    printout = printout + "\t Should be type={}\n\t".format(
                        type(vars[key])
                    )
                    printout = printout + "Please check Autogrow documentation using -h"
                    raise IOError(printout)
    return vars, argv


def load_in_commandline_parameters(argv):
    """
    Load in the command-line parameters

    Inputs:
    :param dict argv: Dictionary of User specified variables

    Returns:
    :returns: dict vars: Dictionary of User variables
    :returns: str printout: a string to be printed to screen and saved to output file
    """

    vars = define_defaults()

    # Load the parameters from the json
    if "json" in argv:
        json_vars = json.load(open(argv["json"]))
        json_vars = convert_json_params_from_unicode(json_vars)
        check_for_required_inputs(json_vars)
        vars, json_vars = check_value_types(vars, json_vars)
        for key in list(json_vars.keys()):
            vars[key] = json_vars[key]

    else:
        check_for_required_inputs(argv)
        argv = handle_custom_inputs_if_argparsed(argv)
        vars, argv = check_value_types(vars, argv)
        for key in list(argv.keys()):
            vars[key] = argv[key]

    vars = multiprocess_handling(vars)

    printout = "(RE)STARTING AUTOGROW 4.0: " + str(datetime.datetime.now())
    printout = printout + program_info()
    printout = (
        printout + "\nUse the -h tag to get detailed help regarding program usage.\n"
    )
    print(printout)
    sys.stdout.flush()
    # Check all Dependencies are installed
    check_dependencies()

    vars = filter_choice_handling(vars)

    ###########################################
    ########## Check variables Exist ##########
    ###########################################

    # Check if the Operating System is Windows, if so turn off Multiprocessing.
    if os.name == "nt" or os.name == "ce":
        # so it's running under windows. multiprocessing disabled
        vars["number_of_processors"] = 1
        printout = (
            printout
            + "\nWARNING: Multiprocessing is disabled on\
            windows machines.\n"
        )

    # More Handling for Windows OS
    # convert path names with spaces if this is windows
    if os.name == "nt" or os.name == "ce":
        # so it's running under windows. multiprocessing disabled
        if " " in vars["root_output_folder"]:
            vars["root_output_folder"] = '"' + vars["root_output_folder"] + '"'

    # output the paramters used
    printout = printout + "\nPARAMETERS" + "\n"
    printout = printout + " ========== " + "\n"


    # Check if the user wants to continue a run or start a new run.
    # Make new run directory if necessary. return the Run folder path
    # The run folder path will be where we place our generations and output files
    vars["output_directory"] = set_run_directory(vars["root_output_folder"])

    # Save variables in vars dict to a .json file for later usage and reference
    # It saves the file to the output_directory + "vars.json"
    # -If AutoGrow has been run multiple times for the same directory it
    # will save the new vars file as append a number to the file name
    # starting with 2. The util scripts will only look at the original "vars.json"
    #     ie) output_directory + "vars_2.json"
    save_vars_as_json(vars)

    return vars, printout


############################################
######### File Handlining Settings #########
############################################
def find_previous_runs(folder_name_path):
    """
    This will check if there are any previous runs in the output directory.
        - If there are it will return the interger of the number label of the last Run folder path.
            - ie if there are folders Run_0, Run_1, Run_2 the function will return int(2)
        - If there are no previous Run folders it returns None.

    Inputs:
    :param str folder_name_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files

    Returns:
    :returns: int last_run_number: the int of the last run number or None if no previous runs.
    """

    path_exists = True
    i = 0
    while path_exists is True:
        folder_path = "{}{}{}".format(folder_name_path, i, os.sep)
        if os.path.exists(folder_path):
            i = i + 1
        else:
            path_exists = False

    if i == 0:
        # There are no previous runs in this directory
        last_run_number = None
        return None

    # A previous run exists. The number of the last run.
    last_run_number = i - 1
    return last_run_number


def set_run_directory(root_folder_path):
    """
    Determine and make the folder for the run directory.
        -If no previous runs exist in the root_folder_path then make a new
            folder named root_folder_path + "Run_0"
        -If there are previous runs in the root_folder_path then make a
            new folder incremental increasing the name by 1 from the last
            run in the same output directory.
    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files
    Returns:
    :returns: str folder_path: the string of the newly created directory for
        puting output folders
    """

    folder_name_path = root_folder_path + os.sep + "Run_"
    print(folder_name_path)

    last_run_number = find_previous_runs(folder_name_path)

    if last_run_number is None:
        # There are no previous simulation runs in this directory
        print("There are no previous runs in this directory.")
        print("Starting a new run named Run_0.")

        # make a folder for the new generation
        run_number = 0
        folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
        os.makedirs(folder_path)

    else:
        run_number = last_run_number + 1
        folder_path = "{}{}{}".format(folder_name_path, run_number, os.sep)
        os.makedirs(folder_path)

    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_path)
    print("")
    return folder_path


############################################
########   Custom Option Settings   ########
############################################
def handle_custom_inputs_if_argparsed(input_params):
    """
    There are several Custom options such as filters. Because Filters can use multiple options
    at once it takes a list of list information.
    This function is used to properly import and parse those user variables if
     using the commandline argparse

    This function will handle those if there are used and return
    the modified input_params dict

    Inputs:
    :param dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    Returns:
    :returns: dict input_params: The parameters. A dictionary of
        {parameter name: value}.
    """

    # Custom Filters
    if "alternative_filter" not in input_params.keys():
        input_params["alternative_filter"] = None
    if (
            input_params["alternative_filter"] is not None
            and input_params["alternative_filter"] != []
    ):
        orginal = input_params["alternative_filter"][0]
        orginal = orginal.replace("[[", "[").replace("]]", "]")
        new_alternative_filter = []
        for custom_filter in orginal.split("]"):
            custom_filter = custom_filter.replace("[", "").replace("]", "")
            custom_filter = [x for x in custom_filter.split(",") if x != ""]
            if len(custom_filter) == 2:
                new_alternative_filter.append(custom_filter)
        input_params["alternative_filter"] = new_alternative_filter

    return input_params


#
def handle_alternative_filters(vars, filter_list):
    """
    This will handle Custom Filters

    Inputs:
    :param dict vars: Dictionary of User variables
    :param list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user vars["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.

    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow the same
        nomenclature and enter that filter in the user vars["alternative_filters"]
        as the name of that class and place that file in the same folder as the
        other filter classes.
    """
    if vars["alternative_filter"] is not None:
        if type(vars["alternative_filter"]) != list:
            raise Exception(
                "If you want to add Custom filters to the filter \
                child classes Must be a list of lists \
                [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
            )
        if type(vars["alternative_filter"][0]) != list:
            print(vars["alternative_filter"])
            raise Exception(
                "If you want to add Custom filters to the filter \
                child classes Must be a list of lists \
                [[name_filter1, Path/to/name_filter1.py],[name_filter2, Path/to/name_filter2.py]]"
            )

        full_children_dict = make_complete_children_dict("filter")
        scripts_to_copy = []
        for custom_class in vars["alternative_filter"]:
            if custom_class[0] not in full_children_dict.keys():
                if os.path.exists(custom_class[1]) is False:
                    # Check that the path to the original script exists.
                    raise Exception(
                        "File can not be found for alternative_filter \
                        {}\n If you want to add Custom filters to the filter child \
                        classes Must be a list of lists \
                        [[name_filter1, Path/to/name_filter1.py],\
                        [name_filter2, Path/to/name_filter2.py]]".format(custom_class[1])
                    )

                new_file = os.sep.join(
                    [
                        os.path.abspath(os.path.dirname(__file__)),
                        "operators",
                        "filter",
                        "filter_classes",
                        "filter_children_classes",
                        os.path.basename(custom_class[0]) + ".py",
                    ]
                )

                if os.path.exists(new_file) is True:
                    # File has been copied to proper dir but is not being found by the code
                    printout = "A copy of the custom script {} has been moved \
                        to {}\n".format(custom_class[1], new_file)
                    printout = (
                        printout
                        + "Unfortunately this could not be  \
                        imported by the filter module."
                    )
                    printout = (
                        printout
                        + "Please check the file naming \
                        corresponding to: {}\n\n".format(
                            custom_class
                        )
                    )
                    print(printout)
                    raise Exception(printout)

                # Add to list of scripts to copy into the filter folder
                scripts_to_copy.append([custom_class[1], new_file])

            filter_list.append(custom_class[0])
        if len(scripts_to_copy) != 0:
            for filter_info in scripts_to_copy:
                print("copying Custom class file into the FilterClasses folder:")
                print(
                    "\t Copying : {}\n\t New file: {}\n".format(
                        custom_class[1], new_file
                    )
                )
                print(
                    "AutoGrow will need to be restarted once all custom scripts \
                    have been copied to their required location."
                )
                print(
                    "This is done once so if the script needs to be changed \
                    please either remove or replace the script within the \
                    FilterClasses folder."
                )
                print(
                    "Please ensure you unit test this code properly before \
                    incorporating.\n"
                )
                copyfile(filter_info[0], filter_info[1])

            print(
                "\n########################################"
                + "#####################################"
            )
            print("AutoGrow has incorporated the custom files into"
                  + " the filter Module.")
            print(
                " AutoGrow needs to be restarted and should now "
                + "be able to run custom scripts."
            )
            print("Please ensure you unit test this code properly before incorporating.")
            print(
                "#####################################"
                + "########################################\n"
            )
            # Technically Exit intentionally but maybe should be a raise Exception
            sys.exit(0)
    return filter_list

#
def make_complete_children_dict(purpose_of_object):
    """
    This will retrieve all the names of every child class of the parent class
    This can be  filter class
    
    Inputs:
    :param str purpose_of_object: the filter
    Returns:
    :returns: dict child_dict: Dictionary of all the class objects for Filtering
    """
    if purpose_of_object == "filter":
        import SMILESMerge.operators.filter.filter_classes.filter_children_classes
        from SMILESMerge.operators.filter.filter_classes.parent_filter_class import ParentFilter as parent_object
        from SMILESMerge.operators.filter.filter_classes.get_child_filter_class import get_all_subclasses

    children = get_all_subclasses(parent_object)
    child_dict = {}
    for child in children:
        child_object = child()
        child_name = child_object.get_name()
        child_dict[child_name] = child_object

    return child_dict

############################################
######## Filter Handlining Settings ########
############################################
def filter_choice_handling(vars):
    """
    This function handles selecting the user defined Ligand filters.

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: dict vars: Dictionary of User variables with the
        chosen_ligand_filters added
    """
    if "No_Filters" in list(vars.keys()):
        if vars["No_Filters"] is True:
            chosen_ligand_filters = None
        else:
            chosen_ligand_filters, vars = picked_filters(vars)
    else:
        chosen_ligand_filters, vars = picked_filters(vars)
    vars["chosen_ligand_filters"] = chosen_ligand_filters

    import SMILESMerge.operators.filter.execute_filters as Filter


    # get child filter class object function dictionary
    vars["filter_object_dict"] = Filter.make_run_class_dict(chosen_ligand_filters)

    return vars


#
def picked_filters(vars):
    """
    This will take the user vars and return a list of the filters
    which a molecule must pass to move into the next generation.

    Inputs:
    :param dict vars: Dictionary of User variables
    Returns:
    :returns: list filter_list: a list of the class of filter which will be used
        later to check for drug likeliness for a generation.
        If a User adds their own filter they just need to follow
        the same nomenclature and enter that filter in the user
        vars["alternative_filters"] as the name of that class and place
        that file in the same folder as the other filter classes.
    """
    filter_list = []
    vars_keys = list(vars.keys())

    if "LipinskiStrictFilter" in vars_keys:
        if vars["LipinskiStrictFilter"] is True:
            filter_list.append("LipinskiStrictFilter")
    else:
        vars["LipinskiStrictFilter"] = False

    if "LipinskiLenientFilter" in vars_keys:
        if vars["LipinskiLenientFilter"] is True:
            filter_list.append("LipinskiLenientFilter")
    else:
        vars["LipinskiLenientFilter"] = False

    if "GhoseFilter" in vars_keys:
        if vars["GhoseFilter"] is True:
            filter_list.append("GhoseFilter")
    else:
        vars["GhoseFilter"] = False

    if "GhoseModifiedFilter" in vars_keys:
        if vars["GhoseModifiedFilter"] is True:
            filter_list.append("GhoseModifiedFilter")
    else:
        vars["GhoseModifiedFilter"] = False

    if "MozziconacciFilter" in vars_keys:
        if vars["MozziconacciFilter"] is True:
            filter_list.append("MozziconacciFilter")
    else:
        vars["MozziconacciFilter"] = False

    if "VandeWaterbeemdFilter" in vars_keys:
        if vars["VandeWaterbeemdFilter"] is True:
            filter_list.append("VandeWaterbeemdFilter")
    else:
        vars["VandeWaterbeemdFilter"] = False

    if "PAINSFilter" in vars_keys:
        if vars["PAINSFilter"] is True:
            filter_list.append("PAINSFilter")
    else:
        vars["PAINSFilter"] = False

    if "NIHFilter" in vars_keys:
        if vars["NIHFilter"] is True:
            filter_list.append("NIHFilter")
    else:
        vars["NIHFilter"] = False

    if "BRENKFilter" in vars_keys:
        if vars["BRENKFilter"] is True:
            filter_list.append("BRENKFilter")
    else:
        vars["BRENKFilter"] = False

    if "alternative_filter" in vars_keys:
        filter_list = handle_alternative_filters(vars, filter_list)
    else:
        vars["alternative_filter"] = None

    # if there is no user specified ligand filters but they haven't set
    # filters to None ---> set filter to default of LipinskiLenientFilter.
    if len(filter_list) == 0:
        vars["LipinskiLenientFilter"] = True
        filter_list.append("LipinskiLenientFilter")

    return filter_list, vars
