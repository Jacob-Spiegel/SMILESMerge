#!/usr/bin/env python
"""
This script will handle running SMILESMerge in a docker container.
It handles generating the docker image and container, handling the user variables,
executing SMILESMerge, and copying the files from the container to the desired directory.

This script requires a JSON file that contains all the parameters that would be
required to run SMILESMerge on a host system (ie paths on your computer).
Necessary files, such as the source compound file, will be copied into the docker.

To run SMILESMerge from within docker. Launches docker
  image. Accepts the exact same parameters as SMILESMerge, with the following
  exceptions:
    1) User variables must be supplied in JSON format.
        - Please see documentation within the tutorial manual and an example can be found:
          -  ./examples/sample_SMILESMerge_docker_json.json

    Required variables within the JSON file:
    - `-root_output_folder`: folder path on host system that results will be copied to.
    - `-source_compound_file`: Path on host system to the tab-delineate .smi
        file that will seed generation 1.

The resulting SMILESMerge output to the desired root_output_folder.

An example JSON is provided in: ./sample_SMILESMerge_docker_json.json


To run SMILESMerge in a docker, please run the `smilesmerge_in_docker.py` script:
    Example on Linux/MacOS:
        #  cd to this directory in a bash terminal
        1) cd SMILESMerge/docker/
        # Run smilesmerge_in_docker.py with sudo and supply a json file using the
        # normal pathing of your system.
        2) `sudo python smilesmerge_in_docker.py -j ./examples/sample_SMILESMerge_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

    Example on Windows OS:
        1) open a docker enabled and bash enabled terminal with administrative privileges
        #  cd to this directory in a bash terminal
        3) cd SMILESMerge/docker/
        4)  `python smilesmerge_in_docker.py -j ./examples/sample_SMILESMerge_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

"""
import __future__

import os
import shutil
import json
import argparse
import sys

def change_permissions(file_path):
    """
    This will open the permissions for a given file.

    Inputs:
    :param str file_path: Path to a file to open permissions to.
    """
    os.chmod(file_path, 0o777)


def change_permissions_recursively(file_or_folder_path):
    """
    This will open the permissions for a given file/folder.

    Skip permissions change if Windows.

    Inputs:
    :param str file_or_folder_path: Path to a file/folder to open permissions to.
    """
    if os.name == "nt" or os.name == "ce":
        # chmod and os.chmod do not apply to Windows OS so lets skip this.
        pass
    elif sys.platform.lower() in ["linux", "linux2"]:
        # chmod -R recursively open the permissions
        os.system("chmod -R a+rwx {}".format(file_or_folder_path))
    elif sys.platform.lower() == "darwin":
        # chmod -R recursively open the permissions
        os.system("chmod -R a+rwx {}".format(file_or_folder_path))
    else:
        # chmod may not be a valid command on other OS systems.
        #  So let's do this the manual way.
        if os.path.isdir(file_or_folder_path):
            directory_path_list = []
            file_list = []
            for top_dir, dir_list, list_of_files in os.walk(
                    file_or_folder_path, topdown=False
            ):
                for directory in [os.path.join(top_dir, d) for d in dir_list]:
                    directory_path_list.append(directory)
                for file_path in [os.path.join(top_dir, fil) for fil in list_of_files]:
                    file_list.append(file_path)

            # Convert mods on all files within a directory
            file_list = list(set(file_list))
            directory_path_list = list(set(directory_path_list))
            for file_path in file_list:
                change_permissions(file_path)
            for dir_path in directory_path_list:
                change_permissions(dir_path)
        else:
            change_permissions(file_or_folder_path)


def adjust_dockerfile():
    """
    This will open Dockerfile and check if the entrypoint has been switched
    to the windows version (run_smilesmerge_in_container.bash) if not it will
    modify the Dockerfile to use the windows version of the script.

    This only should run on Windows OS.

    Change:
        ENTRYPOINT ["bash", "/smilesmerge/run_smilesmerge_in_container.bash"]
    To:
        # ENTRYPOINT ["bash", "/smilesmerge/run_smilesmerge_in_container.bash"]
        ENTRYPOINT ["bash", "/smilesmerge/run_smilesmerge_in_container_windows.bash"]
    """
    printout = ""
    normal_entry = "/smilesmerge/run_smilesmerge_in_container.bash"
    windows_entry = "/smilesmerge/run_smilesmerge_in_container_windows.bash"
    replacement_line = (
        'ENTRYPOINT ["bash", "/smilesmerge/run_smilesmerge_in_container_windows.bash"]\n'
    )
    print("Modifying the Dockerfile to run for Windows. Changing Entrypoint.")
    with open(os.path.abspath("Dockerfile"), "r") as f:
        for line in f.readlines():
            if "ENTRYPOINT" in line and normal_entry in line:
                if "#" not in line:
                    line = "# " + line + "\n"
                printout = printout + line
                printout = printout + replacement_line
                continue
            if "ENTRYPOINT" in line and windows_entry in line:
                continue
            printout = printout + line
    with open(os.path.abspath("Dockerfile"), "w") as f:
        f.write(printout)


def make_docker():
    """
    This will create the docker to run SMILESMerge.
    This is also where all of the files are copied into the image.

    If docker image can not be created it will raise an exception.
    """
    if os.name == "nt" or os.name == "ce":
        # so it's running under windows. multiprocessing disabled
        adjust_dockerfile()

    print("Creating new docker image for SMILESMerge")
    output_and_log_dir = os.path.abspath("output_and_log_dir") + os.sep
    log_file = "{}log.txt".format(output_and_log_dir)
    printout = (
        "\nAttempting to create the docker container. If 1st time running "
        + "this script it may take a few minutes. Output details are piped to: "
        + "{}\n".format(log_file)
    )

    print(printout)
    try:
        os.system("docker build -t smilesmerge . > {}".format(log_file))
    except:
        printout = (
            "\nCan not create a docker file. Please make sure to run the "
            + "script with sudo/administrative privileges.\nIf Linux/MacOS:\n"
            + "\t'sudo python smilesmerge_in_docker.py -j PATH/JSON.json'\n"
            + "If Windows:\n\tRun from bash and "
            + "docker enabled terminal with administrative privileges.\n"
            + "Please also make sure docker is installed on the system."
        )
        print(printout)
        raise Exception(printout)

    # Remove the temporary SMILESMerge directory
    shutil.rmtree("SMILESMerge")


def check_for_required_inputs(json_vars):
    """
    Confirm all the required inputs were provided.

    Required Variables go here.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name:
                           value}.

    Returns:
    :returns: dict json_vars: The updated json_vars with input-file paths
                              changed to the output dir/inputs/ subdirectory.
    """

    keys_from_input = list(json_vars.keys())

    list_of_required_inputs = [
        "root_output_folder",
        "source_compound_file",
    ]

    missing_variables = []
    for variable in list_of_required_inputs:
        if variable in keys_from_input:
            continue
        missing_variables.append(variable)

    if len(missing_variables) != 0:
        printout = "\nRequired variables are missing from the input. A description \
            of each of these can be found by running python ./RunSMILESMerge -h"
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
    json_vars["root_output_folder"] = os.path.abspath(json_vars["root_output_folder"])
    json_vars["source_compound_file"] = os.path.abspath(
        json_vars["source_compound_file"]
    )

    # Check root_output_folder exists
    if os.path.exists(json_vars["root_output_folder"]) is False:
        # If the output directory doesn't exist, then make ithe output
        # directory doesn't exist, then make it
        try:
            os.makedirs(json_vars["root_output_folder"])
            os.makedirs(json_vars["root_output_folder"] + os.sep + "inputs")
            change_permissions_recursively(json_vars["root_output_folder"])
        except:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

        if os.path.exists(json_vars["root_output_folder"]) is False:
            raise NotImplementedError(
                "root_output_folder could not be found and could not be created. \
                Please manual create desired directory or check input parameters"
            )

    if os.path.isdir(json_vars["root_output_folder"]) is False:
        raise NotImplementedError(
            "root_output_folder is not a directory. \
            Check your input parameters."
        )

    # Check source_compound_file exists
    if os.path.isfile(json_vars["source_compound_file"]) is False:
        raise NotImplementedError(
            "source_compound_file must be a tab delineated .smi file. \
            source_compound_file can not be found: \
            {}.".format(
                json_vars["source_compound_file"]
            )
        )
    if ".smi" not in json_vars["source_compound_file"]:
        raise NotImplementedError(
            "source_compound_file must be a \
            tab delineated .smi file."
        )

    # You need to copy the input files to the output directory, so it can be
    # easily access from the docker container. (JDD addition.)
    shutil.copy2(
        json_vars["source_compound_file"],
        json_vars["root_output_folder"]
        + os.sep + "inputs" + os.sep
        + os.path.basename(json_vars["source_compound_file"]),
    )
    json_vars["source_compound_file"] = "/Outputfolder/inputs/" + \
            os.path.basename(json_vars["source_compound_file"])

    return json_vars

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


def get_run_number(root_folder_path):
    """
    Determine run number for the new directory.
        Always Start a fresh new run.
            -If no previous runs exist in the root_folder_path
            -If there are previous runs in the root_folder_path
                incremental increasing the name by 1 from the last
                run in the same output directory.
    Inputs:
    :param str root_folder_path: is the path of the root output folder. We will
        make a directory within this folder to store our output files
    Returns:
    :returns: str run_num: the string of the run number "Run_*"
    """

    folder_name_path = root_folder_path + "Run_"
    print(folder_name_path)

    last_run_number = find_previous_runs(folder_name_path)

    if last_run_number is None:
        # There are no previous simulation runs in this directory
        print("There are no previous runs in this directory.")
        print("Starting a new run named Run_0.")

        # make a folder for the new generation
        run_number = 0

    else:
        # Start a new fresh simulation
        # Make a directory for the new run by increasing run number by +1
        # from last_run_number
        run_number = last_run_number + 1

    folder_name_path = root_folder_path + "Run_" + str(run_number) + os.sep
    print("The Run number is: ", run_number)
    print("The Run folder path is: ", folder_name_path)
    print("")
    return "Run_{}".format(run_number)


def get_output_folder(json_vars):
    """
    Find the folder for where to place output runs on host system.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    Returns:
    :returns: str root_output_folder: the string of the directory for
        puting output folders
    :returns: str run_num: the string of the run number "Run_*"
    """

    root_output_folder = os.path.abspath(json_vars["root_output_folder"]) + os.sep
    change_permissions_recursively(root_output_folder)
    run_num = get_run_number(root_output_folder)

    return root_output_folder, run_num


def move_files_to_temp_dir(json_vars):
    """
    This will move all files needed to a temp_user_files directory and will created a modified
    json_vars dict called docker_json_vars which will be used for pathing within
    the docker.

    Inputs:
    :param dict json_vars: The parameters. A dictionary of {parameter name: value}.
    """

    docker_json_vars = {}
    # make or remove and make the temp_user_files dir
    temp_dir_path = os.path.abspath("temp_user_files") + os.sep
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)
    change_permissions_recursively(temp_dir_path)

    # make or remove and make an output_and_log_dir
    output_and_log_dir = os.path.abspath("output_and_log_dir") + os.sep
    if os.path.exists(output_and_log_dir):
        shutil.rmtree(output_and_log_dir)
    os.mkdir(output_and_log_dir)
    change_permissions_recursively(output_and_log_dir)

    print("copying files into temp directory: temp_user_files")
    # get files from json_vars
    for var_name in json_vars.keys():
        var_item = json_vars[var_name]
        if str(type(var_item)) not in ["<type 'unicode'>", "<type 'unicode'>"]:
            continue
        var_item = str(var_item)
        # This could be a different variable that is not a path
        if os.path.exists(var_item) is False:
            continue
        if var_name == "root_output_folder":
            continue
        basename = os.path.basename(var_item)
        temp_path = temp_dir_path + basename
        if os.path.isdir(var_item):
            shutil.copytree(var_item, temp_path)
            docker_json_vars[var_name] = "/UserFiles/" + basename + "/"
            continue

        if os.path.isfile(var_item):
            shutil.copyfile(var_item, temp_path)
            docker_json_vars[var_name] = "/UserFiles/" + basename

    for var_name in json_vars.keys():
        if var_name not in docker_json_vars.keys():
            docker_json_vars[var_name] = json_vars[var_name]

    # Set output folder
    docker_json_vars["root_output_folder"] = "/Outputfolder/"

    with open(temp_dir_path + "docker_json_vars.json", "w") as file_item:
        json.dump(docker_json_vars, file_item, indent=4)

    # update permissions so files can be manipulated without sudo/admin
    change_permissions_recursively(temp_dir_path)
    change_permissions_recursively(output_and_log_dir)

    # Copy over SMILESMerge files into a temp directory
    temp_smilesmerge_path = os.path.abspath("SMILESMerge") + os.sep

    script_dir = str(os.path.dirname(os.path.realpath(__file__)))
    smilesmerge_top_dir = str(os.path.dirname(script_dir))
    if os.path.exists(temp_smilesmerge_path):
        shutil.rmtree(temp_smilesmerge_path)
    os.mkdir(temp_smilesmerge_path)
    smilesmerge_top_dir = smilesmerge_top_dir + os.sep
    change_permissions_recursively(temp_smilesmerge_path)

    # Copy all files in SMILESMerge directory into a temp except the Docker folder
    for fol_to_copy in [
            "smilesmerge",
            "source_compounds",
            "accessory_scripts",
            "tutorial",
    ]:
        shutil.copytree(
            smilesmerge_top_dir + fol_to_copy, temp_smilesmerge_path + fol_to_copy
        )
    shutil.copyfile(
        smilesmerge_top_dir + "RunSMILESMerge.py", temp_smilesmerge_path + "RunSMILESMerge.py"
    )
    # Open permissions
    change_permissions_recursively(temp_smilesmerge_path)


def handle_json_info(vars):
    """
    This will open the json file.
        1) check that JSON file has basic info
            - source compound file.
        2) copy files to a temp directory
            -source compound file.
        3) make a JSON file with modified information for within docker
    Inputs:
    :param dict vars: Dictionary of User specified variables

    Returns:
    :param dict json_vars: Dictionary of User specified variables
    :returns: str root_output_folder: the string of the directory for
        puting output folders
    :returns: str run_num: the string of the run number "Run_*"
    """
    print("Handling files")

    json_file = vars["json_file"]
    if os.path.exists(json_file) is False:
        printout = "\njson_file is required. Can not find json_file: {}.\n".format(
            json_file
        )
        print(printout)
        raise Exception(printout)
    json_vars = json.load(open(json_file))
    json_vars = check_for_required_inputs(json_vars)
    move_files_to_temp_dir(json_vars)

    # get output folder
    outfolder_path, run_num = get_output_folder(json_vars)

    return json_vars, outfolder_path, run_num


def run_SMILESMerge_docker_main(vars):
    """
    This function runs the processing to:
        1) check that JSON file has basic info
            -source compound file.
        2) copy files to a temp directory
            -source compound file.
        3) make a JSON file with modified information for within docker
        4) Build docker image and link files to output folder
            -This includes an adjustment to the Dockerfile if
            running it on a Windows OS
        5) execute RunSMILESMerge.py from within the docker container
        6) export the files back to the final end dir

    Inputs:
    :param dict vars: Dictionary of User specified variables
    """
    printout = (
        "\n\nThis script builds a docker for SMILESMerge and runs SMILESMerge "
        + "within the docker. The setup may take a few minutes the first time being run "
        + "and SMILESMerge may take a long time depending on the settings.\n\n"
    )
    print(printout)

    # Check that we are in the correct directory if not raise exception
    script_dir = str(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if os.path.abspath(os.getcwd()) != os.path.abspath(script_dir):
        printout = "\nMust execute this script from this directory: {}\n".format(
            script_dir
        )
        printout = printout + "Before running please 'cd {}'\n".format(script_dir)
        print(printout)
        raise Exception(printout)

    # Run parts 1-3
    # 1) check that JSON file has basic info
    #     -source compound file.
    # 2) copy files to a temp directory
    #     -source compound file.
    # 3) make a JSON file with modified information for within docker
    json_vars, outfolder_path, run_num = handle_json_info(vars)

    # Run build docker image
    make_docker()

    # Run part 5) run SMILESMerge in the container
    print("\nRunning SMILESMerge in Docker")

    command = "docker run --rm -it -v {}:/Outputfolder/".format(outfolder_path)
    command = command + " smilesmerge  --name smilesmerge  --{}".format(run_num)
    # Execute SMILESMerge
    print(command)
    os.system(command)
    change_permissions_recursively(outfolder_path)
    print("SMILESMerge Results placed in: {}".format(outfolder_path))


PARSER = argparse.ArgumentParser()

# Allows the run commands to be submitted via a .json file.
PARSER.add_argument(
    "--json_file",
    "-j",
    metavar="param.json_file",
    required=True,
    help="Name of a json file containing all parameters. \
    Overrides other arguments. This takes all the parameters described in \
    RunSMILESMerge.py.",
)
PARSER.add_argument(
    "--override_sudo_admin_privileges",
    metavar="param.override_sudo_admin_privileges",
    default=False,
    help="Docker normally requires `sudo` (linux/macos) or `Administrator` \
    privileges (windows/cygwin). If an system does not have such privileges, \
    or does not require such privileges, setting this to True, will skip the \
    check for privileges. This variable is provided via commandline, and \
    IS NOT RECOMMENDED for most OS. ",
)
ARGS_DICT = vars(PARSER.parse_args())

print("")
print("BE SURE TO RUN THIS SCRIPT WITH SUDO (LINUX/MACOS) OR ADMINISTRATOR")
print("(WINDOWS) PRIVILEGES!")
print("")

# Check that this is running with appropriate privileges.
# i.e., sudo (linux/macos) or Administrator privileges (Windows/cygwin)
if ARGS_DICT["override_sudo_admin_privileges"] == False:
    if sys.platform.lower() in ["darwin", "linux", "linux2"]:
        if os.getuid() != 0:
            printout = "\n\nMust run this script with `sudo` privileges.\n\t"
            printout = printout + "Please retry running with `sudo` privileges.\n\n"
            print(printout)
            raise Exception(printout)
    elif sys.platform.lower() == "win32" or sys.platform.lower() == "cygwin":
        import ctypes
        if ctypes.windll.shell32.IsUserAnAdmin() != 1:
            printout = "\n\nMust run this script from a terminal with `Administrator` privileges.\n\t"
            printout = printout + "Please retry running with `Administrator` privileges.\n\n"
            print(printout)
            raise Exception(printout)
    else:
        print("")
        print("BE SURE TO RUN THIS SCRIPT WITH SUDO (LINUX/MACOS) OR ADMINISTRATOR")
        print("(WINDOWS) PRIVILEGES!")
        print("")
else:

    print("\n##############################################################")
    print("WARNING: Skipping check for privileges.")
    print("\tBE SURE TO RUN THIS SCRIPT WITH APPROPRIATE PRIVILEGES:")
    print("\tSUDO (LINUX/MACOS) OR ADMINISTRATOR (WINDOWS) PRIVILEGES!")
    print("\tFailure to do so may result in Docker failures.")
    print("##############################################################\n")


run_SMILESMerge_docker_main(ARGS_DICT)
