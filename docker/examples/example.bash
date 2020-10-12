#!/bin/bash

# This script is an example for running SMILESMerge within a docker.
# To modify the protein, pocket dimensions, and GA paremeters... please 
# create a JSON file with the desired user variables.
# An example JSON is provided at: /SMILESMerge/docker/examples/sample_submit_SMILESMerge_docker.json


# Make sure we are in the /SMILESMerge/docker/ directory


# sudo should only be run in Linux or MacOS
# If Windows please instead just open the terminal with admin privileges
#   and omit the 'sudo'

sudo python ./smilesmerge_in_docker.py -j ./examples/sample_submit_SMILESMerge_docker.json
