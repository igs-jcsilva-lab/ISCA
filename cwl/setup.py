#!/usr/bin/env python3

# Script to initialize the CWL files so that they can be run in the user's 
# particular environment. 
#
# Requires that conf.py in this same directory is set with the correct values. 
#
# Run the script using a command like this:
# setup.py --path_to_cwl /path/to/targeted_assembly/cwl
#
# Author: James Matsumura

import argparse,os,fileinput
from conf import TARGETED_ASSEMBLY_BIN,PYTHON3_EXE,GMAP_GSNAP_BIN,SMALT_BIN

def main():

    parser = argparse.ArgumentParser(description='Script to set up CWL files with the correct paths.')
    parser.add_argument('--path_to_cwl', '-p', type=str, default='.', required=False, help='Path to the the base directory where all the CWL files are stored.')

    args = parser.parse_args()

    replacements = {
        'TARGETED_ASSEMBLY_BIN': TARGETED_ASSEMBLY_BIN,
        'PYTHON3_EXE': PYTHON3_EXE,
        'GMAP_GSNAP_BIN': GMAP_GSNAP_BIN,
        'SMALT_BIN': SMALT_BIN
    }

    for file in os.listdir(args.path_to_cwl):
        if os.path.isfile(file) and file.endswith('cwl'):
            with fileinput.FileInput(file, inplace=True, backup='.bak') as cwl_file:
                # iterate over each line, and in each line replace any of the
                # bin values with those set in the conf
                for line in cwl_file:
                    for k,v in replacements.items():
                        line = line.replace(k,v)
                    print(line, end='')


if __name__ == '__main__':
    main()