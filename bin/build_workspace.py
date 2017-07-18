

# This script establishes the necessary directory structures needed to run the
# pipeline. It requires the destination for where to place these directories. 
# Please make sure you have the necessary permissions to create directories at
# the specified location. 
#
# Run the script using a command like this:
# python3 build_workspace.py --workspace_location /path/to/build_directories
#
# Author: James Matsumura

import argparse, urllib.request
from shared_fxns import make_directory

def main():

    parser = argparse.ArgumentParser(description='Script to establish the necessary directory structures for the pipelines output.')
    parser.add_argument('--workspace_location', '-wl', type=str, default='.', help='Path to build directories at.')
    args = parser.parse_args()

    # directories for the grid to output logs and errors
    make_directory("{0}/grid_out".format(args.workspace_location))
    make_directory("{0}/grid_err".format(args.workspace_location))

    make_directory("{0}/gsnap_idx".format(args.workspace_location)) # stores GSNAP index
    make_directory("{0}/smalt_idx".format(args.workspace_location)) # stores SMALT index
    make_directory("{0}/sam".format(args.workspace_location)) # alignment output

    # Aligner 1
    make_directory("{0}/first_reads".format(args.workspace_location)) # individual read sets per locus
    make_directory("{0}/first_spades_assemblies".format(args.workspace_location)) # assembly method 1 results
    make_directory("{0}/first_hga_assemblies".format(args.workspace_location)) # assembly method 2 results
    make_directory("{0}/first_alignments".format(args.workspace_location)) # alignment results
    make_directory("{0}/first_intermediary_end_results".format(args.workspace_location)) # results from assembly_verdict.py FIRST assembly method
    make_directory("{0}/first_final_end_results".format(args.workspace_location)) # results from assembly_verdict.py SECOND assembly method

    # Aligner 2
    make_directory("{0}/second_reads".format(args.workspace_location)) 
    make_directory("{0}/second_spades_assemblies".format(args.workspace_location)) 
    make_directory("{0}/second_hga_assemblies".format(args.workspace_location)) 
    make_directory("{0}/second_alignments".format(args.workspace_location)) 
    make_directory("{0}/second_intermediary_end_results".format(args.workspace_location))
    make_directory("{0}/second_final_end_results".format(args.workspace_location))

    # Pull HGA and SB from their own repos, note these are modified from their original implementations for this pipeline
    hga_url = 'https://raw.githubusercontent.com/jmatsumura/Hierarchical-Genome-Assembly-HGA/master/HGA.py'
    with urllib.request.urlopen(hga_url) as response, open("{0}/HGA.py".format(args.workspace_location), 'wb') as out_file:
        data = response.read() # a `bytes` object
        out_file.write(data)

    sb_url = 'https://raw.githubusercontent.com/jmatsumura/Scaffold_builder/master/scaffold_builder.py'
    with urllib.request.urlopen(sb_url) as response, open("{0}/scaffold_builder.py".format(args.workspace_location), 'wb') as out_file:
        data = response.read() # a `bytes` object
        out_file.write(data)

if __name__ == '__main__':
    main()