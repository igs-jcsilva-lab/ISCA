

# This script establishes the necessary directory structures needed to run the
# pipeline. It requires the destination for where to place these directories. 
# Please make sure you have the necessary permissions to create directories at
# the specified location. 
#
# Run the script using a command like this:
# python3 setup.py -output_location /path/to/build_directories
#
# Author: James Matsumura

from shared_fxns import make_directory

def main():

    parser = argparse.ArgumentParser(description='Script to establish the necessary directory structures for the pipelines output.')
    parser.add_argument('-output_location', type=str, help='Path to build directories at.')
    args = parser.parse_args()

    # directories for the grid to output logs and errors
    make_directory("{0}/grid_out".formatargs.output_location)
    make_directory("{0}/grid_err".formatargs.output_location)

    make_directory("{0}/gsnap_idx".formatargs.output_location) # stores GSNAP index
    make_directory("{0}/smalt_idx".formatargs.output_location) # stores SMALT index
    make_directory("{0}/sam".formatargs.output_location) # alignment output
    make_directory("{0}/reads".formatargs.output_location) # individual read sets per locus
    make_directory("{0}/spades_assemblies".formatargs.output_location) # assembly method 1 results
    make_directory("{0}/hga_assemblies".formatargs.output_location) # assembly method 2 results
    make_directory("{0}/alignments".formatargs.output_location) # alignment results

if __name__ == '__main__':
    main()