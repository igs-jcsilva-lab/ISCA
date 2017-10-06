
# Holds functions used by at least two scripts in the targeted assembly pipeline.

import os,errno

# Function to ensure that a directory is only created if it does not
# exist. Eliminates the race condition of a simple call for whether
# the directory exists.
# Argument:
# path = path to a directory to check its existence 
def make_directory(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Function to generate new FASTA files after trimming extensions.
# Arguments:
# file: name/path of file to be written
# header: header ID
# seq: sequence to be written
def write_fasta(file,header,seq):
    with open(file,'a') as out:
        out.write(">{0}\n".format(header))
        for j in range(0, len(seq), 60):
            out.write(seq[j:j+60] + "\n")
            