#!/usr/bin/env python3

"""
A script to setup numerous grid runs of the targeted assembly pipeline for an 
SGE environment. Copy the preparation_file_examples directory and modify the 
contents of those files accordingly except for those that start with 
PLACEHOLDER_*. 

Input: 
    1. A CSV file with each line containing: /path/to/reads1,/path/to/reads2
    2. An input directory which mirrors the contents of preparation_file_examples
    directory found in this repository location. Contents of the files should be 
    modified except for those elements that start with PLACEHOLDER_*.
    3. Location of an output directory to have all these runs output to.

Output: 
    1. A set of directories at the specified path to populate with the results
    of the pipeline. Each will also contain its own *.sh and *.yml file for
    running the pipeline.
    2. A preparation_map.csv file which tells you which directory belongs to 
    which set of input reads. Columns are ID,qsub,reads1,reads2.

Usage:
    prepare_grid_runs.py -c read_info_file.csv -y master.yml -o /path/to/out_dir

Author: 
    James Matsumura
"""

import argparse
import errno
import os

def main():

    parser = argparse.ArgumentParser(description='Script to setup numerous targeted assembly runs in an SGE environment.')
    parser.add_argument('-c', type=str, required=True, help='Path to a two-column csv file with first column being the location of the reads1 and second being reads2.')
    parser.add_argument('-i', type=str, required=True, help='Input directory with cwl.sh,qsub,targeted_assembly.yml files present.')
    parser.add_argument('-o', type=str, required=True, help='Location to generate output directories.')
    args = parser.parse_args()

    make_directory(args.o)

    with open("{}/cwl.sh".format(args.i),'r') as sh_infile:
        sh_content = sh_infile.read()

    with open("{}/qsub".format(args.i),'r') as qsub_infile:
        qsub_content = qsub_infile.read().strip()

    with open("{}/targeted_assembly.yml".format(args.i),'r') as yml_infile:
        yml_content = yml_infile.read()

    with open("{}/preparation_map.csv".format(args.o),'w') as prep_map: 

        prep_map.write("{}\n".format((',').join(['ID','qsub','reads1','reads2'])))

        with open(args.c,'r') as reads_file: 

            run_id = 1
            
            for line in reads_file:
                reads = line.strip().split(',')
                reads1,reads2 = reads[0],reads[1]

                cur_out_dir = "{}/{}".format(args.o,run_id) 
                make_directory(cur_out_dir)
                    
                with open("{}/cwl.sh".format(cur_out_dir),'w') as sh_outfile:

                    content = sh_content                   
                    content = content.replace('PLACEHOLDER_OUTDIR',cur_out_dir)
                    content = content.replace('PLACEHOLDER_YML',"{}/targeted_assembly.yml".format(cur_out_dir))

                    sh_outfile.write(content)

                with open("{}/targeted_assembly.yml".format(cur_out_dir),'w') as yml_outfile:

                    content = yml_content
                    content = content.replace('PLACEHOLDER_READS1',reads1)
                    content = content.replace('PLACEHOLDER_READS2',reads2)

                    yml_outfile.write(content)

                qsub_cmd = qsub_content
                qsub_cmd = qsub_cmd.replace('PLACEHOLDER_CWL_SH',"{}/cwl.sh".format(cur_out_dir))
                qsub_cmd = qsub_cmd.replace('PLACEHOLDER_CWL_OUT',"{}/cwl.out".format(cur_out_dir))
                qsub_cmd = qsub_cmd.replace('PLACEHOLDER_CWL_ERR',"{}/cwl.err".format(cur_out_dir))

                prep_map.write("{}\n".format((',').join([str(run_id),qsub_cmd,reads1,reads2])))

                run_id += 1 


def make_directory(path):
    """ Takes in a path and tries to build that directory """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


if __name__ == '__main__':
    main()