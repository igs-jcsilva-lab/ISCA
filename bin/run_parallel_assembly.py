

# This script acts similar to an SGE array job where it queues tasks as 
# resources on the system as they become available. 
#
# PLEASE BE CAREFUL WHEN RUNNING THIS SCRIPT. 
#
# It is certainly possible to consume all the resources on your machine 
# if you request too much in terms of threads and memory. How many threads
# you need for each assembly is going to vary depending on your 
# dataset. Please start conservative with your estimates and if SPAdes 
# is prohibitively slow or gets OOM errors then increment up. 
#
# Generally, with default read inputs SPAdes can succeed with 1:1
# GB:thread count. The amount of memory per thread will increase
# to 3:1 when running HGA. These are only estimates though, YMMV. 
#
# ******
# total number of cores/cpus used by this script is determined by: (--number_of_jobs) + (--number_of_jobs) * (--threads_per_job)
# rough total memory used by this script is determined by: (--number_of_jobs) * (--memory_per_job)
# ******
#
# Run the script using a command like this:
# python3 run_parallel_assembly.py --assmb_step (SPAdes|HGA|SB) --number_of_jobs 8 --assmb_map /path/to/map.tsv \
# --threads_per_job 5 --memory_per_job 8 --reads_dir /path/to/reads_dir --spades_install /path/to/spades/bin \
# -rd /path/to/reads_dir -ap /path/to/assemblies
#
# Author: James Matsumura

import re,argparse,os,subprocess
import multiprocessing as mp

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    parser.add_argument('--assmb_map', '-am', type=str, required=True, help='Path to map.tsv output from format_for_assembly.py or assembly_verdict.py.')
    parser.add_argument('--assmb_step', '-as', type=str, required=True, help='Either "SPAdes","HGA", or "SB".q')
    parser.add_argument('--assmb_dir', '-ap', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('--number_of_jobs', '-n', required=True, help='Number of cores to call for individual processes.')
    parser.add_argument('--threads_per_job', '-t', required=True, help='Number of threads to invoke for each job.')
    parser.add_argument('--memory_per_job', '-m', required=True, help='Number of threads to invoke for each job.')
    parser.add_argument('--reads_dir', '-rd', type=str, required=True, help='Base path to where the reads are stored.')
    parser.add_argument('--spades_install', '-si', type=str, required=True, help='Path to the the location of the SPAdes install.')

    args = parser.parse_args()

    # Set up the multiprocessing manager, pool, and queue
    manager = mp.Manager()
    pool = mp.Pool(args.number_of_jobs)
    jobs = []

    # Need to iterate over the map generated from SPAdes step.
    with open(args.assmb_map,'r') as loc_map:
        for line in loc_map:
            line = line.rstrip()
            ele = line.split('\t')
            locus = ele[1] # no matter the map (SPAdes or HGA), the second column is what we want

            reads = "{0}/{1}/reads.fastq".format(args.reads_dir,locus)
            assembly_out = "{0}/{1}".format(args.assmb_dir,locus)

            if args.assmb_type == "SPAdes":
                jobs.append(pool.apply_async(spades_assemble, (args.spades_install,reads,args.memory_per_job,args.threads_per_job,assembly_out)))

    # Get all the returns from the apply_async function.
    for job in jobs:
        job.get()
    
    pool.close() #  Tell the queue it's done getting new jobs
    pool.join() # Make sure these new jobs are all finished

# Worker for assembling via SPAdes
# Arguments:
# spades = path to SPAdes bin
# reads = paired reads file to assemble
# memory = number of GB to reserve for SPAdes
# threads = number of threads to reserve for SPAdes 
# assmb_dir = where to place these assemblies
def spades_assemble(spades,reads,memory,threads,assmb_dir):
    
    spades_exe = "{0}/spades.py".format(spades)

    command = "{0} -m {1} -t {2} --careful --pe1-12 {3} -o {4}".format(spades_exe,memory,threads,reads,assembly_out)

    subprocess.call(command,shell=True)


if __name__ == '__main__':
    main()
