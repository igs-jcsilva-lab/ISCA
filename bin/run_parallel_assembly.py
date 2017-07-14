

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

import re,argparse,os,subprocess,shutil
import multiprocessing as mp

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    parser.add_argument('--assmb_map', '-am', type=str, required=True, help='Path to map.tsv output from format_for_assembly.py or assembly_verdict.py.')
    parser.add_argument('--assmb_step', '-as', type=str, required=True, help='Either "SPAdes","HGA", or "SB".q')
    parser.add_argument('--assmb_path', '-ap', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('--number_of_jobs', '-n', type=int, required=True, help='Number of cores to call for individual processes.')
    parser.add_argument('--threads_per_job', '-t', type=int, required=True, help='Number of threads to invoke for each job.')
    parser.add_argument('--memory_per_job', '-m', required=True, help='Memory, in GB, to use per job.')
    parser.add_argument('--reads_dir', '-rd', type=str, required=True, help='Base path to where the reads are stored.')
    parser.add_argument('--spades_install', '-si', type=str, required=True, help='Path to the the location of the SPAdes install.')
    # making all HGA options optional as to not flood SPAdes only assembly with reqs
    parser.add_argument('--HGA_install', '-hi', type=str, required=False, help='Path to the the location of HGA.py.')
    parser.add_argument('--python2_install', '-pi', type=str, required=False, help='Path to the the location of the Python2 executable.')
    parser.add_argument('--velvet_install', '-vi', type=str, required=False, help='Path to the the location of the velvet install.')
    parser.add_argument('--partitions', '-p', type=int, required=True, help='Number of partitions to split on during HGA.')

    args = parser.parse_args()

    if (args.number_of_jobs * args.threads_per_job) >= mp.cpu_count():
        print("Number of CPUs*threads requested will consume everything. You will want to allow at least one core to be used by processes outside this program.")
        sys.exit(0)

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

            assembly_out = "{0}/{1}".format(args.assmb_path,locus)

            if args.assmb_step == "SPAdes":
                reads = "{0}/{1}/reads.fastq.gz".format(args.reads_dir,locus)
                jobs.append(pool.apply_async(spades_assemble, (args.spades_install,reads,args.memory_per_job,args.threads_per_job,assembly_out)))

            elif args.assmb_step == "HGA":
                reads = "{0}/{1}/reads.fastq".format(args.reads_dir,locus)
                jobs.append(pool.apply_async(hga_assemble, (args.HGA_install,reads,assembly_out,args.python2_install,args.velvet_install,args.spades_install,args.threads_per_job,args.partitions,args.memory_per_job)))

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
    
    spades_exe = "{0}/bin/spades.py".format(spades)

    command = "{0} -m {1} -t {2} --careful --pe1-12 {3} -o {4}".format(spades_exe,memory,threads,reads,assmb_dir)

    subprocess.call(command.split())

    for path,dirs,files in os.walk(assmb_dir):
        for name in dirs: # directories need to go first
            shutil.rmtree(os.path.join(path,name))
        for name in files:
            if 's.fasta' not in name:
                if 'spades.log' not in name:
                    os.remove(os.path.join(path,name))

# Worker for assembling via HGA
# Arguments:
# HGA = path to HGA.py
# reads = paired reads file to assemble
# assmb_dir = location to output the assemblies
# python2 = location of python2 executable
# velvet = location of velvet executable
# spades = location of spades executable
# threads = number of threads to run SPAdes with
# partitions = number of partitions to randomly assign reads from 
# memory = how much memory to limit SPAdes to 
def hga_assemble(HGA,reads,assmb_dir,python2,velvet,spades,threads,partitions,memory):

    spades_exe = "{0}/bin".format(spades)

    # By default, the reads are gunzipped so need to uncompress. 
    subprocess.call("gunzip {0}.gz".format(reads).split())

    # Cleanup a bit of the unnecessary files. 
    remove_fastqs = "{0}/*.fastq".format(assmb_dir)
    remove_partitions = "{0}/part*assembly".format(assmb_dir)

    #
    # NOTE HARDCODED -PA SPAdes #
    #
    command = ("{0} {1} -velvet {2} -spades {3} -PA SPAdes -P12 {4} -R12 {5}"
         " -ins 150 -std 50 -Pkmer 31 -Rkmer 81 -t {6} -P {7} -out {8} -m {9}"
         " && rm {10} && rm -rf {11}"
         .format(python2,HGA,velvet,spades_exe,reads,reads,threads,partitions,assmb_dir,memory,remove_fastqs,remove_partitions)
    )

    subprocess.call(command.split())


if __name__ == '__main__':
    main()
