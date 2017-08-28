#!/usr/bin/env python3

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
# run_parallel_assembly.py --assmb_step (SPAdes|HGA|SB) --number_of_jobs 8 --assmb_map /path/to/map.tsv \
# --threads_per_job 5 --memory_per_job 8 --reads_dir /path/to/reads_dir --spades_install /path/to/spades/bin \
# -rd /path/to/reads_dir -ap /path/to/assemblies
#
# Author: James Matsumura

import re,argparse,os,subprocess,shutil,sys
import multiprocessing as mp
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(description='Script to assess EMBOSS Needle alignments, follows global_alignment.py.')
    # All are required, but different combinations for SPAdes, HGA, or SB.
    parser.add_argument('--assmb_map', '-am', type=str, required=True, help='Path to map.tsv output from format_for_assembly.py or assembly_verdict.py.')
    parser.add_argument('--assmb_step', '-as', type=str, required=True, help='Either "SPAdes","HGA", or "SB"')
    parser.add_argument('--assmb_path', '-ap', type=str, required=True, help='Path to the the directory preceding all the ref directories (e.g. for "/path/to/ref123" put "/path/to" as the input).')
    parser.add_argument('--number_of_jobs', '-n', type=int, required=True, help='Number of cores to call for individual processes.')
    parser.add_argument('--threads_per_job', '-t', default=1, type=int, required=False, help='Number of threads to invoke for each job.')
    parser.add_argument('--memory_per_job', '-m', required=False, help='Memory, in GB, to use per job.')
    parser.add_argument('--reads_dir', '-rd', type=str, required=False, help='Base path to where the reads are stored.')
    parser.add_argument('--spades_install', '-si', type=str, required=False, help='Path to the the location of the SPAdes install.')
    parser.add_argument('--HGA_install', '-hi', type=str, required=False, help='Path to the the location of HGA.py.')
    parser.add_argument('--SB_install', '-sbi', type=str, required=False, help='Path to the the location of scaffold_builder.py.')
    parser.add_argument('--python2_install', '-pi', type=str, required=False, help='Path to the the location of the Python2 executable.')
    parser.add_argument('--velvet_install', '-vi', type=str, required=False, help='Path to the the location of the velvet install.')
    parser.add_argument('--partitions', '-p', type=int, required=False, help='Number of partitions to split on during HGA.')
    parser.add_argument('--ea_map', '-eam', type=str, required=False, help='Path to the output from extract_alleles.py.')
    parser.add_argument('--original_fsa', '-of', type=str, required=False, help='Path to where the unbuffered FASTA from extract_sequences.py is.')

    args = parser.parse_args()

    # SB runs poorly in parallel, until source code is addressed run serially
    if args.assmb_step == "SB": 
        with open(args.assmb_map,'r') as loc_map:
            for line in loc_map:
                line = line.rstrip()
                ele = line.split('\t')
                locus = ele[1] # no matter the map (SPAdes or HGA), the second column is what we want
                locus_id = ele[0] # for SB, we care about the locus ID not the mapped number

                assembly_out = "{0}/{1}".format(args.assmb_path,locus) 

                sb_align(args.SB_install,assembly_out,locus_id,args.python2_install,args.ea_map,args.original_fsa)

        return # done at this point if at SB step              

    if (args.number_of_jobs * args.threads_per_job) >= mp.cpu_count():
        print("Number of CPUs*threads requested will consume everything. You will want to allow at least one core to be used by processes outside this program.")
        sys.exit(1)

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
            locus_id = ele[0] # for SB, we care about the locus ID not the mapped number

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

    keep_us = set(['contigs.fasta','spades.log'])
    clean_up(assmb_dir,keep_us)

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
    # Also, -ins and -std are set arbitrarily. This is because SPAdes is called
    # for partition assembly and velvet is only called for building inputs for
    # merged/combined assemblies. Thus, right now -ins/-std are not used. 
    #
    command = ("{0} {1} -velvet {2} -spades {3} -PA SPAdes -P12 {4} -R12 {5}"
         " -ins 150 -std 50 -Pkmer 31 -Rkmer 81 -t {6} -P {7} -out {8} -m {9}"
         " && rm {10} && rm -rf {11}"
         .format(python2,HGA,velvet,spades_exe,reads,reads,threads,partitions,assmb_dir,memory,remove_fastqs,remove_partitions)
    )

    subprocess.call(command.split())

    keep_us = set(['contigs.fasta','spades.log','HGA.log','HGA_combined','HGA_merged'])
    clean_up(assmb_dir,keep_us)

# Worker for assembling via HGA
# Arguments:
# SB = path to scaffold_builder.py
# assmb_dir = location to output the assemblies
# locus = reference locus to build scaffolds based on 
# python2 = location of python2 executable
# ea_map = output from extract_alleles.py
# fasta = file that contains the sequences referenced by ea_map
def sb_align(SB,assmb_dir,locus,python2,ea_map,fasta):

    # Now that we know where to do work, build the reference file for just the
    # relevant locus. First need to get all the correct sequences.
    allele_list = []
    with open(ea_map,'r') as i:
        for line in i:
            line = line.rstrip()
            elements = line.split('\t')

            # Need to handle the case where the reference locus is 
            # split into multiple like ABC123.1,ABC123.2,etc.
            if '.' in elements[0]:
                ea_locus = elements[0].split('.')[0]
            else:
                ea_locus = elements[0]

            if ea_locus == locus:
                for x in range(1,len(elements)): # grab all the alleles
                    allele_list.append(elements[x].split('|')[4])

    # Extract and write out the sequences. 
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta,"fasta"))

    # First check if HGA worked
    query = "{0}/HGA_combined/contigs.fasta".format(assmb_dir)

    # If HGA built contigs, build both F/R scaffolds to align to
    if os.path.isfile(query):

        f_fasta = "{0}/f_locus.fasta".format(assmb_dir)
        with open(f_fasta,'w') as out:
            for allele in allele_list:
                seq = seq_dict[allele]
                SeqIO.write(seq,out,"fasta")

        r_fasta = "{0}/r_locus.fasta".format(assmb_dir)
        with open(r_fasta,'w') as out:
            for allele in allele_list:
                seq = seq_dict[allele]
                seq.seq = seq.seq.reverse_complement()
                SeqIO.write(seq,out,"fasta")

        f_out = "{0}/f".format(assmb_dir)
        r_out = "{0}/r".format(assmb_dir)

        command = "{0} {1} -q {2} -r {3} -p {4}".format(python2,SB,query,f_fasta,f_out)
        subprocess.call(command.split())
        command = "{0} {1} -q {2} -r {3} -p {4}".format(python2,SB,query,r_fasta,r_out)
        subprocess.call(command.split())

# Function to delete some directories/files made during assembly that are
# extraneous and don't need to be captured in the output. 
# Arguments:
# relevant_dir = base directory to traverse and clean up
# keep_set = set of strings/values to skip over, this should contain both 
# directory and file names.
def clean_up(relevant_dir,keep_set):

    for path,dirs,files in os.walk(relevant_dir):
        for name in dirs: # directories need to go first
            if name not in keep_set:
                shutil.rmtree(os.path.join(path,name))
        for name in files:
            if name not in keep_set:
                os.remove(os.path.join(path,name))


if __name__ == '__main__':
    main()
