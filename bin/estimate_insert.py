

# This script automates a few steps in order to estimate insert size for 
# running the second assembly method.  
#
# Run the script using a command like this:
# python3 estimate_insert.py -samtools /path/to/samtools -java /path/to/java -picard /path/to/picard -sam_file /path/to/original/sam 
#
# Author: James Matsumura

import re,argparse,subprocess,sys

def main():

    parser = argparse.ArgumentParser(description='Script to automate samtools+picard for insert size estimation of paired reads.')
    parser.add_argument('--samtools', '-s', type=str, required=True, help='Path to SAMtools install.')
    parser.add_argument('--java', '-j', type=str, required=True, help='Path to Java install.')
    parser.add_argument('--picard', '-p', type=str, required=True, help='Path to Picard tools install.')
    parser.add_argument('--sam_file', '-sf', type=str, required=True, help='Alignment file (must end in .sam).')
    parser.add_argument('--out_dir', '-o', type=str, required=True, help='Location to place the insert.stats results file')
    args = parser.parse_args()

    # compress and then remove SAM for the sake of space
    bam_file = args.sam_file.replace(".sam",".bam")
    sam_to_bam = '{0} view -bS {1} > {2}'.format(args.samtools,args.sam_file,bam_file)
    subprocess.call(sam_to_bam,shell=True)
    subprocess.call("rm {0}".format(args.sam_file),shell=True)

    # sort the BAM file and remove non sorted
    sorted_bam_file = bam_file.replace(".bam","_sorted")
    sort_bam = '{0} sort {1} {2}'.format(args.samtools,bam_file,sorted_bam_file)
    subprocess.call(sort_bam,shell=True)
    subprocess.call("rm {0}".format(bam_file),shell=True)

    # use Picard tools to estimate insert
    estimate_insert = '{0} -jar {1}/CollectInsertSizeMetrics.jar I={2}.bam O={3}/insert.stats VALIDATION_STRINGENCY=LENIENT HISTOGRAM_file={3}/histo.out'.format(args.java,args.picard,sorted_bam_file,args.out)
    subprocess.call(estimate_insert,shell=True)

if __name__ == '__main__':
    main()