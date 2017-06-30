#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Targeted Assembly -- assign reads to individual locus directories
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


inputs:
  reads_dir:
    label: Reads directory made by build_workspace.cwl
    type: Directory
    inputBinding:
      prefix: "--reads_dir"

  filter:
    label: Either "yes" or "no" for removing discrepancies + multi-locus mapping reads
    type: string
    inputBinding:
      prefix: "--filter"

  paired_suffixes:
    label: Either "yes" or "no" for whether the reads are mapped to one another with suffixes like .1 and .2 and one wants to assess for concordancy. This is dependent on the aligner. Check the *read_map.tsv file and see if the first elements are by read pair (so no suffix) or individual read (each read has  suffix) and answer accordingly
    type: string
    inputBinding:
      prefix: "--paired_suffixes"

  ab_read_map:
    label: Path to *_read_map.tsv output from analyze_bam.cwl
    type: File
    inputBinding:
      prefix: "--ab_read_map"

  reads1:
    label: Path to where the first paired fastq.gz file
    type: File
    inputBinding:
      prefix: "--fastq1"

  reads2:
    label: Path to where the second paired fastq.gz file
    type: File
    inputBinding:
      prefix: "--fastq2"


outputs:
  assigned_reads_dir:
    type: Directory
    outputBinding:
      outputEval: $(inputs.reads_dir)


baseCommand: ["/Library/Frameworks/Python.framework/Versions/3.5/bin/python3","/Users/jmatsumura/dev/targeted_assembly/bin/fastq_reads_to_fastq_alleles.py"]