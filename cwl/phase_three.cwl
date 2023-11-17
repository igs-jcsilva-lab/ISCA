#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Phase three of the workflow consists of alignment.cwl, alignment_assessment.cwl, assembly_verdict.cwl, and get_final_sequences.cwl
class: Workflow


requirements:
  - class: InlineJavascriptRequirement


inputs:
  ea_map:
    label: Path to the output from extract_alleles.py
    type: File

  assmb_map:
    label: Path to map from format_for_assembly.cwl
    type: File

  emboss_tool:
    label: Path to install directory of EMBOSS needle/water executable (e.g. /path/to/packages/emboss/bin/[needle|water])
    type: File

  original_fsa:
    label: Path to where the initial FASTA file generated from the pipeline is
    type: File

  original_buffered_fsa:
    label: Path to where the initial, with any padding, FASTA file generated from the pipeline is
    type: File

  assmb_path:
    label: Path to the the directory to initialize directories for all the assembly output
    type: Directory

  align_path:
    label: Path to output directory for all these alignments.
    type: Directory

  prefix:
    label: Prefix for the output FASTA files to generate in current or existing directory
    type: string

  assmb_type:
    label: Either "SPAdes" or "HGA". Determines how many assembled sequences are aligned to
    type: string

  ivc_outfile:
    label: Name of output file
    type: File

  sequences_outfile:
    label: Name of output file
    type: string

  best_only:
    label: Either "yes" or "no" for whether to report stats of only the best alignment or all alignments
    type: string

  groupby:
    label: Get sequences by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds"
    type: string

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string

  number_of_jobs:
    label: Number of alignment jobs to spawn
    type: int

  threshold:
    label: Minimum threshold of %ID that needs to be met to pass final assembly
    type: int

  min_align_len:
    label: Minimum alignment length ratio of assembled sequence as cutoff to pull sequence or not
    type: double

  assmb_stdout:
    label: Prior step stdout used as workaround for outputEval.
    type: File

outputs:

  ids_v_cov:
    type: File
    outputSource: alignment_assessment/ids_v_cov

  leftovers:
    type: File
    outputSource: assembly_verdict/leftovers

  hga_assmb_map:
    type: File
    outputSource: assembly_verdict/hga_assmb_map

  final_sequences:
    type: File
    outputSource: get_final_sequences/final_sequences


steps:
  alignment:
    run: alignment.cwl
    in:
      ea_map: ea_map
      assmb_map: assmb_map
      original_fsa: original_fsa
      assmb_path: assmb_path
      emboss_tool: emboss_tool
      align_path: align_path
      number_of_jobs: number_of_jobs
      assmb_type: assmb_type
      python3_lib: python3_lib
    out: [stdout]

  alignment_assessment:
    run: alignment_assessment.cwl
    in:
      assmb_map: assmb_map
      align_path: align_path
      number_of_jobs: number_of_jobs
      assmb_type: assmb_type
      ivc_outfile: ivc_outfile
      best_only: best_only
      python3_lib: python3_lib
      alignment_stdout: alignment/stdout
    out: [ids_v_cov]

  assembly_verdict:
    run: assembly_verdict.cwl
    in:
      ivc: alignment_assessment/ids_v_cov
      threshold: threshold
      original_buffered_fsa: original_buffered_fsa
      original_assmb_map: assmb_map
      prefix: prefix
      python3_lib: python3_lib
    out: [
      leftovers,
      hga_assmb_map
      ]
    
  get_final_sequences:
    run: get_final_sequences.cwl
    in:
      ivc: alignment_assessment/ids_v_cov
      align_path: align_path
      groupby: groupby
      threshold: threshold
      original_fsa: original_fsa
      outfile: sequences_outfile
      min_align_len: min_align_len      
      ea_map: ea_map
      python3_lib: python3_lib
    out: [final_sequences]
