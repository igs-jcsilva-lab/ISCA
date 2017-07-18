#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Encapsulates the entirety of the targeted assembly pipeline. This means phase_one -> gsnap/smalt_workflow -> phase_two -> run_parallele_assembly -> phase_three
class: Workflow


requirements:
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: LD_LIBRARY_PATH
        envValue: $(inputs.python3_lib)

inputs:
  ea_input:
    label: Path to a TSV list for references and isolates
    type: File
  subset_list:
    label: File with loci to subset the sequences by, simply include all if all are desired 
    type: File
  reads1:
    label: Path to the first read pair file
    type: File
  reads2:
    label: Path to the second read pair file
    type: File

  python3_lib:
    label: Path tp allow Python3 to be found in the ENV
    type: string
  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string
  outfile:
    label: Name of the output FASTA file to generate in current or existing directory
    type: string
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int
  aligner_threads:
    label: Number of threads to use for alignment
    type: int


outputs:
  HGA:
    type: File
    outputSource: phase_one/HGA
  scaffold_builder:
    type: File
    outputSource: phase_one/scaffold_builder
  ea_map:
    type: File
    outputSource: phase_one/ea_map
  sequences:
    type: File
    outputSource: phase_one/sequences
  gsnap_idx:
    type: Directory
    outputSource: phase_one/gsnap_idx
  smalt_idx:
    type: Directory
    outputSource: phase_one/smalt_idx
  sam_dir:
    type: Directory
    outputSource: phase_one/smalt_idx
  first_reads:
    type: Directory
    outputSource: phase_one/first_reads
  first_spades_assemblies:
    type: Directory
    outputSource: phase_one/first_spades_assemblies
  first_hga_assemblies:
    type: Directory
    outputSource: phase_one/first_hga_assemblies
  first_alignments:
    type: Directory
    outputSource: phase_one/first_alignments
  first_end_results:
    type: Directory
    outputSource: phase_one/first_end_results
  second_reads:
    type: Directory
    outputSource: phase_one/second_reads
  second_spades_assemblies:
    type: Directory
    outputSource: phase_one/second_spades_assemblies
  second_hga_assemblies:
    type: Directory
    outputSource: phase_one/second_hga_assemblies
  second_alignments:
    type: Directory
    outputSource: phase_one/second_alignments
  second_end_results:
    type: Directory
    outputSource: phase_one/second_end_results

  gsnap_sam:
    type: File
    outputSource: gsnap/gsnap_sam


steps:
  phase_one:
    run: phase_one.cwl
    in:
      ea_input: ea_input
      gene_or_exon: gene_or_exon
      buffer: buffer
      outfile: outfile
      subset_list: subset_list
      python3_lib: python3_lib
    out: [
      gsnap_idx,
      smalt_idx,
      sam_dir,
      first_reads,
      first_spades_assemblies,
      first_hga_assemblies,
      first_alignments,
      first_end_results,
      second_reads,
      second_spades_assemblies,
      second_hga_assemblies,
      second_alignments,
      second_end_results,
      HGA,
      scaffold_builder,
      ea_map,
      sequences
    ]

  gsnap:
    run: gsnap_workflow.cwl
    in:
      threads: aligner_threads
      reads1: reads1
      reads2: reads2
      gsnap_genome: gsnap_genome
      sequences: phase_one/sequences
      gsnap_dir: phase_one/gsnap_idx
      python3_lib: python3_lib
    out: [
      gsnap_sam
    ]