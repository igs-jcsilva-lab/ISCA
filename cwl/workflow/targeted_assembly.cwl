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
    label: Path to allow Python3 to be found in the ENV
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
  filter:
    label: Either "yes" or "no" for removing discrepancies + multi-locus mapping reads
    type: string
  first_paired_suffixes:
    label: Either "yes" or "no" for whether the reads are mapped to one another with suffixes like .1 and .2 and one wants to assess for concordancy. This is dependent on the aligner. Check the *read_map.tsv file and see if the first elements are by read pair (so no suffix) or individual read (each read has  suffix) and answer accordingly
    type: string
  first_prefix:
    label: Name of the prefix to yield the two maps (one read-based and one reference-based)
    type: string
  first_assmb_map:
    label: Name of the map to create which maps a reference locus to an int ID
    type: string

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int
  aligner_threads:
    label: Number of threads to use for alignment
    type: int
  recruitment_threshold:
    label: Percent cut-off for how many matches a read must hit (uses whichever is longer the read or reference as denominator)
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

  first_phase_two_assmb_map:
    type: File
    outputSource: first_phase_two/assmb_map
  first_phase_two_read_map:
    type: File
    outputSource: first_phase_two/read_map
  first_phase_two_ref_map:
    type: File
    outputSource: first_phase_two/ref_map
  first_phase_renamed_reads_dir:
    type: Directory
    outputSource: first_phase_two/renamed_reads_dir
  first_phase_renamed_assmb_dir:
    type: Directory
    outputSource: first_phase_two/renamed_assmb_dir

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

  first_phase_two:
    run: phase_two.cwl
    in:
      threshold: recruitment_threshold
      prefix: first_prefix
      filter: filter
      paired_suffixes: first_paired_suffixes
      reads1: reads1
      reads2: reads2
      outfile: first_assmb_map
      ea_map: phase_one/ea_map
      reads_dir: phase_one/first_reads
      assmb_path: phase_one/first_spades_assemblies
      sam: gsnap/gsnap_sam
      python3_lib: python3_lib
    out: [
      read_map,
      ref_map,
      renamed_reads_dir,
      renamed_assmb_dir,
      assmb_map
    ]