#!/usr/bin/env cwl-runner
cwlVersion: v1.0
label: Encapsulates the entirety of the targeted assembly pipeline. This means phase_one -> (gsnap/smalt_workflow -> phase_two -> run_parallel_assembly -> phase_three  -> run_parallel_assembly -> phase_three) -- section wrapped in parentheses performs two iterations with different aligner at the start
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
  emboss_tool:
    label: Path to install directory of EMBOSS needle/water executable (e.g. /path/to/packages/emboss/bin/[needle|water])
    type: File
  python2_exe:
    label: Location of the Python2 executable
    type: File

  HGA_exe:
    label: Location of the HGA executable
    type: File
  SB_exe:
    label: Location of the scaffold_builder executable
    type: File

  spades_install:
    label: Location of the SPAdes installation
    type: Directory
  velvet_install:
    label: Location of the Velvet installation
    type: Directory
  samtools_install:
    label: Location of samtools
    type: Directory

  python3_lib:
    label: Path to allow Python3 to be found in the ENV
    type: string
  gene_or_exon:
    label: Either "gene" or "exon" for which sequences to pull
    type: string
  prefix:
    label: Prefix of the output FASTA file to generate in current or existing directory
    type: string
    default: "original_seqs"
  gsnap_genome:
    label: Name of the "genome" for GSNAP, really just a unique identifier for this index
    type: string
    default: "gsnap"
  smalt_prefix:
    label: SMALT prefix for the indices built in smalt_index.cwl
    type: string
    default: "smalt"
  filter:
    label: Either "yes" or "no" for removing discrepancies + multi-locus mapping reads
    type: string
  first_paired_suffixes:
    label: Either "yes" or "no" for whether the reads are mapped to one another with suffixes like .1 and .2 and one wants to assess for concordancy. This is dependent on the aligner. Check the *read_map.tsv file and see if the first elements are by read pair (so no suffix) or individual read (each read has  suffix) and answer accordingly. GSNAP always trims so leave this as "no". 
    type: string
    default: "no"
  second_paired_suffixes:
    label: Either "yes" or "no" for whether the reads are mapped to one another with suffixes like .1 and .2 and one wants to assess for concordancy. This is dependent on the aligner. SMALT will not trim this, so check your WGS input reads and see whether they end in ".1"/".2". 
    type: string
  first_prefix:
    label: Name of the prefix to yield the two maps (one read-based and one reference-based)
    type: string
    default: "first"
  second_prefix:
    label: Name of the prefix to yield the two maps (one read-based and one reference-based)
    type: string
    default: "second"
  first_assmb_map:
    label: Name of the map to create which maps a reference locus to an int ID
    type: string
    default: "first_assmb_map.tsv"
  second_assmb_map:
    label: Name of the map to create which maps a reference locus to an int ID
    type: string
    default: "second_assmb_map.tsv"
  spades_str:
    label: Static string to relate to the parallel assembly and alignment steps which algorithm is being handled
    type: string
    default: "SPAdes"
  hga_str:
    label: Static string to relate to the parallel assembly and alignment steps which algorithm is being handled
    type: string
    default: "HGA"
  sb_str:
    label: Static string to relate to the parallel assembly and alignment steps which algorithm is being handled
    type: string
    default: "SB"
  kmer:
    label: Comma separated list of odd, ascending numbers lesser than 128 to be used for SPAdes and HGA step as a backup
    type: string
    default: "7,11,21"
  best_only:
    label: Either "yes" or "no" for whether to report stats of only the best alignment or all alignments
    type: string
  groupby:
    label: Get sequences by loci, alleles/exons, or CDS (could also be exons if extracted at ea_map step), choose either "l", "ae", or "cds"
    type: string
  first_intermediary_sequences:
    label: Name of first round SPAdes assembled seqs file
    type: string
    default: "first_intermediary_sequences.fsa"
  first_final_sequences:
    label: Name of first round HGA assembled seqs file
    type: string
    default: "first_final_sequences.fsa"
  first_intermediary_prefix:
    label: Prefix for the leftover FASTA seqs at this point
    type: string
    default: "first_intermediary"
  first_final_prefix:
    label: Prefix for the leftover FASTA seqs at this point
    type: string
    default: "first_final"
  second_intermediary_sequences:
    label: Name of first round SPAdes assembled seqs file
    type: string
    default: "second_intermediary_sequences.fsa"
  second_final_sequences:
    label: Name of second round HGA assembled seqs file
    type: string
    default: "second_final_sequences.fsa"
  second_intermediary_prefix:
    label: Prefix for the leftover FASTA seqs at this point
    type: string
    default: "second_intermediary"
  second_final_prefix:
    label: Prefix for the leftover FASTA seqs at this point
    type: string
    default: "second_final"

  buffer:
    label: How much of a buffer to add to each end of the gene/exon, defaults to 0
    type: int
  aligner_threads:
    label: Number of threads to use for alignment
    type: int
  recruitment_threshold:
    label: Percent cut-off for how many matches a read must hit (uses whichever is longer the read or reference as denominator)
    type: int
  aligner_threshold:
    label: Minimum threshold of %ID that needs to be met to pass final assembly
    type: int
  threads_per_job:
    label: Number of threads to use for each assembly job
    type: int
  memory_per_job:
    label: How much memory to limit for each individual assembly job
    type: int
  number_of_jobs:
    label: Number of assembly jobs to spawn
    type: int
  partitions:
    label: Number of partitions to use in HGA
    type: int

  min_align_len:
    label: Minimum alignment length ratio of assembled sequence as cutoff to pull sequence or not
    type: double


outputs:
  first_ivc:
    type: File
    outputSource: phase_one/first_ivc
  second_ivc:
    type: File
    outputSource: phase_one/second_ivc
  ea_map:
    type: File
    outputSource: phase_one/ea_map
  buffered_sequences:
    type: File
    outputSource: phase_one/buffered_sequences
  unbuffered_sequences:
    type: File
    outputSource: phase_one/unbuffered_sequences
  gsnap_idx:
    type: Directory
    outputSource: phase_one/gsnap_idx
  smalt_idx:
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

  first_phase_two_assmb_map:
    type: File
    outputSource: first_phase_two/assmb_map
  first_phase_two_read_map:
    type: File
    outputSource: first_phase_two/read_map
  first_phase_two_ref_map:
    type: File
    outputSource: first_phase_two/ref_map

  first_intermediary_phase_three_am:
    type: File
    outputSource: first_intermediary_phase_three/hga_assmb_map
  first_intermediary_phase_three_fs:
    type: File
    outputSource: first_intermediary_phase_three/final_sequences

  first_final_phase_three_am:
    type: File
    outputSource: first_final_phase_three/hga_assmb_map
  first_final_phase_three_fs:
    type: File
    outputSource: first_final_phase_three/final_sequences

  second_phase_two_assmb_map:
    type: File
    outputSource: second_phase_two/assmb_map
  second_phase_two_read_map:
    type: File
    outputSource: second_phase_two/read_map
  second_phase_two_ref_map:
    type: File
    outputSource: second_phase_two/ref_map

  second_intermediary_phase_three_am:
    type: File
    outputSource: second_intermediary_phase_three/hga_assmb_map
  second_intermediary_phase_three_fs:
    type: File
    outputSource: second_intermediary_phase_three/final_sequences

  second_final_phase_three_am:
    type: File
    outputSource: second_final_phase_three/hga_assmb_map
  second_final_phase_three_fs:
    type: File
    outputSource: second_final_phase_three/final_sequences

  assembled_seqs:
    type: File
    outputSource: gather_sequences/all_seqs

  assembled_below_threshold:
    type: File
    outputSource: extract_reconstructions_step/assembled_below_threshold

steps:
  phase_one:
    run: phase_one.cwl
    in:
      ea_input: ea_input
      gene_or_exon: gene_or_exon
      buffer: buffer
      prefix: prefix
      subset_list: subset_list
      python3_lib: python3_lib
    out: [
      gsnap_idx,
      smalt_idx,
      first_reads,
      first_spades_assemblies,
      first_hga_assemblies,
      first_alignments,
      second_reads,
      second_spades_assemblies,
      second_hga_assemblies,
      second_alignments,
      first_ivc,
      second_ivc,
      ea_map,
      buffered_sequences,
      unbuffered_sequences
    ]

  gsnap:
    run: gsnap_workflow.cwl
    in:
      threads: aligner_threads
      reads1: reads1
      reads2: reads2
      gsnap_genome: gsnap_genome
      sequences: phase_one/buffered_sequences
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
      samtools_install: samtools_install
    out: [
      read_map,
      ref_map,
      assmb_map
    ]

  first_spades_assmb:
    run: run_parallel_assembly.cwl
    in:
      python2_exe: python2_exe
      spades_install: spades_install
      assmb_step: spades_str
      number_of_jobs: number_of_jobs
      threads_per_job: threads_per_job
      memory_per_job: memory_per_job
      kmer: kmer
      assmb_map: first_phase_two/assmb_map
      reads_dir: phase_one/first_reads
      assmb_path: phase_one/first_spades_assemblies
      python3_lib: python3_lib
    out: [
      stdout
    ]

  first_intermediary_phase_three:
    run: phase_three.cwl
    in:
      emboss_tool: emboss_tool
      min_align_len: min_align_len
      threshold: aligner_threshold
      assmb_type: spades_str
      number_of_jobs: aligner_threads
      best_only: best_only
      groupby: groupby
      ivc_outfile: phase_one/first_ivc
      sequences_outfile: first_intermediary_sequences
      prefix: first_intermediary_prefix
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      original_buffered_fsa: phase_one/buffered_sequences
      align_path: phase_one/first_alignments
      assmb_path: phase_one/first_spades_assemblies
      assmb_map: first_phase_two/assmb_map
      python3_lib: python3_lib
      assmb_stdout: first_spades_assmb/stdout
    out: [
      ids_v_cov,
      hga_assmb_map,
      final_sequences
    ]

  first_hga_assmb:
    run: run_parallel_assembly.cwl
    in:
      spades_install: spades_install
      velvet_install: velvet_install
      assmb_step: hga_str
      python2_exe: python2_exe
      HGA_exe: HGA_exe
      SB_exe: SB_exe
      number_of_jobs: number_of_jobs
      threads_per_job: threads_per_job
      memory_per_job: memory_per_job
      partitions: partitions
      kmer: kmer
      reads_dir: phase_one/first_reads
      assmb_map: first_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/first_hga_assemblies
      python3_lib: python3_lib
    out: [
      stdout
    ]

  first_sb_assmb:
    run: run_parallel_assembly.cwl
    in:
      assmb_step: sb_str
      python2_exe: python2_exe
      HGA_exe: HGA_exe
      SB_exe: SB_exe
      number_of_jobs: number_of_jobs
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      kmer: kmer
      reads_dir: phase_one/first_reads
      assmb_map: first_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/first_hga_assemblies
      python3_lib: python3_lib
      assmb_stdout: first_hga_assmb/stdout
    out: [
      stdout
    ]

  first_final_phase_three:
    run: phase_three.cwl
    in:
      emboss_tool: emboss_tool
      min_align_len: min_align_len
      threshold: aligner_threshold
      assmb_type: hga_str
      number_of_jobs: aligner_threads
      best_only: best_only
      groupby: groupby
      ivc_outfile: phase_one/first_ivc
      sequences_outfile: first_final_sequences
      prefix: first_final_prefix
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      original_buffered_fsa: phase_one/buffered_sequences
      align_path: phase_one/first_alignments
      assmb_map: first_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/first_hga_assemblies
      python3_lib: python3_lib
      assmb_stdout: first_sb_assmb/stdout
    out: [
      ids_v_cov,
      leftovers,
      hga_assmb_map,
      final_sequences
    ]

  smalt:
    run: smalt_workflow.cwl
    in:
      threads: aligner_threads
      reads1: reads1
      reads2: reads2
      smalt_prefix: smalt_prefix
      smalt_dir: phase_one/smalt_idx
      sequences: first_final_phase_three/leftovers
      python3_lib: python3_lib
    out: [
      smalt_sam
    ]

  second_phase_two:
    run: phase_two.cwl
    in:
      threshold: recruitment_threshold
      prefix: second_prefix
      filter: filter
      paired_suffixes: second_paired_suffixes
      reads1: reads1
      reads2: reads2
      outfile: second_assmb_map
      ea_map: phase_one/ea_map
      reads_dir: phase_one/second_reads
      assmb_path: phase_one/second_spades_assemblies
      sam: smalt/smalt_sam
      python3_lib: python3_lib
      samtools_install: samtools_install
    out: [
      read_map,
      ref_map,
      assmb_map
    ]

  second_spades_assmb:
    run: run_parallel_assembly.cwl
    in:
      python2_exe: python2_exe
      spades_install: spades_install
      assmb_step: spades_str
      number_of_jobs: number_of_jobs
      threads_per_job: threads_per_job
      memory_per_job: memory_per_job
      kmer: kmer
      assmb_map: second_phase_two/assmb_map
      reads_dir: phase_one/second_reads
      assmb_path: phase_one/second_spades_assemblies
      python3_lib: python3_lib
    out: [
      stdout
    ]

  second_intermediary_phase_three:
    run: phase_three.cwl
    in:
      emboss_tool: emboss_tool
      min_align_len: min_align_len
      threshold: aligner_threshold
      assmb_type: spades_str
      number_of_jobs: aligner_threads
      best_only: best_only
      groupby: groupby
      ivc_outfile: phase_one/second_ivc
      sequences_outfile: second_intermediary_sequences
      prefix: second_intermediary_prefix
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      original_buffered_fsa: phase_one/buffered_sequences
      align_path: phase_one/second_alignments
      assmb_path: phase_one/second_spades_assemblies
      assmb_map: second_phase_two/assmb_map
      python3_lib: python3_lib
      assmb_stdout: second_spades_assmb/stdout
    out: [
      ids_v_cov,
      hga_assmb_map,
      final_sequences
    ]

  second_hga_assmb:
    run: run_parallel_assembly.cwl
    in:
      spades_install: spades_install
      velvet_install: velvet_install
      assmb_step: hga_str
      python2_exe: python2_exe
      HGA_exe: HGA_exe
      SB_exe: SB_exe
      number_of_jobs: number_of_jobs
      threads_per_job: threads_per_job
      memory_per_job: memory_per_job
      partitions: partitions
      kmer: kmer
      reads_dir: phase_one/second_reads
      assmb_map: second_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/second_hga_assemblies
      python3_lib: python3_lib
    out: [
      stdout
    ]

  second_sb_assmb:
    run: run_parallel_assembly.cwl
    in:
      assmb_step: sb_str
      python2_exe: python2_exe
      HGA_exe: HGA_exe
      SB_exe: SB_exe
      number_of_jobs: number_of_jobs
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      kmer: kmer
      reads_dir: phase_one/second_reads
      assmb_map: second_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/second_hga_assemblies
      python3_lib: python3_lib
      assmb_stdout: second_hga_assmb/stdout
    out: [
      stdout
    ]

  second_final_phase_three:
    run: phase_three.cwl
    in:
      emboss_tool: emboss_tool
      min_align_len: min_align_len
      threshold: aligner_threshold
      assmb_type: hga_str
      number_of_jobs: aligner_threads
      best_only: best_only
      groupby: groupby
      ivc_outfile: phase_one/second_ivc
      sequences_outfile: second_final_sequences
      prefix: second_final_prefix
      ea_map: phase_one/ea_map
      original_fsa: phase_one/unbuffered_sequences
      original_buffered_fsa: phase_one/buffered_sequences
      align_path: phase_one/second_alignments
      assmb_map: second_intermediary_phase_three/hga_assmb_map
      assmb_path: phase_one/second_hga_assemblies
      python3_lib: python3_lib
      assmb_stdout: second_sb_assmb/stdout
    out: [
      ids_v_cov,
      leftovers,
      hga_assmb_map,
      final_sequences
    ]

  gather_sequences:
    run: combine_assembled.cwl
    in:
        first_final: first_final_phase_three/final_sequences
        second_final: second_final_phase_three/final_sequences
        python3_lib: python3_lib
    out: [
        all_seqs
    ]

  extract_reconstructions_step:
    run: extract_reconstructions.cwl
    in:
      python3_lib: python3_lib
      first_align_path: phase_one/first_alignments
      second_align_path: phase_one/second_alignments
      first_idvc : first_intermediary_phase_three/ids_v_cov
      first_ids_v_cov : first_final_phase_three/ids_v_cov
      second_idvc : second_intermediary_phase_three/ids_v_cov
      second_ids_v_cov : second_final_phase_three/ids_v_cov
      original_fsa: phase_one/unbuffered_sequences
    out: [
      assembled_below_threshold
    ]
