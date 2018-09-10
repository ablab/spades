#!/bin/bash

SEQTKPATH=/home/tdvorkina/soft/seqtk/
ideal_path_script=/home/tdvorkina/tmp/algorithmic-biology/assembler/src/projects/segal/aligner_stats/build_idealpaths.py
ideal_path_gen=/home/tdvorkina/tmp/algorithmic-biology/assembler/build/release/bin/idealreads_aligner

bwa index raw_data/E.coli/ref/MG1655-K12.fasta

mkdir tmp3
python2 $ideal_path_script raw_data/E.coli/real_pacbio/real_pacbio.fasta raw_data/E.coli/ref/MG1655-K12.fasta \
                   /Sid/tdvorkina/gralign/E.coli_synth/benchmarking/graph_art/saves/simplification 77 pacbio
cut -f1 raw_data/E.coli/real_pacbio/refseq_bwamem_real_pacbio.fasta_mapping.tsv > raw_data/E.coli/real_pacbio/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/E.coli/real_pacbio/realseq_bwamem_real_pacbio.fasta raw_data/E.coli/real_pacbio/refseq_bwamem_mapping_names.txt > \
                        raw_data/E.coli/real_pacbio/realpb_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/E.coli/real_pacbio/realpb_mapped.fasta 10000 > benchmarkpaths_data/E.coli/input/realpb_mapped_10000.fasta
cp raw_data/E.coli/real_pacbio/refseq_bwamem_real_pacbio.fasta_mapping.tsv  benchmarkpaths_data/E.coli/input/realpb_true_mapping.tsv

python2 $ideal_path_script raw_data/E.coli/sim_pacbio/sd_0001.fasta raw_data/E.coli/ref/MG1655-K12.fasta \
                   /Sid/tdvorkina/gralign/E.coli_synth/benchmarking/graph_art/saves/simplification 77 pacbio
cut -f1 raw_data/E.coli/sim_pacbio/refseq_bwamem_sd_0001.fasta_mapping.tsv > raw_data/E.coli/sim_pacbio/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/E.coli/sim_pacbio/realseq_bwamem_sd_0001.fasta raw_data/E.coli/sim_pacbio/refseq_bwamem_mapping_names.txt > \
                        raw_data/E.coli/sim_pacbio/simpb_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/E.coli/sim_pacbio/simpb_mapped.fasta 10000 > benchmarkpaths_data/E.coli/input/simpb_mapped_10000.fasta
cp raw_data/E.coli/sim_pacbio/refseq_bwamem_sd_0001.fasta_mapping.tsv benchmarkpaths_data/E.coli/input/simpb_true_mapping.tsv

python2 $ideal_path_script raw_data/E.coli/real_np/nanopore_R9.fasta raw_data/E.coli/ref/MG1655-K12.fasta \
                   /Sid/tdvorkina/gralign/E.coli_synth/benchmarking/graph_art/saves/simplification 77 nanopore
cut -f1 raw_data/E.coli/real_np/refseq_bwamem_nanopore_R9.fasta_mapping.tsv > raw_data/E.coli/real_np/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/E.coli/real_np/realseq_bwamem_nanopore_R9.fasta raw_data/E.coli/real_np/refseq_bwamem_mapping_names.txt > \
                        raw_data/E.coli/real_np/realnp_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/E.coli/real_np/realnp_mapped.fasta 10000 > benchmarkpaths_data/E.coli/input/realnp_mapped_10000.fasta
cp raw_data/E.coli/real_np/refseq_bwamem_nanopore_R9.fasta_mapping.tsv benchmarkpaths_data/E.coli/input/realnp_true_mapping.tsv

python2 $ideal_path_script raw_data/E.coli/sim_np/simulated_reads_aligned.fasta raw_data/E.coli/ref/MG1655-K12.fasta \
                   /Sid/tdvorkina/gralign/E.coli_synth/benchmarking/graph_art/saves/simplification 77 nanopore
cut -f1 raw_data/E.coli/sim_np/refseq_bwamem_simulated_reads_aligned.fasta_mapping.tsv > raw_data/E.coli/sim_np/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/E.coli/sim_np/realseq_bwamem_simulated_reads_aligned.fasta raw_data/E.coli/sim_np/refseq_bwamem_mapping_names.txt > \
                        raw_data/E.coli/sim_np/simnp_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/E.coli/sim_np/simnp_mapped.fasta 10000 > benchmarkpaths_data/E.coli/input/simnp_mapped_10000.fasta
cp raw_data/E.coli/sim_np/refseq_bwamem_simulated_reads_aligned.fasta_mapping.tsv benchmarkpaths_data/E.coli/input/simnp_true_mapping.tsv

rm -rf tmp3


