#!/bin/bash

SEQTKPATH=/home/tdvorkina/soft/seqtk/
ideal_path_script=/home/tdvorkina/tmp/algorithmic-biology/assembler/src/projects/segal/aligner_stats/build_idealpaths.py
ideal_path_gen=/home/tdvorkina/tmp/algorithmic-biology/assembler/build/release/bin/idealreads_aligner

bwa index raw_data/C.elegans/ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa

mkdir tmp3
python2 $ideal_path_script raw_data/C.elegans/real_pacbio/real_pacbio.fasta raw_data/C.elegans/ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
                   /Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/graph_art/saves/simplification 77 pacbio
cut -f1 raw_data/C.elegans/real_pacbio/refseq_bwamem_real_pacbio.fasta_mapping.tsv > raw_data/C.elegans/real_pacbio/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/C.elegans/real_pacbio/realseq_bwamem_real_pacbio.fasta raw_data/C.elegans/real_pacbio/refseq_bwamem_mapping_names.txt > \
                        raw_data/C.elegans/real_pacbio/realpb_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/C.elegans/real_pacbio/realpb_mapped.fasta 10000 > benchmarkpaths_data/C.elegans/input/realpb_mapped_10000.fasta

python2 $ideal_path_script raw_data/C.elegans/sim_pacbio/sim_pacbio.fasta raw_data/C.elegans/ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
                   /Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/graph_art/saves/simplification 77 pacbio
cut -f1 raw_data/C.elegans/sim_pacbio/refseq_bwamem_sim_pacbio.fasta_mapping.tsv > raw_data/C.elegans/sim_pacbio/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/C.elegans/sim_pacbio/realseq_bwamem_sim_pacbio.fasta raw_data/C.elegans/sim_pacbio/refseq_bwamem_mapping_names.txt > \
                        raw_data/C.elegans/sim_pacbio/simpb_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/C.elegans/sim_pacbio/simpb_mapped.fasta 10000 > benchmarkpaths_data/C.elegans/input/simpb_mapped_10000.fasta

python2 $ideal_path_script raw_data/C.elegans/real_np/real_nanopore.fasta raw_data/C.elegans/ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
                   /Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/graph_art/saves/simplification 77 nanopore
cut -f1 raw_data/C.elegans/real_np/refseq_bwamem_real_nanopore.fasta_mapping.tsv > raw_data/C.elegans/real_np/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/C.elegans/real_np/realseq_bwamem_real_nanopore.fasta raw_data/C.elegans/real_np/refseq_bwamem_mapping_names.txt > \
                        raw_data/C.elegans/real_np/realnp_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/C.elegans/real_np/realnp_mapped.fasta 10000 > benchmarkpaths_data/C.elegans/input/realnp_mapped_10000.fasta

python2 $ideal_path_script raw_data/C.elegans/sim_np/simulated_reads_aligned.fasta raw_data/C.elegans/ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
                   /Sid/tdvorkina/gralign/C.elegans_synth/benchmarking/graph_art/saves/simplification 77 nanopore
cut -f1 raw_data/C.elegans/sim_np/refseq_bwamem_simulated_reads_aligned.fasta_mapping.tsv > raw_data/C.elegans/sim_np/refseq_bwamem_mapping_names.txt
$SEQTKPATH/seqtk subseq raw_data/C.elegans/sim_np/realseq_bwamem_simulated_reads_aligned.fasta raw_data/C.elegans/sim_np/refseq_bwamem_mapping_names.txt > \
                        raw_data/C.elegans/sim_np/simnp_mapped.fasta
$SEQTKPATH/seqtk sample -s100 raw_data/C.elegans/sim_np/simnp_mapped.fasta 10000 > benchmarkpaths_data/C.elegans/input/simnp_mapped_10000.fasta

rm -rf tmp3

cp raw_data/C.elegans/real_pacbio/refseq_bwamem_real_pacbio.fasta_mapping.tsv benchmarkpaths_data/C.elegans/input/realpb_true_mapping.tsv
cp raw_data/C.elegans/sim_pacbio/refseq_bwamem_sim_pacbio.fasta_mapping.tsv benchmarkpaths_data/C.elegans/input/simpb_true_mapping.tsv
cp raw_data/C.elegans/real_np/refseq_bwamem_real_nanopore.fasta_mapping.tsv benchmarkpaths_data/C.elegans/input/realnp_true_mapping.tsv
cp raw_data/C.elegans/sim_np/refseq_bwamem_simulated_reads_aligned.fasta_mapping.tsv benchmarkpaths_data/C.elegans/input/simnp_true_mapping.tsv
