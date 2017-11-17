#!/bin/sh

FILE_ECOLI=$1
FILE_CELEGANS=$2
FILE_RUMEN=$3
PREFIX=$4

echo "--------------Ecoli--------------"
python2 ~/scripts/aligner_stats/aligner_stats.py ~/gralign/bwamem_aligned/ecoli/ecoli_bwaes_origin.bam  /Sid/tdvorkina/gralign/E.coli_synth/reads/pacbio/sd_0001.fasta $FILE_ECOLI
echo "--------------Celegans--------------"
python2 ~/scripts/aligner_stats/aligner_stats.py ~/gralign/bwamem_aligned/celegans/celegans_bwaes_origin.bam /Sid/tdvorkina/gralign/C.elegans/reads/pacbio/celegans.pacbio.sub10k.fasta  $FILE_CELEGANS
echo "--------------Rumen--------------"
python2 ~/scripts/aligner_stats/aligner_stats_noref.py /Sid/tdvorkina/gralign/Rumen/reads/rumen.nanopore.fasta $FILE_RUMEN


echo "--------------Ecoli--------------"
python2 ~/scripts/aligner_stats/count_graphalnscore.py ~/gralign/bwamem_aligned/ecoli/$PREFIX $FILE_ECOLI ~/gralign/bwamem_aligned/ecoli/ecoli.fasta /Sid/tdvorkina/gralign/E.coli_synth/reads/pacbio/sd_0001.fasta ~/gralign/bwamem_aligned/ecoli/ecoli_bwaes_origin.bam
echo "--------------Celegans--------------"
python2 ~/scripts/aligner_stats/count_graphalnscore.py ~/gralign/bwamem_aligned/celegans/$PREFIX $FILE_CELEGANS ~/gralign/bwamem_aligned/celegans/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa /Sid/tdvorkina/gralign/C.elegans/reads/pacbio/celegans.pacbio.sub10k.fasta ~/gralign/bwamem_aligned/celegans/celegans_bwaes_origin.bam

