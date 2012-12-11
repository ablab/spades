#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/BFAECIUM_QUAKE
./spades.py -k 21,33,55 -m 16 --only-assembler -1 /tmp/data/input/B.faecium/std_left.cor.fastq.gz -2 /tmp/data/input/B.faecium/std_right.cor.fastq.gz -s /tmp/data/input/B.faecium/std_right.cor_single.fastq.gz -s /tmp/data/input/B.faecium/std_left.cor_single.fastq.gz -o spades_output/BFAECIUM_QUAKE
rm -rf ~/quast-1.3/BFAECIUM_QUAKE/
python ~/quast-1.3/quast.py -R /tmp/data/input/B.faecium/ref.fasta.gz -o ~/quast-1.3/BFAECIUM_QUAKE/ spades_output/BFAECIUM_QUAKE/contigs.fasta
#python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_IS220_QUAKE_1K/transposed_report.tsv 1000 0
popd
