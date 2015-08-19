#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/ECOLI_IS220_QUAKE_1K
./spades.py -k 21,55 --only-assembler -1 /tmp/data/input/E.coli/is220/cropped/s_6.first1000_1.fastq.gz -2 /tmp/data/input/E.coli/is220/cropped/s_6.first1000_2.fastq.gz -o spades_output/ECOLI_IS220_QUAKE_1K
rm -rf ~/quast-1.3/ECOLI_IS220_QUAKE_1K/
python ~/quast-1.3/quast.py -R /tmp/data/input/E.coli/MG1655-K12.fasta.gz -G /tmp/data/input/E.coli/genes/genes.gff -O /tmp/data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_IS220_QUAKE_1K/ spades_output/ECOLI_IS220_QUAKE_1K/contigs.fasta
python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_IS220_QUAKE_1K/transposed_report.tsv 1000 0
popd
