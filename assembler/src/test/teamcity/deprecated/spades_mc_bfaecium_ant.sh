#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf /tmp/data/output/spades_output/BFAECIUM_QUAKE
./spades.py -k 21,33,55 -m 16 --only-assembler -1 data/input/B.faecium/std_left.cor.fastq.gz -2 data/input/B.faecium/std_right.cor.fastq.gz -s data/input/B.faecium/std_right.cor_single.fastq.gz -s data/input/B.faecium/std_left.cor_single.fastq.gz -o /tmp/data/output/spades_output/BFAECIUM_QUAKE
rm -rf ~/quast-1.3/BFAECIUM_QUAKE/
python ~/quast-1.3/quast.py -R /tmp/data/input/B.faecium/ref.fasta.gz -o ~/quast-1.3/BFAECIUM_QUAKE/ /tmp/data/output/spades_output/BFAECIUM_QUAKE/contigs.fasta
python src/test/teamcity/assess.py ~/quast-1.3/BFAECIUM_QUAKE/transposed_report.tsv 330000 0
popd
