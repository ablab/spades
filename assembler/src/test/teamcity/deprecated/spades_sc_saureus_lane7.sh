#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/SAUREUS_LANE_7

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 -1 data/input/S.aureus/sc_lane_7/bacteria_mda_lane7_left.fastq -2 data/input/S.aureus/sc_lane_7/bacteria_mda_lane7_right.fastq -o /tmp/data/output/spades_output/SAUREUS_LANE_7

rm -rf ~/quast-1.3/SAUREUS_LANE_7/

python ~/quast-1.3/quast.py -R data/input/S.aureus/USA300_FPR3757.fasta -G data/input/S.aureus/genes/bacteria_genes.txt -o ~/quast-1.3/SAUREUS_LANE_7/ /tmp/data/output/spades_output/SAUREUS_LANE_7/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/SAUREUS_LANE_7/transposed_report.tsv 85000 7 2500 99.8 6 6
exitlvl=$?
popd

