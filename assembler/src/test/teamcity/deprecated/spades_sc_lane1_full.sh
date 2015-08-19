#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/ECOLI_SC_LANE_1

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/E.coli/sc_lane_1/ecoli_mda_lane1.fastq.gz -o /tmp/data/output/spades_output/ECOLI_SC_LANE_1

rm -rf ~/quast-1.3/ECOLI_SC_LANE_1/

python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_SC_LANE_1/ /tmp/data/output/spades_output/ECOLI_SC_LANE_1/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_SC_LANE_1/transposed_report.tsv 108000 2 4080 96.4 4 0.7
exitlvl=$?
popd

