#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/ECOLI_SC_JGI_LANE_7

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/E.coli/jgi_lane_7/Eco_no7.fastq -o /tmp/data/output/spades_output/ECOLI_SC_JGI_LANE_7

rm -rf ~/quast-1.3/ECOLI_SC_JGI_LANE_7/

python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_SC_JGI_LANE_7/ /tmp/data/output/spades_output/ECOLI_SC_JGI_LANE_7/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_SC_JGI_LANE_7/transposed_report.tsv 115000 4 4260 99.8 3.2 0.7
exitlvl=$?
popd

