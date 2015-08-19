#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/ECOLI_SC_LANE_1_BH_woHUMAN

#./spades_compile.sh
./spades.py --sc --debug --only-assembler -m 160 -k 21,33,55 -1 data/input/E.coli/sc_lane_1/bh20111014/human_filtered_paired_left.fastq.gz -2 data/input/E.coli/sc_lane_1/bh20111014/human_filtered_paired_right.fastq.gz -s data/input/E.coli/sc_lane_1/bh20111014/human_filtered_single_left.fastq.gz -s data/input/E.coli/sc_lane_1/bh20111014/human_filtered_single_right.fastq.gz -o /tmp/data/output/spades_output/ECOLI_SC_LANE_1_BH_woHUMAN

rm -rf ~/quast-1.3/ECOLI_SC_LANE_1_BH_woHUMAN/

python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_SC_LANE_1_BH_woHUMAN/ /tmp/data/output/spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_SC_LANE_1_BH_woHUMAN/transposed_report.tsv 85000 7
exitlvl=$?
popd

if [ $exitlvl -ne 0 ]; then
    exit $exitlvl
else
    etalon=/smallnas/teamcity/etalon_output/ECOLI_SC_LANE_1_BH_woHUMAN/etalon
    ./detect_diffs.sh /tmp/data/output/spades_output/ECOLI_SC_LANE_1_BH_woHUMAN $etalon
    exit $?
fi
