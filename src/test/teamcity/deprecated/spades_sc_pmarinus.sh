#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/PMARINUS

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 -1 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R1_001.fastq.gz -1 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R1_002.fastq.gz -1 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R1_003.fastq.gz -1 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R1_004.fastq.gz -2 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R2_001.fastq.gz -2 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R2_002.fastq.gz -2 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R2_003.fastq.gz -2 data/input/P.marinus/AC-193-C02/lane6-index21-GTTTCG-AC-193-C02_GTTTCG_L006_R2_004.fastq.gz -o /tmp/data/output/spades_output/PMARINUS

rm -rf ~/quast-1.3/PMARINUS/

python ~/quast-1.3/quast.py -R data/input/P.marinus/ref.fasta -o ~/quast-1.3/PMARINUS/ /tmp/data/output/spades_output/PMARINUS/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/PMARINUS/transposed_report.tsv 292000 0 1612 95.3 3.4 1.5
exitlvl=$?
popd

