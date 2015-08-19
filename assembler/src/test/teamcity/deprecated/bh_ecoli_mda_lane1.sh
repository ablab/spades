#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf /tmp/data/output/spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH
./spades.py --sc --only-error-correction --12 ./data/input/E.coli/sc_lane_1/ecoli_mda_lane1.fastq.gz -m 160 -o /tmp/data/output/spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH
python src/tools/reads_utils/reads_quality.py -r data/input/E.coli/MG1655-K12.fasta.gz -o /tmp/data/output/spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/reads_quality_results /tmp/data/output/spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/corrected/dataset.info
python src/test/teamcity/assess_reads.py /tmp/data/output/spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/reads_quality_results/report.horizontal.tsv 90 97 # 97 97
popd
