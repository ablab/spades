#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/ECOLI_1K_TEAMCITY_BH
./spades.py --only-error-correction -1 ./data/input/E.coli/is220/cropped/s_6.first1000_1.fastq.gz -2 ./data/input/E.coli/is220/cropped/s_6.first1000_2.fastq.gz -o spades_output/ECOLI_1K_TEAMCITY_BH
python src/tools/reads_utils/reads_quality.py -r test_dataset/reference_1K.fa.gz -o spades_output/ECOLI_1K_TEAMCITY_BH/reads_quality_results spades_output/ECOLI_1K_TEAMCITY_BH/corrected/dataset.info
python src/test/teamcity/assess_reads.py spades_output/ECOLI_1K_TEAMCITY_BH/reads_quality_results/report.horizontal.tsv 100 100
popd
