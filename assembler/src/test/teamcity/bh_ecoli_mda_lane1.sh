#!/bin/bash
set -e
pushd ../../../
rm -rf spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH
./spades.py src/test/teamcity/bh_config_ecoli_mda_lane1.info
python src/tools/reads_utils/reads_quality.py -r data/input/E.coli/MG1655-K12.fasta.gz -o spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/reads_quality_results spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/corrected/ECOLI_MDA_LANE1_TEAMCITY_BH.dataset
python src/test/teamcity/assess_reads.py spades_output/ECOLI_MDA_LANE1_TEAMCITY_BH/reads_quality_results/all.tsv 90 97 # 97 97
popd
