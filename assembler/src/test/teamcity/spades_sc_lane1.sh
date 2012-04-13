#!/bin/bash
set -e
./prepare.sh
pushd ../../../
rm -f spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/latest
./spades.py src/test/teamcity/spades_config_sc_lane1.info
python src/test/teamcity/assess.py spades_output/ECOLI_SC_LANE_1_BH_woHUMAN/latest/quality_results/all.tsv 47000 1
popd
