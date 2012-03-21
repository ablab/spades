#!/bin/bash
set -e
pushd ../../../
rm -f data/debruijn/ECOLI_SC_LANE_1_BH_woHUMAN/K55/latest
rm -rf data/quality
./cpcfg
./spades.py src/test/teamcity/spades_config_sc_lane1.info
src/tools/quality/run_Ecoli.sh -o data/quality data/debruijn/ECOLI_SC_LANE_1_BH_woHUMAN/K55/latest/final_contigs.fasta
python src/test/teamcity/assess.py data/quality/all.tsv 47000 1
popd
