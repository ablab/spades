#!/bin/bash
set -e
pushd ../../../
rm -f data/debruijn/ECOLI_IS220_QUAKE/latest
./cpcfg
./spades.py src/test/teamcity/spades_config_mc_is220.info
python src/test/teamcity/assess.py data/debruijn/ECOLI_IS220_QUAKE/latest/quality_results/all.tsv 80000 3
popd
