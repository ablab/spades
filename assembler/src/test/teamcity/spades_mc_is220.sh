#!/bin/bash
set -e
pushd ../../../
rm -f data/debruijn/ECOLI_IS220_QUAKE/K55/latest
rm -rf data/quality
./cpcfg
./spades.py src/test/teamcity/spades_config_mc_is220.info
python src/test/teamcity/assess.py data/quality/all.tsv 80000 3
popd
