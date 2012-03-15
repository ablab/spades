#!/bin/bash
set -e
pushd ../../../
rm -f data/debruijn/ECOLI_IS220_QUAKE/K55/latest
rm -rf data/quality
./prepare_cfg
make clean
./cpcfg
./spades.py src/test/teamcity/spades_config_rectangles_is220.info
src/tools/quality/run_Ecoli.sh -o data/quality data/debruijn/ECOLI_IS220_QUAKE/K55/latest/saves/rectangle_after.fasta
python src/test/teamcity/assess.py data/quality/all.tsv 80000 3
popd
