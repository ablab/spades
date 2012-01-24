#!/bin/bash
set -e
pushd ../../../
[ -e data/debruijn/ECOLI_IS220_QUAKE/K55/latest ] && rm data/debruijn/ECOLI_IS220_QUAKE/K55/latest
[ -e data/quality ] && rm -rf data/quality
make clean
./cpcfg
./spades.py src/test/teamcity/spades_config.full.info
src/tools/quality/run_Ecoli.sh -o data/quality data/debruijn/ECOLI_IS220_QUAKE/K55/latest/final_contigs.fasta
python src/test/teamcity/assess.py data/quality/all.tab
