#!/bin/bash
set -e
pushd ../../../
[ -e data/debruijn/ECOLI_IS220_QUAKE/K55/latest ] && rm data/debruijn/ECOLI_IS220_QUAKE/K55/latest
[ -e src/tools/quality/results ] && rm -rf src/tools/quality/results
make clean
./cpcfg
./spades.py src/test/teamcity/spades_config.full.info
src/tools/quality/run_Ecoli.sh -m data/debruijn/ECOLI_IS220_QUAKE/K55/latest/resolved_and_cleared_*.fasta
popd
python assess.py
