#!/bin/bash
pushd ../../../
rm data/debruijn/ECOLI_IS220_QUAKE/K55/latest
./prepare_cfg
popd
./debruijn_make.sh
pushd ../../../
./cpcfg
cp ./src/test/teamcity/config.info.FULL ./src/debruijn/config.info
./run rd
cd src/tools/quality/
./run_EColi.sh -m ../../../data/debruijn/ECOLI_IS220_QUAKE/K55/latest/resolved_and_cleared_*.fasta
popd

python teamcity.py
