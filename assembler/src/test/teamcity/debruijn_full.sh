#!/bin/bash
. ../../../prepare_cfg
./debruijn_make.sh
cd ../../../
./cpcfg
cp ./src/test/teamcity/config.info.FULL ./src/debruijn/config.info
./run rd
cd src/tools/quality/
./run_EColi_morality.sh -m ../../../data/debruijn/QUAKE_FULL/K55/latest/resolved_and_cleared_*.fasta
cd ../../test/teamcity/

python teamcity.py
