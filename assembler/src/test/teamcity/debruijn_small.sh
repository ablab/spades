#!/bin/bash
./debruijn_make.sh
pushd ../../../
./cpcfg
./run rd
#cd src/tools/quality/
#./run_EColi_morality.sh -m ../../../data/debruijn/QUAKE_CROPPED_1K/K55/latest/resolved_and_cleared_*.fasta
#cd ../../test/teamcity/

#python teamcity.py
