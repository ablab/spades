#!/bin/bash
set -e
pushd ../../../
rm -f data/debruijn/ECOLI_IS220_QUAKE/K55/latest
rm -rf data/quality
./prepare_cfg
if [ ! -e "./data/input" ]
then
  pushd data
  ./link_morality.sh
  popd
fi
make clean
./cpcfg
sed -i -r 's/^resolving_mode[ \t]*split/resolving_mode rectangle/g' configs/debruijn/config.info
sed -i -r 's/ECOLI_IS220_QUAKE_1K/ECOLI_IS220_QUAKE/g' configs/debruijn/config.info
make rd
./run
src/tools/quality/run_Ecoli.sh -o data/quality data/debruijn/ECOLI_IS220_QUAKE/K55/latest/saves/rectangle_after.fasta
python src/test/teamcity/assess.py data/quality/all.tsv 80000 3
popd
