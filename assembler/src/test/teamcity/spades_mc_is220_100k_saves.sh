#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -f spades_output/ECOLI_IS220_QUAKE_100K_SAVES/latest
#etalon=spades_output/ECOLI_IS220_QUAKE_100K_SAVES/etalon
etalon=../../../etalon_output/ECOLI_IS220_QUAKE_100K_SAVES/etalon
if [ ! -d $etalon ]; then
    echo "Error: no etalon saves at $etalon"
    exit 9
fi
./spades.py src/test/teamcity/spades_config_mc_is220_100k_saves.info

pushd spades_output/ECOLI_IS220_QUAKE_100K_SAVES/latest
diffs=0
for i in link_*
do
    echo "Exploring diff of saves for K=$i"
    for f in $i/saves/*
    do
        set +e
        diff $f ../../../$etalon/$f >> diff_with_etalon.txt
        errlvl=$?
        if [ $errlvl -ne 0 ]; then
            if [ $errlvl -eq 1 ]; then
                echo "^^^^^^^ it was $f" >> diff_with_etalon.txt
                echo "BAD: difference found in $f"
            elif
                echo "BAD: unable to compare with $f"
            fi
            (( diffs += 1 ))
        fi
        set -e
    done
done
popd
popd

echo $diffs differences with etalon saves found
exit $diffs
