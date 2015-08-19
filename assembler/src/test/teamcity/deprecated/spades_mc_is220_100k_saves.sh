#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

target=spades_output/ECOLI_IS220_QUAKE_100K_SAVES
if [ -d $target ]; then
	rm -rf $target #etalon=spades_output/ECOLI_IS220_QUAKE_100K_SAVES/etalon
fi	

etalon=/smallnas/teamcity/etalon_output/ECOLI_IS220_QUAKE_100K_SAVES/etalon_saves

#if [ ! -d $etalon ]; then
#    echo "Error: no etalon saves at $etalon"
#    exit 9
#fi

#./spades_compile.sh
./spades.py --only-assembler --debug -t 4 -k 55 -1 data/input/E.coli/is220/cropped/s_6.first100000_1.fastq.gz -2 data/input/E.coli/is220/cropped/s_6.first100000_2.fastq.gz -o $target

#pushd spades_output/ECOLI_IS220_QUAKE_100K_SAVES
#diffs=0
#    for f in saves/*
#    do
#        set +e
#        diff $f $etalon/link_K55/$f >> diff_with_etalon.txt
#        errlvl=$?
#        if [ $errlvl -ne 0 ]; then
#            if [ $errlvl -eq 1 ]; then
#                echo "^^^^^^^ it was $f" >> diff_with_etalon.txt
#                echo "BAD: difference found in $f"
#            else
#                echo "BAD: unable to compare with $f"
#            fi
#            (( diffs += 1 ))
#        fi
#        set -e
#    done
#popd
popd

#echo $diffs differences with etalon saves found
./detect_diffs.sh ../../../$target $etalon
exit $?
