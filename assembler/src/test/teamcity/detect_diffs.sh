#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
	
if [ "$#" -ne 2 ]; then
    echo "Usage: detect_diffs.sh <target folder> <etalons folder>"
    exit
fi

target_folder=$1
etalons_folder=$2
if [ ! -d $etalons_folder ]; then
    echo "Error: no etalon saves at $etalons_folder"
    exit 9
fi

diffs=0
    for f in $etalons_folder/*/saves/* $etalons_folder/*/saves/*/* $etalons_folder/before_rr.fasta $etalons_folder/assembly_graph* $etalons_folder/scaffolds* $etalons_folder/contigs*
    do
        if [[ -f $f ]]; then
            echo "Checking diffs in " $f
            set +e
            diff $f $target_folder/${f#$etalons_folder} >> diff_with_etalon.txt
            errlvl=$?
            if [ $errlvl -ne 0 ]; then
                if [ $errlvl -eq 1 ]; then
                    echo "^^^^^^^ it was $f" >> diff_with_etalon.txt
                    echo "BAD: difference found in $f"
                else
                    echo "BAD: unable to compare with $f"
                fi
                (( diffs += 1 ))
            fi
            set -e
        fi
    done

echo $diffs differences with etalon saves found
#returning $diffs is bad idea since return code should be less than 255
if [ $diffs -ne 0 ]; then
	exit 1
else
	exit 0
fi
