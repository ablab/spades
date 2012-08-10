#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -f spades_output/ECOLI_SC_BH_SPADES
read line < /storage/data/input/E.coli.K12/mda_lane_1/spades.options 
echo $line

python ./spades.py $line -o spades_output/ECOLI_SC_BH_SPADES
read line < /storage/data/input/E.coli.K12/mda_lane_1/quast.options 
echo $line
python ./smallnas/dima/quast/quast.py $line -o spades_output/ECOLI_SC_BH_SPADES/quality_results/
python src/test/teamcity/assess.py spades_output/ECOLI_SC_BH_SPADES/quality_results/transposed_report.tsv 105000 3
popd
