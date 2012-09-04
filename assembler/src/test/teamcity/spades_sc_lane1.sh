#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf spades_output/ECOLI_SC_LANE_1_BH_woHUMAN
rm -rf ~/quast/ECOLI_SC_LANE_1_BH_woHUMAN/
python ~/quast/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast/ECOLI$
python src/test/teamcity/assess.py ~/quast/ECOLI_SC_LANE_1_BH_woHUMAN/transposed_report.tsv 85000 6
popd
