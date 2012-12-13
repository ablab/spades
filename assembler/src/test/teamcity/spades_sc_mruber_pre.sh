#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/MRUBER_PRENORM

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/M.ruber/prenoemalized/Mru_no9.subsample.contam.artifact.clean.fastq -o /tmp/data/output/spades_output/MRUBER_PRENORM

rm -rf ~/quast-1.3/MRUBER_PRENORM/

python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/MRUBER_PRENORM/ /tmp/data/output/spades_output/MRUBER_PRENORM/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/MRUBER_PRENORM/transposed_report.tsv 15100 24 2148 81.0 100 100
exitlvl=$?
popd

