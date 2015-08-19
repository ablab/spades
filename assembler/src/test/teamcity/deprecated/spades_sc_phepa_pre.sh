#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/PHEPARINUS_PRENORM

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/P.heparinus/prenoemalized/Phe_no7.subsample.contam.artifact.clean.fastq -o /tmp/data/output/spades_output/PHEPARINUS_PRENORM

rm -rf ~/quast-1.3/PHEPARINUS_PRENORM/

python ~/quast-1.3/quast.py -R data/input/P.heparinus/ref.fasta -G data/input/P.heparinus/PHE_genes.txt -o ~/quast-1.3/PHEPARINUS_PRENORM/ /tmp/data/output/spades_output/PHEPARINUS_PRENORM/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/PHEPARINUS_PRENORM/transposed_report.tsv 160000 9 3900 95.5 10.0 1.0
exitlvl=$?
popd

