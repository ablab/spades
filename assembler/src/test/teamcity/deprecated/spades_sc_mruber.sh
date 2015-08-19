#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/MRUBER

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/M.ruber/Mru_no9.fastq -o /tmp/data/output/spades_output/MRUBER

rm -rf ~/quast-1.3/MRUBER/

python ~/quast-1.3/quast.py -R data/input/M.ruber/ref.fasta -G data/input/M.ruber/MRU_genes.txt -o ~/quast-1.3/MRUBER/ /tmp/data/output/spades_output/MRUBER/contigs.fasta

python src/test/teamcity/assess.py ~/quast-1.3/MRUBER/transposed_report.tsv 20000 20 1500 50.0 30 10
exitlvl=$?
popd

