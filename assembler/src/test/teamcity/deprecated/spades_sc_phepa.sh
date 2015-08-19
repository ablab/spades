#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../

rm -rf /tmp/data/output/spades_output/PHEPARINUS

#./spades_compile.sh
./spades.py --sc -m 160 -k 21,33,55 --12 data/input/P.heparinus/Phe_no7.fastq.gz -o /tmp/data/output/spades_output/PHEPARINUS

rm -rf ~/quast-1.3/PHEPARINUS/

python ~/quast-1.3/quast.py -R data/input/P.heparinus/ref.fasta -G data/input/P.heparinus/PHE_genes.txt -o ~/quast-1.3/PHEPARINUS/ /tmp/data/output/spades_output/PHEPARINUS/contigs.fasta
python src/test/teamcity/assess.py ~/quast-1.3/PHEPARINUS/transposed_report.tsv 140000 6 4140 98.0 8.0 2.0
exitlvl=$?
popd

