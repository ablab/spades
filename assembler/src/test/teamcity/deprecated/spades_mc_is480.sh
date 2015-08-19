
#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf /tmp/data/output/spades_output/ECOLI_IS480_QUAKE
./spades.py --only-assembler -k 21,33,55 -m 12 -1 data/input/E.coli/is480/ERR022075_1.cor.fastq.gz -2 data/input/E.coli/is480/ERR022075_2.cor.fastq.gz -s data/input/E.coli/is480/ERR022075_1.cor_single.fastq.gz -s data/input/E.coli/is480/ERR022075_2.cor_single.fastq.gz -o /tmp/data/output/spades_output/ECOLI_IS480_QUAKE
rm -rf ~/quast-1.3/ECOLI_IS480_QUAKE/
python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_IS480_QUAKE/  /tmp/data/output/spades_output/ECOLI_IS480_QUAKE/contigs.fasta
python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_IS480_QUAKE/transposed_report.tsv 133000 4
popd
