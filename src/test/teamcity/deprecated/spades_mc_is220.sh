
#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
rm -rf /tmp/data/output/spades_output/ECOLI_IS220_QUAKE
./spades.py --generate-sam-file --only-assembler -m 12 -k 21,33,55 -1 data/input/E.coli/is220/s_6_1.cor.fastq.gz -2 data/input/E.coli/is220/s_6_2.cor.fastq.gz -s data/input/E.coli/is220/s_6_1.cor_single.fastq.gz -s data/input/E.coli/is220/s_6_2.cor_single.fastq.gz -o /tmp/data/output/spades_output/ECOLI_IS220_QUAKE
rm -rf ~/quast-1.3/ECOLI_IS220_QUAKE/
python ~/quast-1.3/quast.py -R data/input/E.coli/MG1655-K12.fasta.gz -G data/input/E.coli/genes/genes.gff -O data/input/E.coli/genes/operons.gff -o ~/quast-1.3/ECOLI_IS220_QUAKE/  /tmp/data/output/spades_output/ECOLI_IS220_QUAKE/contigs.fasta
python src/test/teamcity/assess.py ~/quast-1.3/ECOLI_IS220_QUAKE/transposed_report.tsv 80000 3
popd
