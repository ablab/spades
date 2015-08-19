#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# input_folder references_folder output_folder threads

mkdir $3
mkdir $3/tmp
mkdir $3/all_reports
ls $1 | xargs -P $4 -I {} -t quast --min-contig 1000 --contig-thresholds 5000,8000,12000 -e -R $2/{}.fasta $1/{}/scaffolds.fasta -o $3/tmp/{}
ls $1 | xargs -P $4 -I {} -t cp $3/tmp/{}/report.tsv $3/all_reports/{}.tsv
python quast_all.py $3/all_reports $3
grep "( inversion )" $3/tmp/*/contigs_reports/contigs_report_*.stdout > $3/inversions.txt
grep "( translocation )" $3/tmp/*/contigs_reports/contigs_report_*.stdout > $3/translocations.txt
grep "( relocation," $3/tmp/*/contigs_reports/contigs_report_*.stdout > $3/relocations.txt
wc -l $3/*.txt > $3/mis_classification.info
