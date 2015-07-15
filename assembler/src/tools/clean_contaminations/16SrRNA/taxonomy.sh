#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


if [ -z "$1" ]
then
	echo Usage: $0 reads.fastq [threads=1]
	echo Example: $0 my_reads.fastq 32
	exit
fi

if [ -z "$2" ]
then
	threads="1"
else
	threads=$2
fi

echo Analyzing taxonomy for $1 ...

bowtie --threads $threads -k 10 --best --chunkmbs 1024 SSURef_108_NR_tax_silva_v2.fasta $1 | python parse_tax_from_bowtie.py
