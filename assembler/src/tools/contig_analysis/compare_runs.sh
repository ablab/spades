#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


for file in $(ls $1/saves/*.sqn)
do
	echo $(basename $file)
	python compare_fasta.py $1/saves/$(basename $file)  $2/saves/$(basename $file)
done
