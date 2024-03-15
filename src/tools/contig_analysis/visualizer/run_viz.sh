#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


rm list.txt
touch list.txt

for file in $(ls *.stdout)
do
	python parser.py $file >> list.txt
done

python misassemblies.py -o plot -a list.txt --arcs

rm list.txt
