#!/bin/bash

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


for file in $(ls *.stdout)
do
	python parser.py $file >> list.txt
done
