#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13
do
	echo Copying to ant$i...
	srun -w ant$i ./copy_to_ant.sh $i
	echo Done
done
