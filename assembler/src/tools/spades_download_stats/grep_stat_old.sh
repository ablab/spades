#!/bin/bash

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


rm ./spades2$1.txt

grep spades-2.$1.0 /var/log/apache2/spades-access.log >> ./spades2$1.txt
grep spades-2.$1.0 /var/log/apache2/spades-access.log.1 >> ./spades2$1.txt

for i in $(seq 2 1 52)
do
	echo $i
	zgrep spades-2.$1.0 /var/log/apache2/spades-access.log.$i.gz >> ./spades2$1.txt
done

