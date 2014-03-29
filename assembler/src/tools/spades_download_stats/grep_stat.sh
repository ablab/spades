#!/bin/bash

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


rm ./spades$1$2$3.txt

grep SPAdes-$1.$2.$3 /var/log/apache2/spades-access.log >> ./spades$1$2$3.txt
grep SPAdes-$1.$2.$3 /var/log/apache2/spades-access.log.1 >> ./spades$1$2$3.txt

for i in $(seq 2 1 50)
do
	echo $i
	zgrep SPAdes-$1.$2.$3 /var/log/apache2/spades-access.log.$i.gz >> ./spades$1$2$3.txt
done

