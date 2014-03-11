#!/bin/bash

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


rm ./spades2r.txt

grep spades.rpm /var/log/apache2/spades-access.log >> ./spades2r.txt
grep spades.rpm /var/log/apache2/spades-access.log.1 >> ./spades2r.txt


for i in $(seq 2 1 52)
do
        echo $i
        zgrep spades.rpm /var/log/apache2/spades-access.log.$i.gz >> ./spades2r.txt
done

