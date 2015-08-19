#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################



while read -r LINE; do 
    echo $LINE;
    LINE1=$(echo "$LINE" | sed 's/\([0-9]*\).*/\1/' );
    LINE2=$(echo "$LINE" | sed 's/[0-9]* \([0-9]*\).*/\1/' );
    echo "$LINE1 $LINE2";
    echo `grep -A1 $LINE2 distance_estimation.sqn` > "$LINE1""_""$LINE2.fasta"
    echo `grep -A1 $LINE1 distance_estimation.sqn` >> "$LINE1""_""$LINE2.fasta"
    sed -i 's/\s/\n/' "$LINE1""_""$LINE2.fasta"
    nucmer --maxmatch --coords -c 55 /acestorage/data/input/E.coli/MG1655-K12.fasta "$LINE1""_""$LINE2.fasta"
    mv out.coords "$LINE1""_""$LINE2.coords"

done < $1;
