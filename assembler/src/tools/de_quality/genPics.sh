#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

i=0

if [ "x$3" != "x" ]
then
    i=$3
    i=`expr $i - 1`
fi
STRR=`wc -l $1`
SIZZE=`echo $STRR | cut --delimiter=' ' -f 1`
PROJ_DIR=`pwd`;
OUTPUT_DIR=$PROJ_DIR/data/output/

mkdir -p $OUTPUT_DIR/pics
while [ $i -lt $SIZZE ]; do
    i=`expr $i + 1`;
    echo "Iteration "$i" of "$SIZZE
    LINE=`head -$i $1 | tail -1`
    LINE1=$(echo "$LINE" | sed 's/\([0-9]*\).*/\1/' );
    LINE2=$(echo "$LINE" | sed 's/[0-9]* \([0-9]*\).*/\1/' );
    echo "$LINE1 $LINE2";
    f=`find $OUTPUT_DIR/hists/ -name $LINE1"_"$LINE2*`
    echo "Dir "$f
    if [ "x$2" != "x" ]
    then
        grep -P $LINE1'.*'$LINE2'.*:' $2
    fi
    cd $f
    pwd
    gnuplot png_plot.conf > /dev/null 2> /dev/null || { echo "Gnuplot exited "; exit 1; }
    name=`basename $f`
    cd $OUTPUT_DIR
    mv $f/$name.png pics/$i'_'$name.png    
    ls -hl pics | wc -l
    cd $PROJ_DIR
done;
