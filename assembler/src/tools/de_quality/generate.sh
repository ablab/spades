#!/bin/bash
i=0

if [ "x$3" != "x" ]
then
    i=$3
    i=$i-1
fi
STRR=`wc -l $1`
SIZZE=`echo $STRR | cut --delimiter=' ' -f 1`
DIRR=`pwd`;
mkdir overall
while [ $i -lt $SIZZE ]; do
    i=$[i+1];
    echo "Iteration "$i" of "$SIZZE
    LINE=`head -$i $1 | tail -1`
    LINE1=$(echo "$LINE" | sed 's/\([0-9]*\).*/\1/' );
    LINE2=$(echo "$LINE" | sed 's/[0-9]* \([0-9]*\).*/\1/' );
    echo "$LINE1 $LINE2";
    f=`find data/out/ -name $LINE1"_"$LINE2*`
    echo "Dir "$f
    grep -P $LINE1'.*'$LINE2'.*:' $2
    cd $f
    pwd
    gnuplot png_plot.conf > /dev/null 2> /dev/null || { echo "Gnuplot exited "; exit 1; }
    name=`basename $f`
    mv ../../../overall/$name.png ../../../overall/$i'_'$name.png    
    ls -hl ../../../overall | wc -l
    cd $DIRR
done;
