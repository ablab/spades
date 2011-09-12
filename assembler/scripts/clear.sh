#!/bin/bash


sort -rnk 12,12 $1 > $1.temp
while read -r LINE; do 
#    echo $LINE;
    NUM1=$(echo "$LINE" | sed 's/[A-Z]* \([0-9]*\).*/\1/' );
    NUM2=$(echo "$LINE" | sed 's/[A-Z]* [0-9]* \([0-9]*\).*/\1/' );
    NUM3=$(echo "$LINE" | sed 's/[A-Z]* [0-9]* [0-9]* \(-\).*/\1/' );
#    echo "$NUM1 $NUM2";
#    echo $NUM3;
    let "rem = $NUM1 % 2"
    if [ $rem -eq 0 -a "x$NUM3" != "x-" ]; then
        echo $LINE ;           
    fi
done < $1.temp;
rm $1.temp
