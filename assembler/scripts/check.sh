#!/bin/bash
cut -d ' ' -f -4 $1 > temp.prd;
while read -r line; do
#echo $line;
#echo $(cat saves/2_repeats_resolved_before.prd | grep "$line" | wc -l);
line1=$(echo $line | cut -d ' ' -f -2);
NUM=$(cat "saves/repeats_resolved_before.prd" | grep "$line1" | wc -l)
if [ $NUM -ge 2 ]; 
then  
    echo "FOUND" $line ", the number of points is " $NUM;
fi;
done < temp.prd;
rm temp.prd;
