#!/bin/bash

#cp 2_simplified_graph.prd simple.prd;
#cp a_simplified_graph.prd advanced.prd;
cut -d ' ' -f -3 2_simplified_graph_et.prd > etalon.prd;
cut -d ' ' -f -3 $1 > temp.prd;

cp etalon.prd $2;
echo $1 $2;
while read -r LINE; do 
    echo $LINE;
    sed -i '/'"$LINE"'/d' $2 
done < temp.prd;
