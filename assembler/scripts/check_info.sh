#!/bin/bash

cp 2_simplified_graph.prd simple.prd
cp a_simplified_graph.prd advanced.prd
cut -d ' ' -f -3 2_simplified_graph_et.prd > etalon.prd
cut -d ' ' -f -4 2_simplified_graph.prd > simple.prd
cut -d ' ' -f -4 a_simplified_graph.prd > advanced.prd
sed '1d' etalon.prd > etalon1.prd
sed '/^$/d' etalon1.prd > etalon.prd

while read -r LINE; do 
    echo $LINE;
    sed -i /"$LINE"'/d' simple.prd; 
    sed -i /"$LINE"'/d' advanced.prd;
done < etalon.prd;
