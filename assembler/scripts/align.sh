#!/bin/bash


while read -r LINE; do 
    echo $LINE;
    LINE1=$(echo "$LINE" | sed 's/[A-Z]* \([0-9]*\).*/\1/' );
    LINE2=$(echo "$LINE" | sed 's/[A-Z]* [0-9]* \([0-9]*\).*/\1/' );
    echo "$LINE1 $LINE2";
    \> > "$LINE1""_""$LINE2.fasta";
    cut -d ' ' -f -2 2_simplified_graph.sqn | sed  \/$LINE1'/!d' >>"$LINE1""_""$LINE2.fasta"
    \> > "$LINE1""_""$LINE2.fasta";
    cut -d ' ' -f -2 2_simplified_graph.sqn | sed  \/$LINE2'/!d' >> "$LINE1""_""$LINE2.fasta"
    sed -i 's/\s/\n/' "$LINE1""_""$LINE2.fasta"
    nucmer --maxmatch --coords /home/ftp/data/input/E.Coli.K12.MG1655/MG1655-K12.fasta "$LINE1""_""$LINE2.fasta"
    cp out.coords "$LINE1""_""$LINE2.coords"

    break;
done < $1;
