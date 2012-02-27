#!/bin/bash

rootDir=`pwd`

for dir in $@; do
    cd $rootDir/$dir
    reads_file=`ls *.fasta`
    echo "Stats for $dir/$reads_file ##############################################"

# reads
    echo -ne "Number of reads:\t"
    grep -c '^>' $reads_file

# vertices and edges
    grep ^Number_of_ transformed/$reads_file.bgraph | perl -lane 'print $F[0],":\t",$F[1];'
    
# components
    echo -ne "Number of components:\t"
    $EUSRC/assembly/x86_64/printGraphSummary  transformed/$reads_file -components | wc -l
    
# contigs
    ./fasta2rows.pl $reads_file.contig | perl -lane '@line = split /\t/; print length($line[1]);' | ./moments.pl

# retained reads
    echo -ne "Num. retained reads:\t";
    $EUSRC/assembly_utils/PrintRetainedReads.pl $reads_file transformed/$reads_file.path | grep -c '^>'

# coverage
    $EUSRC/assembly_utils/IntervalsToCoverage.pl transformed/$reads_file.intv | tail -n 1

    echo; echo;
done;