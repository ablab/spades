#!/bin/bash

if [ $# -eq 0 ] 
then
	echo "Usage $0 forward_read.fa reverse_reads.fa outfile.fa"
	echo "\tforward_reads.fa / reverse_reads.fa : paired reads to be merged"
	echo "\toutfile.fa : outfile to be created"
	exit
fi

if [ ! -e $1 ] 
then
	echo "$1 does not exist"
	exit
fi
 
if [ ! -e $2 ] 
then
	echo "$2 does not exist"
	exit
fi

exec 4<$1
exec 5<$2
exec >$3
 
read lineA <&4
read lineB <&5

while [ $lineA ] 
do
	echo $lineA 
	read lineA <&4
	while [[ -n "$lineA"  && ! "$lineA" =~ \>* ]] 
	do
		echo $lineA 
		read lineA <&4
	done

	echo $lineA 
	read lineA <&4
	while [[ -n "$lineA"  && ! "$lineA" =~ \>* ]] 
	do
		echo $lineA 
		read lineA <&4
	done

	echo $lineB 
	read lineB <&5
	while [[ -n "$lineB"  && ! "$lineB" =~ \>* ]] 
	do
		echo $lineB 
		read lineB <&5
	done

	echo $lineB 
	read lineB <&5
	while [[ -n "$lineB"  && ! "$lineB" =~ \>* ]] 
	do
		echo $lineB 
		read lineB <&5
	done
done
