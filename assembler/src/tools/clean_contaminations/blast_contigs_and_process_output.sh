#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Blasts contigs against blast database and processes the output: finds the sequences to which the contigs match most often (100 most often occurring sequences). The script considers 30 best matches of each contig to the blast database.

function usage ()
{
    echo "usage: $1 -c <contigs.fasta> [-o <output folder>] [-d <blast database>]"
}


##### Main

contigs=
output_folder="."
blast_database="/acestorage/blast/db/nt"

script=$0
scriptname=${script##*/}
path_to_script=${script%$scriptname}

while [ "$1" != "" ]; do
    case $1 in
        -c | --contigs )           shift
                                   contigs=$1
                                   ;;
        -o | --output )            shift
                                   output_folder=$1
                                   ;;
        -d | --database )          shift
                                   blast_database=$1
                                   ;;
        -h | --help )              usage script
                                   exit
                                   ;;
        * )                        usage script
                                   exit 1
    esac
    shift
done

if [ -z "$contigs" ]
then
        usage $script
        exit
fi

mkdir -p $output_folder

timestamp=`date "+%m-%d_%T"`
blastfile="${output_folder}/megablast_output_${timestamp}.tmp"
echo "Outputting blast results to file ${blastfile}"
megablast -a 16 -d $blast_database -i $contigs -o $blastfile
python "${path_to_script}/process_blast_output.py" $blastfile 30 100 > "${output_folder}/blast_results_${timestamp}.txt"
echo "Processed blast results are in file ${output_folder}/blast_results_${timestamp}.txt"

