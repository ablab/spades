#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Matches reads by bowtie to genomes / sequences and throws out the reads that match.

function usage ()
{
    echo "usage: $1 ( -a <single_file_with_all_reads.fastq[.gz] | -1 <first_reads.fastq[.gz]> -2 <second_reads.fastq[.gz]> [-u1 <first_unpaired.fastq[.gz]> -u2 <second_unpaired.fastq[.gz]>] ) -i <genome_index[,genome_index2,...]> [-o <output folder>] [-d <blast database>]"
}

function copy_if_not_same ()
{
    if [ ! $1 -ef $2 ]
    then
         cp $1 $2
    fi
}

function link_if_not_same ()
{
    if [ ! $1 -ef $2 ]
    then
         ln -s $1 $2
    fi
}

function clean_unpaired ()
{
    echo "processing unpaired reads"
    reads=$1
    indices=$2
    output_folder=$3
    filename=$(basename $reads)
    extension=${filename##*.}
    targetfile="${output_folder}/${filename}"
    if [ "$extension" = "gz" ]
    then
      copy_if_not_same $reads $targetfile
      gunzip $targetfile
      filename=${filename%.*}
      extension=${filename##*.}
    else
      link_if_not_same $reads $targetfile
    fi
    myreads="${output_folder}/${filename}"
    filename=${filename%.*}
    cleanedreads="${output_folder}/${filename}_cleaned.${extension}"
    mv $myreads $cleanedreads
    readstmp=$cleanedreads.tmp
    for i in "${indices[@]}"
    do
        echo "running bowtie on index $i"
	bowtie -c -q -p 16 --suppress 6,7,8 $i $cleanedreads --un $readstmp > /dev/null 2> /dev/null
        mv $readstmp $cleanedreads
    done
}

function clean_paired ()
{
    echo "processing paired reads"
    left=$1
    right=$2
    indices=$3
    output_folder=$4
    path_to_script=$5
    filename1=$(basename $left)
    extension1=${filename1##*.}
    if [ "$extension1" = "gz" ]
    then
      copy_if_not_same $left "${output_folder}/${filename1}"
      gunzip "${output_folder}/${filename1}"
      filename1=${filename1%.*}
      extension1=${filename1##*.}
    else
      link_if_not_same $left "${output_folder}/${filename1}"
    fi
    myleft="${output_folder}/${filename1}"
    filename1=${filename1%.*}

    filename2=$(basename $right)
    extension2=${filename2##*.}
    if [ "$extension2" = "gz" ]
    then
      copy_if_not_same $right "${output_folder}/${filename2}"
      gunzip "${output_folder}/${filename2}"
      filename2=${filename2%.*}
      extension2=${filename2##*.}
    else
      link_if_not_same $right "${output_folder}/${filename2}"
    fi
    myright="${output_folder}/${filename2}"
    filename2=${filename2%.*}
    cleanedleft="${output_folder}/${filename1}_cleaned.${extension1}"
    filteredleft="${output_folder}/${filename1}_cleaned_filtered.${extension1}"
    cleanedright="${output_folder}/${filename2}_cleaned.${extension2}"
    filteredright="${output_folder}/${filename2}_cleaned_filtered.${extension2}"
    mv $myleft $cleanedleft
    mv $myright $cleanedright
    timestamp=`date "+%m-%d_%T"`
    for i in "${indices[@]}"
    do
        echo "running bowtie on index $i"
        #bowtie -c -q -p 16 --suppress 6,7,8 --minins 0 --maxins 10000 --fr --tryhard $i -1 $cleanedleft -2 $cleanedright > paired_log.tmp 2> /dev/null
        bowtie -c -q -p 16 --suppress 6,7,8 $i $cleanedleft > "${output_folder}/paired_log1_${timestamp}.tmp" 2> /dev/null
        bowtie -c -q -p 16 --suppress 6,7,8 $i $cleanedright > "${output_folder}/paired_log2_${timestamp}.tmp" 2> /dev/null
        echo "running filter_reads.py"
        python ${path_to_script}/../reads_utils/filter/filter_reads.py $cleanedleft $cleanedright "${output_folder}/paired_log1_${timestamp}.tmp" "${output_folder}/paired_log2_${timestamp}.tmp"
        mv $filteredleft $cleanedleft
        mv $filteredright $cleanedright
    done
    rm "${output_folder}/paired_log1_${timestamp}.tmp" "${output_folder}/paired_log2_${timestamp}.tmp"

}

##### Main

left=
right=
left_unpaired=
right_unpaired=
indices=()
indexnames=()
all_reads=
script=$0
output_folder="."
blast_database="/acestorage/blast/db/nt"

scriptname=${script##*/}
path_to_script=${script%$scriptname}

while [ "$1" != "" ]; do
    case $1 in
        -a | --all_reads )         shift
                                   all_reads=$1
                                   ;;
        -1 | --first )             shift
                                   left=$1
                                   ;;
        -2 | --second )            shift
                                   right=$1
                                   ;;
        -u1 | --first_unpaired )   shift 
                                   left_unpaired=$1
                                   ;;
        -u2 | --second_unpaired )  shift
                                   right_unpaired=$1
                                   ;;
	-o | --output )            shift
                                   output_folder=$1
                                   ;;  
        -d | --database )          shift
                                   blast_database=$1
                                   ;;                      
        -i | --index )             shift
                                   indexnames=$1
                                   ;;
        -h | --help )              usage $script
                                   exit
                                   ;;
        * )                        usage $script
                                   exit 1
    esac
    shift
done

if [ -z "$left" -o -z "$right" ]
then
        if [ -z "$all_reads" ]
        then
        	usage $script
		exit
	fi
fi

if [ ${#indexnames[@]} -eq 0 ]
then
	usage $script
	exit
fi

dbname=${blast_database##*/}
path_to_db=${blast_database%$dbname}

mkdir -p $output_folder

old_IFS=${IFS}
IFS=","
for i in $indexnames #"${indexnames[@]}"
do
	indices[$[${#indices[@]}]]="${path_to_db}/index/${i:0:2}/${i}"
	if [ ! -e "${path_to_db}/index/${i:0:2}/${i}.1.ebwt" ]
	then
		echo "building index for sequence $i"
		mkdir -p "${path_to_db}/index/${i:0:2}"
		python ${path_to_script}/read_fasta_using_index.py $i "${blast_database}.fasta" > "${i}.fasta.tmp"
		bowtie-build "${i}.fasta.tmp" "${path_to_db}/index/${i:0:2}/${i}"
		rm "${i}.fasta.tmp"
	fi
done
IFS=${old_IFS}

if [ "$all_reads" ]
then
    echo "Splitting reads file into multiple files."
    filename=$(basename $all_reads)
    extension=${filename##*.}
    if [ "$extension" = "gz" ]
    then
      cp $all_reads "${output_folder}/${filename}"
      gunzip "${output_folder}/${filename}"
      filename=${filename%.*}
      extension=${filename##*.}
    else
      ln -s $all_reads "${output_folder}/${filename}"
    fi
    python ${path_to_script}/../reads_utils/conversion/split_any_fastq.py "${output_folder}/${filename}"
    filename=${filename%.*}
    left="${output_folder}/${filename}_left.${extension}"
    right="${output_folder}/${filename}_right.${extension}"
    left_unpaired="${output_folder}/${filename}_single_left.${extension}"
    right_unpaired="${output_folder}/${filename}_single_right.${extension}"
fi

clean_paired $left $right $indices $output_folder $path_to_script

if [ "$left_unpaired" ]
then
        clean_unpaired $left_unpaired $indices $output_folder
fi

if [ "$right_unpaired" ]
then
        clean_unpaired $right_unpaired $indices $output_folder
fi

