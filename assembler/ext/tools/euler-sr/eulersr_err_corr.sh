#!/bin/sh

if [ $# -ne 2 ]
then
    echo "Usage: `basename $0` INPUT_READS (in .fasta, .fa, .fastq, .fq, and .gz formats) OUTPUT_READS (in .fasta format)"
    exit
fi

if [ ! -f $1 ]
then
    echo "File with input reads ($1) doesn't exists!"
    exit     
fi

BUILD_DIR=${PWD}/../../../build/ext/euler-sr/

export EUSRC=${BUILD_DIR}/mcsrc/
export MACHTYPE=x86_64

if [ ! -d ${BUILD_DIR} ]
then 
    mkdir -p ${BUILD_DIR}
fi
if [ ! -d ${EUSRC} ]
then
    echo "copying Euler"
    cp -r ./mcsrc ${EUSRC}
fi

echo "== making Euler =="
CUR_DIR=${PWD}
cd ${EUSRC} 
make >${BUILD_DIR}/make.log 2>${BUILD_DIR}/make.err
echo "== making finished =="

input_reads=$1
output_reads=$2

input_reads_ext=${input_reads##*.}
input_reads_base=`basename "$input_reads"`
tmp_folder=`dirname "$output_reads"`
mkdir -p $tmp_folder

tmp_file=""
delete_tmp_file=0

if [ ${input_reads_ext} = "gz" ]
then
    tmp_file_unzip=${input_reads_base%.*}
    tmp_file_unzip=${tmp_folder}/${tmp_file_unzip}    
    
    echo "== unzipping input file to temporary ${tmp_file_unzip} =="
    gunzip -c ${input_reads} >${tmp_file_unzip}
    echo "== unzipping finished =="

    input_reads=$tmp_file_unzip
    input_reads_ext=${input_reads##*.}
    input_reads_base=`basename "$input_reads"`

    tmp_file=$tmp_file_unzip
    delete_tmp_file=1
fi 

if [ ${input_reads_ext} = "fastq" -o ${input_reads_ext} = "fq" ]
then
    tmp_file_fasta=${input_reads_base%.*}
    tmp_file_fasta=${tmp_folder}/${tmp_file_fasta}.fasta    

    echo "== converting fastq file to temporary ${tmp_file_fasta} with Euler Quality Trimmer =="
    ${EUSRC}/assembly/${MACHTYPE}/qualityTrimmer -fastq $input_reads -outFasta ${tmp_file_fasta}
    echo "== converting finished =="

    if [ $delete_tmp_file -eq 1 ]
    then
        echo "removing temporary file $tmp_file"
        rm $tmp_file
        delete_tmp_file=0
    fi

    input_reads=$tmp_file_fasta
    input_reads_ext=${input_reads##*.}
    input_reads_base=`basename "$input_reads"`

    tmp_file=$tmp_file_fasta
    delete_tmp_file=1
fi

if [ ${input_reads_ext} != "fasta" -a ${input_reads_ext} != "fa" ]
then    
    echo "Incorrect format of input file. Exitting.."    
    
    if [ $delete_tmp_file -eq 1 ]
    then
        echo "removing temporary file $tmp_file"
        rm $tmp_file
        delete_tmp_file=0
    fi

    exit 
fi

# error correction itself
echo "== error correction =="
${EUSRC}/assembly/Assemble.pl ${input_reads} 55 -onlyFixErrors
echo "== error correction finished =="

# clearing
if [ $delete_tmp_file -eq 1 ]
then
    echo "removing temporary file $tmp_file"
    rm $tmp_file
    delete_tmp_file=0
fi

resulting_file=${EUSRC}/fixed/$input_reads_base
if [ ! -f $resulting_file ]
then
    echo "Resulting file $resulting_file not found -- error correction failed. Exitting.."
    exit
fi

cd ${CUR_DIR}
echo "== moving temporary file $resulting_file to final destination $output_reads =="
mv $resulting_file $output_reads
echo "== moving finished =="

echo "\n\nSuccesfully finished! Corrected file is $output_reads"
