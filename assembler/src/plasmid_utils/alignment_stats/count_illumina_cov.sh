#!/bin/bash
if [ "$#" -lt 4 ]; then
    echo "script.sh <fragments> <left_ill> <right_ill> <out_dir> [thread count = 8]"
    exit 1
fi

thread_cnt=8
if [ "$#" -ge 5 ]; then
    thread_cnt=$5
fi

fragments=$1
out_dir=$4
exec_path=$(dirname $(readlink -e $0))

if [ ! -e $exec_path/illumina_cov.py ] ; then
    exit 239
fi

mkdir -p $out_dir
if [ ! -e $out_dir/alignment.bam ] ; then
    echo "Constructing index"
    bwa index $fragments
    bowtie2-build $fragments $out_dir/index
    #bowtie2 -p $thread_cnt --met-stderr --no-unal --threads $thread_cnt --no-discordant --no-mixed --fast -a -x $out_dir/index -1 $2 -2 $3 2> $out_dir/bowtie.log | samtools view -Sbh - > $out_dir/alignment.bam
    echo "Aligning reads"
    bowtie2 --met-file $out_dir/align.log --no-unal -X 1000 --threads $thread_cnt -x $out_dir/index -1 $2 -2 $3 2> $out_dir/bowtie.log | samtools view -Sbh - > $out_dir/alignment.bam
    #bwa mem -x ont2d $fragments -t $thread_cnt $2 2> $out_dir/bwa.log | samtools view -Sbh - > $out_dir/alignment.bam
else
    echo "Reads already aligned"
fi

echo "Counting average coverage"
$exec_path/illumina_cov.py $fragments $out_dir/alignment.bam $out_dir/result.txt
#echo "Counting binned coverage"
#$SCRIPTS_PATH/alignment_stats/binned_coverage.py $fragments $out_dir/alignment.bam 1000 > $out_dir/result.txt
