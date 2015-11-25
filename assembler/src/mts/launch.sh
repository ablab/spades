#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "script.sh <in_dir> <out_dir>"
    exit
fi

SPADES_PATH=~/git/algorithmic-biology/assembler/
SCRIPT_PATH=~/git/algorithmic-biology/assembler/src/mts/scripts/
BIN=~/git/algorithmic-biology/assembler/build/release/bin/

K=55
thread_cnt=4

in_dir=$1
out_dir=$2

#rm -rf $out_dir
#mkdir -p $out_dir

export PATH=~/soft/:$PATH
saves_suff=K$K/saves/distance_estimation
#TODO add auto detection
sample_cnt=5
min_contig_length=1000

echo "Generating sample descriptions"
$SCRIPT_PATH/dataset_desc_gen.sh $in_dir/sample $out_dir

descs=$(ls $out_dir/*.desc)

#TODO add md5 check
if [ -f $out_dir/kmers.mpl ];
then
    echo "Kmer multiplicities already counted"
else
    #TODO set number of threads for kmc2
    #TODO make kmc2 launch explicit
    echo "Gathering $(($K + 1))-mer multiplicities from $descs"
    $BIN/kmer_multiplicity_counter -k $(($K + 1)) -o $out_dir/kmers --mult 3 --sample 2 $descs
fi

for s_dir in $in_dir/sample* ; do
    s_name=$(basename $s_dir)    
    assembly_dir=$out_dir/${s_name}_assembly
    contigs=$out_dir/$s_name.fasta

    if [ -f $contigs ];
    then
        echo "Assembly of sample $s_name already done"
    else
#TODO add meta mode
#TODO return BH?
#TODO longer k sequence?
#TODO saves instead of debug
        $SPADES_PATH/spades.py --debug -t $thread_cnt -k $K --only-assembler -1 $s_dir/left.fastq -2 $s_dir/right.fastq -o $assembly_dir

        cp $assembly_dir/scaffolds.fasta $contigs
    fi

    if [ -f $out_dir/${s_name}_abundance.mpl ];
    then
        echo "Abundances for contigs $contigs already counted"
    else
        $BIN/contig_abundance_counter $K $assembly_dir/$saves_suff $out_dir/$s_name.fasta $sample_cnt $out_dir/kmers $out_dir/${s_name}_abundance $min_contig_length
    fi
done

canopy_input=$out_dir/canopy.in
if [ -f $canopy_input ];
then
    echo "Canopy input already prepared"
else
    echo "Preparing canopy input"
    abundance_fns=$(ls $out_dir/*.id)
    $SCRIPT_PATH/make_canopy_input.py $canopy_input $abundance_fns
fi
