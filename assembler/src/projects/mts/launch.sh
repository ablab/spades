#!/bin/bash

function launch {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "Error with $1" >&2
        exit 239
    fi
    return $status
}

if [ "$#" -lt 3 ]; then
    echo "launch.sh <in_dir> <out_dir> <sample_cnt> [thread_cnt = 24]"
    exit 1
fi

export PATH=~/soft/:$PATH
export SPADES_PATH=~/git/algorithmic-biology/assembler/ #/home/snurk/build_30.11/SPAdes-3.7.0-dev-Linux/bin/
export SCRIPT_PATH=~/git/algorithmic-biology/assembler/src/mts/scripts/
export BIN=~/git/algorithmic-biology/assembler/build/release/bin/

export K=55
export min_contig_length=2000

export in_dir=$1
export out_dir=$2

rm -rf tmp
mkdir tmp
#rm -rf $out_dir
mkdir -p $out_dir
rm $out_dir/recount*
#touch $out_dir/recount

#saves_suff=K$K/saves/distance_estimation
export saves_suff=K$K/saves/07_before_repeat_resolution/graph_pack
#TODO add auto detection
export sample_cnt=$3

export thread_cnt=24
if [ "$#" -ge 4 ]; then
    export thread_cnt=$4
fi

export process_cnt=6

#input_dir out_prefix (i.e. $out_dir/sampleX) thread_cnt
function assemble {
    contigs=$2.fasta

    if [ -f $contigs ] && [ ! -f $out_dir/recount ] ;
    then
        echo "Assembly of $contigs already done"
    else
        #TODO longer k sequence?
        $SPADES_PATH/spades.py --meta -t $3 -1 $1/left.fastq -2 $1/right.fastq -o ${2}_assembly
        cp ${2}_assembly/scaffolds.fasta $contigs
        touch $out_dir/recount.tmp
    fi
}

#out_prefix (i.e. $out_dir/sampleX)
function count_cov_vectors {
    if [ -f $1.id ] && [ ! -f $out_dir/recount ] ;
    then
        echo "Abundances for contigs $contigs already counted"
    else
        $BIN/contig_abundance_counter $K ${1}_assembly/$saves_suff $1.fasta $sample_cnt $out_dir/kmers $1 $min_contig_length
        touch $out_dir/recount.tmp
    fi
}

#out_prefix (i.e. $out_dir/sampleX)
function propagate {
    if [ -f $1.ann.prop ] && [ ! -f $out_dir/recount ] ;
    then
        echo "Propagated annotations $1.ann.prop already counted"
    else
        $BIN/annotation_propagator $K ${1}_assembly/$saves_suff $1.fasta $1.ann $1.ann.prop
        touch $out_dir/recount.tmp
    fi
}

#sample_number
function bin_reads {
    $BIN/read_binning $K $out_dir/sample${1}_assembly/$saves_suff $out_dir/sample${1}.fasta $out_dir/sample${1}.ann.prop $in_dir/sample${1}/left.fastq $in_dir/sample${1}/right.fastq $out_dir/binning sample${1}
}

#tmp start
#printf "4\n5\n6\n8\n9" | xargs -I {} -n 1 -P 4 $BIN/annotation_propagator $K $out_dir/sample{}_assembly/$saves_suff $out_dir/sample{}.fasta $out_dir/sample{}.ann $out_dir/sample{}.ann.prop CAG0050
#seq 1 $sample_cnt | xargs -I {} -n 1 -P 10 $BIN/annotation_propagator $K $out_dir/sample{}_assembly/$saves_suff $out_dir/sample{}.fasta $out_dir/sample{}.ann $out_dir/sample{}.ann.prop

#printf "4\n5\n6\n8\n9" | xargs -I {} -n 1 -P 5 $BIN/read_binning $K $out_dir/sample{}_assembly/$saves_suff $out_dir/sample{}.fasta $out_dir/sample{}.ann.prop $in_dir/sample{}/left.fastq $in_dir/sample{}/right.fastq $out_dir/binning sample{} CAG0050

#seq 1 $sample_cnt | xargs -I {} -n 1 -P 10 $BIN/read_binning $K $out_dir/sample{}_assembly/$saves_suff $out_dir/sample{}.fasta $out_dir/sample{}.ann.prop $in_dir/sample{}/left.fastq $in_dir/sample{}/right.fastq $out_dir/binning sample{}
#tmp end


echo "Generating sample descriptions"
$SCRIPT_PATH/dataset_desc_gen.sh $in_dir/sample $out_dir

#descs=$(ls $out_dir/*.desc)
descs=""
for i in `seq 1 $sample_cnt` ;
do
   descs="$descs $out_dir/sample$i.desc"
done

#TODO add md5 check
if [ -f $out_dir/kmers.mpl ] && [ ! -f $out_dir/recount ] ;
then
    echo "Kmer multiplicities already counted"
else
    #TODO set number of threads for kmc2
    #TODO make kmc2 launch explicit
    echo "Gathering $(($K + 1))-mer multiplicities from $descs"
    echo "$BIN/kmer_multiplicity_counter -k $(($K + 1)) -o $out_dir/kmers --mult 2 --sample 3 $descs"
    $BIN/kmer_multiplicity_counter -k $(($K + 1)) -o $out_dir/kmers --mult 2 --sample 3 $descs
    touch $out_dir/recount
fi

export -f assemble
seq 1 $sample_cnt | xargs -I {} -n 1 -P $process_cnt bash -c 'assemble "$@"' _ $in_dir/sample{} $out_dir/sample{} $((thread_cnt / process_cnt))
mv $out_dir/recount{.tmp,} 2>/dev/null

export -f count_cov_vectors
seq 1 $sample_cnt | xargs -I {} -n 1 -P $process_cnt bash -c 'count_cov_vectors "$@"' _ $out_dir/sample{}
mv $out_dir/recount{.tmp,} 2>/dev/null

canopy_input=$out_dir/canopy.in
if [ -f $canopy_input ] && [ ! -f $out_dir/recount ] ;
then
    echo "Canopy input already prepared"
else
    echo "Preparing canopy input"
    abundance_fns=$(ls $out_dir/*.id)
    $SCRIPT_PATH/make_canopy_input.py $canopy_input $abundance_fns
    touch $out_dir/recount
fi

canopy_output=$out_dir/canopy.out
canopy_profile=$out_dir/canopy.prof

#FIXME add canopy launch
#cp $in_dir/../canopy.out $canopy_output
if [ -f $canopy_output ] && [ ! -f $out_dir/recount ] ;
then
    echo "Canopy clustering done already"
else
    $SCRIPT_PATH/canopy_launch.sh $canopy_input $canopy_output $canopy_profile
    touch $out_dir/recount 
fi

if [ -f $out_dir/sample1.ann ] && [ ! -f $out_dir/recount ] ;
then
    echo "Raw annotations already prepared"
else
    echo "Preparing raw annotations"
    $SCRIPT_PATH/parse_canopy_out.py $canopy_output $out_dir
    touch $out_dir/recount 
fi

export -f propagate
seq 1 $sample_cnt | xargs -I {} -n 1 -P $process_cnt bash -c 'propagate "$@"' _ $out_dir/sample{}
mv $out_dir/recount{.tmp,} 2>/dev/null

if [ -d $out_dir/binning ] && [ ! -f $out_dir/recount ] ;
then
    echo "Read binning already done"
else
    echo "Binning reads"

    export -f bin_reads
    seq 1 $sample_cnt | xargs -I {} -n 1 -P $process_cnt bash -c 'bin_reads "$@"' _ {}
    touch $out_dir/recount 
fi

rm -rf tmp
rm -rf $out_dir/recount*
