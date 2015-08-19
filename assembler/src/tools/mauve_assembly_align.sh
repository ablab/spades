#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################



if [ $# != 3 ] 
then
    echo "Usage: 'mauve_dir' 'contigs_dir' 'reference_file'"
    exit
fi

mauve_dir=$1
contigs_dir=$2
reference=$3

export CLASSPATH=$CLASSPATH:$mauve_dir/ext
export PATH=$PATH:$mauve_dir

ref_name=$(basename $reference)
output_dir=$contigs_dir/mauve

#rm -rf $output_dir
mkdir $output_dir

for contigs in $contigs_dir/*.fasta
do
    name=$(basename $contigs)
    echo "Processing $name into $output_dir/${name}"
    java -Xmx500m -cp $mauve_dir/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $output_dir/${name} -ref $reference -draft $contigs
    mv $output_dir/${name}/alignment2/${ref_name} $output_dir/${name}/alignment2/${ref_name}2
done

#./mauveAligner ~/Dropbox/lab/assembly_analysis/bcereus/50x/mauve_out/ray_ordered/alignment2/ray.K35.scf.fasta ~/Dropbox/lab/assembly_analysis/bcereus/50x/mauve_out/ray_ordered/alignment2/ray.K35.scf.fasta.sml ~/Dropbox/lab/assembly_analysis/bcereus/50x/mauve_out/spades3_ordered/alignment3/spades3_0.fasta ~/Dropbox/lab/assembly_analysis/bcereus/50x/mauve_out/spades3_ordered/alignment3/spades3_0.fasta.sml ~/Dropbox/lab/assembly_analysis/bcereus/sequence_full.fasta ~/Dropbox/lab/assembly_analysis/bcereus/sequence_full.fasta.sml --output ~/Dropbox/lab/assembly_analysis/bcereus/50x/mauve_out/alignment

ordered_with_sml=$(ls $output_dir/*/alignment2/*.fasta | xargs -I {} echo {} {}.sml | paste -s -d ' ')

#echo $ordered_with_sml

$mauve_dir/mauveAligner $reference $output_dir/${ref_name}.sml $ordered_with_sml --output $contigs_dir/alignment

rm -rf $output_dir
