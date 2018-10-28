#!/bin/sh

SEQTKPATH=/home/tdvorkina/soft/seqtk/
IDEALALIGNERPATH=/home/tdvorkina/tmp/algorithmic-biology/assembler/build/release/bin/idealreads_aligner

REF=/Sid/tdvorkina/gralign/SPAligner_review/data/raw_data/E.coli/ref/MG1655-K12.fasta
READS=/Sid/tdvorkina/gralign/SPAligner_review/data/raw_data/E.coli/sim_pacbio/sim_pacbio_0001.fasta
GRAPH=/Sid/tdvorkina/gralign/SPAligner_review/data/raw_data/E.coli/real_graph/assembly_graph_with_scaffolds.gfa

OUT=test_reads_mapping
PREFIX=simpb
TYPE=pacbio

# mkdir $OUT
# mkdir $OUT/tmp/
python2 filterbylength.py $READS $OUT/tmp/${PREFIX}2000.fasta 2000
python2 build_idealpaths.py $OUT/tmp/${PREFIX}2000.fasta $REF $GRAPH 77 $TYPE $IDEALALIGNERPATH
# $IDEALALIGNERPATH 77 $GRAPH $OUT/tmp/refseq_${PREFIX}2000.fasta $OUT/tmp/refseq_${PREFIX}2000_mapping --spades > $OUT/tmp/refseq_${PREFIX}2000_mapping.log
cut -f1 $OUT/tmp/refseq_${PREFIX}2000_mapping.tsv > $OUT/tmp/refseq_${PREFIX}2000_mapping_names.txt
$SEQTKPATH/seqtk subseq $OUT/tmp/realseq_${PREFIX}2000.fasta $OUT/tmp/refseq_${PREFIX}2000_mapping_names.txt > $OUT/tmp/${PREFIX}2000_mapped.fasta
$SEQTKPATH/seqtk sample -s115249 $OUT/tmp/${PREFIX}2000_mapped.fasta 10000 > $OUT/input/${PREFIX}2000_mapped.fasta
cp $GRAPH $OUT/input/graph_SPAligner.gfa
cp $OUT/tmp/refseq_${PREFIX}2000_mapping.tsv $OUT/input/${PREFIX}2000_refmapping_SPAligner.tsv
# rm -rf $OUT/tmp/