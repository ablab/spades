#!/bin/sh

ARTSIMPATH=/home/tdvorkina/soft/art_bin_MountRainier/
REF=ref/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa

$ARTSIMPATH/art_illumina -ss HS25 -sam -i $REF -p -l 150 -f 20 -m 200 -s 10 -rs 115249 -o paired_dat
