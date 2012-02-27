#!/usr/bin/env perl

use grimm;
# parse grimm output to find which sequences are inverted without
# any overlapping rearrangements
#
# The current sample input is:
#
#Running:                        /Users/mchaisso/projects/software/GRIMM-2.01/grimm
#Input:                          mouse.human/mgr_macro.txt
#Genome:                         Multichromosomal
#Signs:                          Signed
#Number of Genomes:              2
#Number of Genes:                15 + 2 caps
#Number of Chromosomes:          1 (multichromosomal)
#Multichromosomal Distance:      7
# 
#======================================================================
# 
#Optimal capping:
#>genome1
#16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17
#>genome2
#16 1 -2 3 -4 5 -6 7 -8 9 -10 11 -12 13 -14 15 17
# 
#======================================================================
# 
#An optimal sequence of rearrangements:
#Step 0: (Source)
#  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15   $
#Step 1: Chrom. 1, gene 2 [2] through chrom. 1, gene 2 [2]: Reversal
#  1  -2   3   4   5   6   7   8   9  10  11  12  13  14  15   $
#Step 2: Chrom. 1, gene 4 [4] through chrom. 1, gene 4 [4]: Reversal
#  1  -2   3  -4   5   6   7   8   9  10  11  12  13  14  15   $
#Step 3: Chrom. 1, gene 6 [6] through chrom. 1, gene 6 [6]: Reversal
#  1  -2   3  -4   5  -6   7   8   9  10  11  12  13  14  15   $
#Step 4: Chrom. 1, gene 8 [8] through chrom. 1, gene 8 [8]: Reversal
#  1  -2   3  -4   5  -6   7  -8   9  10  11  12  13  14  15   $
#Step 5: Chrom. 1, gene 10 [10] through chrom. 1, gene 10 [10]: Reversal
#  1  -2   3  -4   5  -6   7  -8   9 -10  11  12  13  14  15   $
#Step 6: Chrom. 1, gene 12 [12] through chrom. 1, gene 12 [12]: Reversal
#  1  -2   3  -4   5  -6   7  -8   9 -10  11 -12  13  14  15   $
#Step 7: Chrom. 1, gene 14 [14] through chrom. 1, gene 14 [14]: Reversal (Destination)
#  1  -2   3  -4   5  -6   7  -8   9 -10  11 -12  13 -14  15   $

if ($#ARGV != 1) {
  print "usage: FindGrimmInversions.pl mgr_macro.txt blocks.txt\n";
  exit(0);
}

$grimmOut = shift @ARGV;
$blocksFile = shift @ARGV;
open(GO, $grimmOut) or die "cannot open $grimmOut\n";
$numGenes = 0;
while(<GO>) {
  if (/Number of Genes:\s+(\d+) .*/) {$numGenes = $1; last;}
}
push @identity, 0;
for ($i = 0; $i < $numGenes; $i++) {
  push @identity, 0;
}
while(<GO>) {
  if (/Multichromosomal Distance:\s+(\d+)/) {$numRearrange = $1; last;}
}
while(<GO>) { 
  if(/An optimal sequence of rearrangements/) {last;}
}
# discard the identity permutation
<GO>;
<GO>;
for ($r = 0; $r < $numRearrange; $r++) {
  $line = <GO>;
  $line =~ /gene (\d+).*gene (\d+)/;
  $start = $1; $end = $2;
  for ($i = $start; $i <= $end; $i++) {
    $identity[$i]+= ($end - $start + 1);
  }
  <GO>; # throw away the resulting permutation
}


# get the blocks
@blocks = ();
grimm::ReadSyntenyBlocks($blocksFile, \@blocks);
for ($i = 1; $i <= $numGenes; $i++) {
  if ($identity[$i] == 1) {
    $rb = $blocks[$i]{"start1"};
    $re = $blocks[$i]{"len1"} + $rb - 1;
    $qb = $blocks[$i]{"start2"};
    $qe = $blocks[$i]{"len2"} + $qb - 1;
    print "$rb $re $qb $qe\n";
  }
}


