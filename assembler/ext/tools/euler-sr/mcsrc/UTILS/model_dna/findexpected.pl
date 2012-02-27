#!/usr/bin/env perl
use POSIX;

if ($#ARGV < 0) {
  print "usage: $0 infile probabilities [order]";
}

$inputFile = shift @ARGV;
$probFile = shift @ARGV;

$order = 1;
if ($#ARGV >= 0) {
  $order = shift @ARGV;
}




%nucs = {};

$nucs{'G'} = 0;
$nucs{'g'} = 0;
$nucs{'A'} = 1;
$nucs{'a'} = 1;
$nucs{'C'} = 2;
$nucs{'c'} = 2;
$nucs{'T'} = 3;
$nucs{'t'} = 3;

@revIndex =  ('g', 'a', 'c', 't');

open(IN, "$inputFile") or die "cannot open $inputFile\n";
$l = <IN>;
@seq = <IN>;
close IN;
chomp @seq;
$dna = join "", @seq;

open(PROB, "$probFile") or die "cannot open $probFile\n";

@probs = <PROB>;

$nprob = $#probs + 1;

$order = 0;
$tp = $nprob;
while ($tp > 1) {
  $tp = $tp / 4;
  $order++;
}


