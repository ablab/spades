#!/usr/bin/env perl

use Bio::SeqIO;

if ($#ARGV != 1) {
  print "usage: $0 in out \n";
  exit(0);
}


$in = shift @ARGV;
$width = shift @ARGV;

$seqIn = new Bio::SeqIO(-file=>$in);
$seqOut = new Bio::SeqIO(-fh=>\*STDOUT, -format=>'fasta');
$seqOut->width($width);
while ($seq = $seqIn->next_seq) {
  $seqOut->write_seq($seq);
}
