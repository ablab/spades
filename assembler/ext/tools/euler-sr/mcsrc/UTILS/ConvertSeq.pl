#!/usr/bin/env perl

use Bio::SeqIO;

if (@ARGV < 2) {
  print "usage: convertSeq to_format infile [-f from_format -w width]\n";
  exit(1);
}
@infileNames = ();
$to = shift @ARGV;
$from = "";
$width = 50;
while ($#ARGV >= 0 ){ 
  $arg = shift @ARGV;
  if ($arg eq "-f") {
    $from = shift @ARGV;
  }
  if ($arg eq "-w" ) {
    $width = shift @ARGV;
  }
  else {
    push @infileNames, $arg;
  }
}

foreach $inFileName (@infileNames) {
  if ($from ne "") {
    $from = shift@ARGV;
    $seqio_obj = Bio::SeqIO->new(-file => $inFileName,
				 -format => $from);
  }
  else {
    $seqio_obj = Bio::SeqIO->new(-file => $inFileName);
  }
  $seq_obj = $seqio_obj->next_seq;
  $inFileName =~ /(.*)\.[^\.]+/;
  $inStem = $1;
  $outFileName = $inStem . "." . $to;
  $seqOut_obj = Bio::SeqIO->new(-file => ">$outFileName",
				-format => $to);
  $seqOut_obj->width($width);
  $seqOut_obj->write_seq($seq_obj);
}
