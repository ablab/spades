#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 infile outfile\n";
  exit(0);
}

$infile = shift @ARGV;
$outfile = shift @ARGV;

open (IN, $infile) or die "cannot open $infile\n";
open (OUT, ">$outfile") or die "cannot open $outfile\n";


<IN>;
while (<IN>) {
  $line = $_;
  $line =~ /\d+\s+(\d+)\s+\d+\s+(\d+)\s+\d+\s+\d+\s+\d/;
  print OUT "$1 $2\n";
}
