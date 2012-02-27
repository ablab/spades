#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: $0 $infile $outfile \n";
  print "infile is a locations file (start end *)";
  print "outfile has one line for each in infile, with start'=end-1000\n";
  print "and end'=end+1000\n";
}
$infile= shift @ARGV;
$outfile = shift @ARGV;
$dist = 1000;
if ($#ARGV >= 0) {
  $dist = shift @ARGV;
}
open(IN, "$infile") or die "Cannot open $infile\n";
open(OUT, ">$outfile") or die "Cannot open $outfile\n";


while (<IN>) {
  @vals = split;
  $v1 = @vals[1] - $dist;
  $v2 = @vals[1] + $dist;
  if ($v1 > 0 && $v2 > 0) {
    print OUT "$v1 $v2\n";
  }
}
close OUT;

  
