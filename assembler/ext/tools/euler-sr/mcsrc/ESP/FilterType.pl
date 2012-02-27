#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 locations_file type. \n";
  print "filter mate-pair locations so that they are of theh same clone-type\n";
  print "locations_file is in the format: \n";
  print "start end type mate_index1 mate_index2\n";
}

$infile = shift @ARGV;
$type   = shift @ARGV;

open(IN, $infile) or die "cannot open $infile\n";

while(<IN>) {
  @vals = split(/\s+/, $_);
  if ($vals[4] eq $type) {
    print $_;
  }
}
