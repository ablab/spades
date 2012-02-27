#!/usr/bin/env perl
open(IN, $ARGV[0]) or die "must supply an input file\n";

use POSIX;

while(<IN>) {  
  $line = $_;
  $line =~ tr/:\-/  /;
#  $line =~ /tr/\-/ /g;
  @vals = split(/\s+/, $line);
  $len = $vals[2] - $vals[1];
  $offset = POSIX::floor($len / 2);
  $start = $vals[1] - $offset;
  $end   = $vals[2] + $offset;
  print "$vals[0] $start $end\n";
}
