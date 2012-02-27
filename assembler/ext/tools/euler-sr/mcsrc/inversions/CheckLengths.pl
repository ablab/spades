#!/usr/bin/env perl
if ($#ARGV != 1) {
  print "usage: $0 projection_free_binning_coords length\n";
  exit(0);
}
$in = shift @ARGV;
$len = shift @ARGV;
$lineNum = 0;
open (IN, $in) or die "cannot open $in\n";
while (<IN>) { 
  $line = $_;
  $lineNum++;
  if ($line =~ /bin:\s+(\d+)/) {
    $bin = $1;
    print "bin: $bin\n";
  }
  else {
    $line =~ /\s+(\S+)\s+(\-?\d+)\s+(\-?\d+)/;
    $spec = $1;
    $start = $2; 
    $end = $3;
    
    $dist = abs($end - $start);
    
    if ($dist > $len) {
      print "$lineNum $spec $dist\n";
    }
  }
} 
