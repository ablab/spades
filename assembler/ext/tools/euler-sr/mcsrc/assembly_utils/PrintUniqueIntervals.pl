#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: $0 intervalFile length\n";
		exit(1);
}
open(IN, "$ARGV[0]") or die "cannot open $ARGV[0]\n";
$min = $ARGV[1];
while(<IN>) {
  $line = $_;
  if ($line =~/EDGE.*Multiplicity (\d+)/) {
    $mult = $1;
    if ($mult == 1) { 
      $printIntv = 1;
    }
    else { 
      $printIntv = 0;
    }
  }
  else {
     if ($printIntv) {
       $line =~ /INTV (\d+) (\d+) (\d+) (\d+)/;
       $strand = $1;
       $start   = $2;
       $length = $3;
       $end = $start + $length;
       if ($strand == 0 && $length > $min) {
         print "$start $end $length\n";
       }
     }
  }
}
