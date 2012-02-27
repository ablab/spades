#!/usr/bin/env perl

$in = shift @ARGV;
open (IN, "$in") or die "cannot open $in\n";
@species = ();
@lengths = ();
while (<IN>) {
  $line = $_;
  if ($line =~ /bin: (\d+)/) {
  #  print "bin: $1\n";
   $bin = $1; 
  }
  else {
    $line =~ /\s*(\S+)\s+(\-?\d+)\s+(\-?\d+)/;
    $spec = $1;
    $diff = $3 - $2;
    push @species, $spec;
    push @lengths, $diff;
  }
}




