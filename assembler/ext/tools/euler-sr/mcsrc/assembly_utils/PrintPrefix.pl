#!/usr/bin/env perl

$in = shift @ARGV;
$length = shift @ARGV;
open(IN, $in) or die "cannot open $in\n";
$start = 0;
while(<IN>) {
  if ($_ =~ /^>/) {
    print $_;
    $start = 1;
  }
  else {
   if ($start == 1) {
     $read = substr $_, 0, $length;
     print "$read\n";
   }
   $start = 0;
  }
}
