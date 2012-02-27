#!/usr/bin/env perl
$in = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
while(<IN>) {
  $line = $_;
  if ($line =~ /^>/) {
				print $line;
  }
  else { $pl = 0;
    $line =~ tr/N/A/;
    print $line;
 }
}
 
