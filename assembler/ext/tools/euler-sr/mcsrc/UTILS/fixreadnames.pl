#!/usr/bin/env perl


$filename = shift @ARGV;

open (IN, "$filename") or die "cannot open $filename\n";


while(<IN>) {
  if ($_=~ ".*(>.*)") {
    print "$1\n";
  }
  else {
    print $_;
  }
}
