#!/usr/bin/env perl

$value = 0;

%nucs = ('a', 0, 'c', 1, 'g', 2, 't', 3, 'A', 0, 'C', 1, 'G', 2, 'T', 3);

$inStr = shift @ARGV;

@chars = split(//, $inStr);
while ($#chars >= 0) {
  $char = shift @chars;
  if (exists $nucs{$char}) {
    $value <<= 2;
    $value += $nucs{$char};
  }
  else {
    print "invalid nucleotide: $char\n";
    exit 0;
  }
}
print "$value\n";
