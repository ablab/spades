#!/usr/bin/env perl
use strict;

if ($#ARGV < 0) {
  printf("usage: calcerrors.pl infile\n");
}
my $infile = @ARGV[0];


open(IN, "$infile") or die "cannot open $infile\n";
my $score = 0;
while (<IN>) {
  if (/.*got score: (\d+)*/) {
    print "found error at $_ with score: $1\n";
    $score = $score + int($1);
  }
}

print "net errors: $score\n";

