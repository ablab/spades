#!/usr/bin/env perl
use strict;

if ($#ARGV != 2) { 
  print "usage: $0 seq  startpos  endpos\n";
  exit(-1);
}


open(SEQ, @ARGV[0]) or die "cannot open @ARGV[0]";

my ($startpos, $endpos) = (@ARGV[1], @ARGV[2]);


my $seq = "";
my $line;
<SEQ>;
while (<SEQ>) {
  $line = $_;
  chomp($line);
  $seq = $seq . $line;
}

if ((length($seq) < $endpos ) || (length($seq) < $startpos)) { 
  printf("invalid bounds!\n");
  exit(-1);
}


my $subseq = substr($seq, $startpos, $endpos - $startpos + 1);

printf(">seq%d_%d\n", $startpos, $endpos);
while (length($subseq) > 0) { 
  $line = substr($subseq, 0, 50);
  $subseq = substr($subseq, 50, length($subseq));
  print "$line\n";
}
