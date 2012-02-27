#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 bed file unmapped file\n";
}
open(DATA, $ARGV[0]) or die "cannot open $ARGV[0]\n";
open(UNM, $ARGV[1]) or die "cannot open $ARGV[1]\n";

@unmasked = ();
while(<UNM>) {
  $line = $_;
  chomp($line);
  if ($line !~ /^\#/) {
    @vals = split(/\s+/, $line);
    $line = join("", @vals);
#    print "line: $line\n";
    push @unmasked, $line;
  }
}

while (<DATA>) {
  $line = $_;
  chomp $line;
  @vals = split(/\s+/, $line);
  $searchkey = join("",@vals);
  @search = grep(/$searchkey/,@unmasked);
  if ($#search < 0) {
    print "$line\n";
  }
}
