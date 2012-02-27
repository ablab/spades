#!/usr/bin/env perl
use strict;

if ($#ARGV ne 2) {
  print "usage: filtermatepairs discardsFile sourceMates destMates\n";
  exit(-1);
}


my $discards = shift @ARGV;
my $sourceMates = shift @ARGV;
my $destMates  = shift @ARGV;


my %discardedReads = ();
# somehow get the discarded reads into a hash

open (DISCARDS, "$discards") or die "cannot open $discards\n";
while (<DISCARDS>) {
  /\S+\D(\d+)/;
  $discardedReads{$1}  = 1;
}
close(DISCARDS);

open (SOURCE, "$sourceMates") or die "cannot open $sourceMates";
open (DEST,   ">$destMates") or die "cannot write to $destMates\n";
while (<SOURCE>) {
  /\S+ (\d+) (\d+)/;
  if (not (exists($discardedReads{$1}) || exists($discardedReads{$2}))) {
    print DEST $_;
  }
}

close(SOURCE);
close(DEST);
