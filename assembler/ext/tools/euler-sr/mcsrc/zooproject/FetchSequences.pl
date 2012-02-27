#!/usr/bin/env perl
use Bio::DB::GenBank;

if ($#ARGV < 3) {
  print "usage: FetchSequences sequencename target accid [accid2 ...]\n";
  exit(0);
}

$species = shift @ARGV;
$target  = shift @ARGV;
@accids = ();
while ($#ARGV >= 0) {
  push @accids, shift @ARGV;
}

$db_obj = new Bio::DB::GenBank;
$numOk = 0;
