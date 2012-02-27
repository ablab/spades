#!/usr/bin/env perl

use common;

if ($#ARGV != 3) {
  print "usage: $0 binFastaFile database region output\n";
  print "prints the human coordinates of each species\n";
  exit(0);
}
$binFile = shift @ARGV;
$database = shift @ARGV;
$region   = shift @ARGV;
$outName  = shift @ARGV;

print "opening $binFile\n";
open (BINFILE, $binFile) or die "cannot open binfile $binFile\n";

$arch = $ENV{"MACHTYPE"};
if (exists $ENV{"GCFG"} ) {
  $cfg = $ENV{"GCFG"};
  $configOpt = " -f $cfg ";
}

print "opening $outName\n";
open (OUT, ">$outName") or die "cannot open output $outName\n";
while (<BINFILE>) {
  $line = $_;
  if ($line=~/>(\w+)\.(\d+)\.(\d+)/) {
    $species = $1;
    $start = $2;
    $end   = $3;
    $humanStart = common::GetOrthoPosOnly($species, "human", $region, $start, 0, $database, 1);
    $humanEnd   = common::GetOrthoPosOnly($species, "human", $region, $end, 0, $database, 1);
    print OUT "$species $start $end $humanStart $humanEnd\n";
  }
}

close OUT;
close BINFILE;
