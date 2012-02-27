#!/usr/bin/env perl

$infile = shift @ARGV;
$spec  = shift @ARGV;
$offset = shift @ARGV;
$chrom = shift @ARGV;

open(IN, "$infile") or die "cannot open $infile\n";

print "track name=inversions description=\"inversions in ENCODE\"\n";
while (<IN>) {
  $line = $_;
  if ($line =~ /(\w+)\.\w+\.fa\.(\w+)\..*/) {
    $refSpec = $1;
    $qrySpec = $2;
  }
  else {
    if ($refSpec eq $spec) {
      if ($line =~ /(\d+) (\d+) (\d+) (\d+).*/) {
	$refStart = $1;
	$refEnd   = $2;
	$qryStart = $3;
	$qryEnd   = $4;
	$refStart += $offset;
	$refEnd += $offset;
	print "$chrom\tInvFinder\tinversion\t$refStart\t$refEnd\t.\t+\t.\t$qrySpec\n";
      }
    }
  }
}
