#!/usr/bin/env perl

$infile = shift @ARGV;
$outfile = shift @ARGV;

open(IN, "$infile") or die "cannot open $infile \n";
open(OUT, ">$outfile") or die "cannot open $outfile \n";


while (<IN>) {
  $line = $_;
  if ($line !~ /\>/) {
    $line = uc $line;
  }
  print OUT $line;
}

close(IN);
close(OUT);
