#!/usr/bin/env perl
if ($#ARGV < 0) {
  print "usage: $0 fasta_file\n";
  print "output: #actg #ACTG #unknown total\n"; 
  exit(0);
}
$infile = shift @ARGV;

open (IN, "$infile\n");

$low = 0;
$high = 0;
$unk =0;
$na = 0;
$nc = 0;
$nt = 0;
$ng = 0;
while(<IN>) {
  $line = $_;
  if ($line !~ /^>/) { 
    $lineCopy = $line;
    $l = tr/actg//;
    $h  =tr/ACTG//;
    $n  = tr/Nn//;
    $na += ($lineCopy =~ tr/aA//);
    $nc += ($lineCopy =~ tr/cC//);
    $ng += ($lineCopy =~ tr/gG//);
    $nt += ($lineCopy =~ tr/tT//);
    $low = $low + $l;
    $high = $high + $h;
    $unk= $unk + $n;
  }
}
$total = $low + $high + $unk;

print "$low $high $unk $total $na $nt $nc $ng\n";
