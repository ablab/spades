#!/usr/bin/env perl
if ($#ARGV < 0) {
  print "usage: $0 fasta_file\n";
  print "output: #actg #ACTG #unknown total\n"; 
  exit(0);
}
$infile = shift @ARGV;

open (IN, "$infile\n");

$low = 0;


$gc = 0;
$gt = 0;
$ct = 0;
while(<IN>) {
  $line = $_;
  if ($line !~ /^>/) { 
			$t = tr/Tt//;
			$a = tr/Aa//;
    $g  =tr/Gg//;
    $c  = tr/Cc//;
		$gt += $g;
		$ct += $c;
		$gc += $g+$c;
			$at += $a + $t;
  }
}
if ($gc + $at > 0) {
		$pctGC = $gc / ($gc + $at);
		print "C: $ct G: $gt GC: $gc  $pctGC\n";
}

