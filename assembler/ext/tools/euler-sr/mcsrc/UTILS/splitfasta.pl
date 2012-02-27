#!/usr/bin/env perl

if ($#ARGV < 2) {
  print "usage: $0 infile stride base [ext]\n";
  print "Split a fasta file infile of n reads into n/stride files\n";
	print "The number of reads in each file is stride, except the last file\n";
	print "may be shorter.\n";
	print "Files will be named base.# or base.#.ext (#=1,2,3,...)\n";
  exit(0);
}
$ext = "fasta";
$infile = shift @ARGV;
$stride = shift @ARGV;
$base   = shift @ARGV;
if ($#ARGV >= 0) {
  $ext    = shift @ARGV;
}

open(IN, $infile) or die "cannot open $infile\n";
$index = 0;
$file  = 0;
$first = 1;
while(<IN>) {
  $line = $_;
  if ($line =~ /^>/) {
    if ($index % $stride == 0) {
      if ($first == 0) {
	close OUT;
	$first = 0;
      }
      $file++;
      open(OUT, ">$base.$file.$ext");
      print "opening $base.$file.$ext\n";
    }
    $index++;
  }
  print OUT $line;
}
