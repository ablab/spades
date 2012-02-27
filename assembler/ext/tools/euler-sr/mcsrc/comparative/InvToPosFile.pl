#!/usr/bin/env perl
use POSIX;

if (@ARGV !=  3) {
  print "usage: InvToPosFile.pl inffile targetspecies seq\n";
  exit(0);
}

$in = shift @ARGV;
$mapTo = shift @ARGV;
$seq = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";



$ref = "";
$qry = "";
while (<IN>) { 
  $line = $_;
  if ($line =~ /EN/) {
  # line is a lav title
    if ($line =~ /(\w+)\.\w+\.fa\.(\w+)\..*/) {
      $ref = $1;
      $qry = $2;
    }
    else {
      print "malformed title line $line\n";
      exit(0);
    }
  }
  else {
    if ($line =~ /(\d+) (\d+) (\d+) (\d+).*/) {
      $refStart = $1;
      $refEnd   = $2;
      $qryStart = $3;
      $qryEnd   = $4;
      $pos = POSIX::floor(($refStart + $refEnd)/2);
      print "$mapTo $ref $qry $refStart $seq 0\n";
      print "$mapTo $ref $qry $refEnd   $seq 0\n";
    }
    else {
      print "malformed coordinates line $line\n";
      exit(1);lformed
    }
  }
}
