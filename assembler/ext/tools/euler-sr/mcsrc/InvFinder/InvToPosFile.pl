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
    if ($line =~ /(\d+) (\d+) (\d+) (\d+) (\d)/) { 
      $refStart = $1;
      $refEnd   = $2;
      $refPos   = ($refEnd + $refStart)/ 2;
      $qryStart = $3;
      $qryEnd   = $4;
      $qryPos = POSIX::floor(($qryEnd + $qryStart) / 2);
      $strand = $5;
      print "$mapTo $ref $qry $refStart $seq $strand\n";
      print "$mapTo $ref $qry $refEnd $seq $strand\n";
    }
    else {
      print "malformed coordinates line $line\n";
      exit(1);lformed
    }
  }
}
