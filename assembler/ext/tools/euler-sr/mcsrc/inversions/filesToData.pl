#!/usr/bin/env perl

use common;

if ($#ARGV != 1) {
  print "usage: $0 invFile outFile\n";
  print "   creates a data file from inv file (tabbed-format version of the inv file\n";
  exit(0);
}

$invFileName = shift @ARGV;
$outFile     = shift @ARGV;



%pairs = ();
@refSpecies = ();

common::ReadInvFile($invFileName, \@refSpecies, \%pairs);

open (OUT, ">$outFile") or die "cannot open $outFile\n";

foreach $ref (keys %pairs) {
  foreach $qry (keys %{$pairs{$ref}}) {
    $size = scalar @{$pairs{$ref}{$qry}};
    print "$ref $qry size:$size @{$pairs{$ref}{$qry}}   \n";
    for ($i = 0; $i <= $#{$pairs{$ref}{$qry}}; $i++) {
      $refStart = $pairs{$ref}{$qry}[$i][1];
      $refEnd   = $pairs{$ref}{$qry}[$i][2];
      $refName = "$ref\_$refStart\_$refEnd\_$qry";
      $qryStart = $pairs{$ref}{$qry}[$i][4];
      $qryEnd   = $pairs{$ref}{$qry}[$i][5];
      $qryName = "$ref\_$qry\_$qryStart\_$qryEnd";
      print OUT "$refName $ref $qry $refStart $refEnd $qryName\n";
      print OUT "$qryName $qry $ref $qryStart $qryEnd  $refName\n";
    }
  }
}
