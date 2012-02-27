#!/usr/bin/env perl

$textfile = shift @ARGV;
$dbname = shift @ARGV;
$region = shift @ARGV;
use common;


$specFile = shift @ARGV;

@species = common::ReadSpeciesFile($specFile);


foreach $s1 (@species) {
  foreach $s2 (@species) {
    $chainFile = "bzchain/$s1.ENm001.fa.$s2.ENm001.fa.chain";
    system("~/projects/mcsrc/inversions/storeChainTables.pl $textfile $dbname $s1 $s2 $region");
  }
}
