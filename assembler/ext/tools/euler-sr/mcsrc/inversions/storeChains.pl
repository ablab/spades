#!/usr/bin/env perl

$dbname = shift @ARGV;
$region = shift @ARGV;
$specFile = shift @ARGV;
use common;



@species = common::ReadSpeciesFile($specFile);


foreach $s1 (@species) {
  foreach $s2 (@species) {
    $chainFile = "bzchain/$s1.fasta.$s2.fasta.chain.txt";
    system("~/projects/mcsrc/inversions/storeChainTables.pl $chainFile $dbname $s1 $s2 $region");
  }
}
