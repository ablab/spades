#!/usr/bin/env perl
use common;
use inversions;

if ($#ARGV < 1) {
  print "usage: $0 invFile outbase\n";
  print "infile has a list species pairs, and each pair is followed by a \n";
  print "list of coordinates of each inversion found between the pair of\n";
  print "species.\n";
  print "outbase is the base file name (usually a directory) that each overlapping\n";
  print "inversion is printed to\n";
  exit(0);
}

$invFile = shift @ARGV;
$outFile = shift @ARGV;


# Grab the sequences
@refSpecies;
my %invPairs;
inversions::ReadInversionList($invFile, \%invPairs);

my @specList;
inversions::GetSpeciesList(\%invPairs, \@specList);

# initialize the list of species sites
my %invSites;
foreach $spec (@specList) {
  $invSites{$spec} = ();
}

inversions::CreateVertices(\%invPairs, \@specList, \%invVertices);

inversions::InitOrthInvGraph(\%invVertices, \@specList, \%invEdges);

inversions::CreateEdges(\%invPairs, \@specList, \%invVertices, \%invEdges);

inversions::PrintConnectedComponents(\%invVertices, \%invEdges, \@specList, $outFile);




